"""Tests for stargazer.structure: estimator bias, calibration of the surrogate tests, power on artificial
signals, time-aware lags, speed, and the behaviour ported from SeaQuest."""
from __future__ import annotations

import math
import time

import numpy as np
import pytest

from stargazer.core.rng import rng_for
from stargazer.structure import (NGramModel, benjamini_hochberg, block_counts, block_entropy, block_shuffle,
                                 compression_bits, compression_ratio, conditional_entropy, entropy,
                                 entropy_grassberger, entropy_miller_madow, entropy_plugin, entropy_rate_curve,
                                 fit_decay, gain_in_bits, grammar_gain, held_out_bits_per_symbol, kt_code_bits,
                                 lagged_mutual_information, markov_surrogate, mutual_information_decay,
                                 permutation_pvalue, repair, shuffle, tp_segmentation, tp_segmentation_null,
                                 transition_matrix)
from stargazer.structure.surrogates import doublet_shuffle

SEED = 20261122


# ---------------------------------------------------------------------------------------------------------------
# generators
# ---------------------------------------------------------------------------------------------------------------

def markov_chain(T: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    """A first-order Markov chain with transition matrix T (rows sum to 1), started from state 0."""
    cdf = np.cumsum(T, axis=1)
    u = rng.random(n)
    out = np.empty(n, dtype=np.int64)
    state = 0
    for i in range(n):
        out[i] = state
        state = int(np.searchsorted(cdf[state], u[i], side="right"))
        state = min(state, T.shape[0] - 1)
    return out


def unary_primes(limit: int = 600) -> np.ndarray:
    """The primes below `limit` as runs of ones separated by single zeros (2 -> 11, 3 -> 111, ...)."""
    primes = [p for p in range(2, limit) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    bits = []
    for p in primes:
        bits.extend([1] * p + [0])
    return np.array(bits[:-1], dtype=np.int64)


def noisy_word(rng: np.random.Generator, repeats: int = 250, flip: float = 0.05) -> np.ndarray:
    """A random 16-bit word repeated `repeats` times with a fraction `flip` of the bits flipped."""
    word = rng.integers(0, 2, 16)
    bits = np.tile(word, repeats)
    return bits ^ (rng.random(bits.size) < flip)


def max_lag_mi(seq: np.ndarray, max_lag: int = 32) -> float:
    """The statistic of the lagged-MI test: the largest plug-in MI (bits) over lags 1..max_lag. Taking the maximum
    over lags and comparing it with the same maximum on shuffles makes the test correct for the number of lags."""
    return float(np.nanmax(lagged_mutual_information(seq, 2, np.arange(1, max_lag + 1), correction="none")))


# ---------------------------------------------------------------------------------------------------------------
# estimators
# ---------------------------------------------------------------------------------------------------------------

def test_estimators_on_small_histograms():
    assert entropy_plugin([5, 5]) == pytest.approx(1.0)
    assert entropy_plugin([0, 7]) == 0.0
    assert entropy_miller_madow([5, 5]) == pytest.approx(1.0 + 1 / (20 * math.log(2)))
    # all singletons: Grassberger gives ln N - G(1) = ln N + gamma + ln 2 nats
    n = 50
    expected = (math.log(n) + 0.5772156649015329 + math.log(2)) / math.log(2)
    assert entropy_grassberger(np.ones(n, dtype=int)) == pytest.approx(expected, rel=1e-9)
    for method in ("plugin", "miller_madow", "grassberger", "chao_shen"):
        assert entropy([], method) == 0.0
        assert entropy([1000, 1000, 1000, 1000], method) == pytest.approx(2.0, abs=1e-3)
    with pytest.raises(ValueError):
        entropy([1, 2], "nsb")


def test_block_counts_and_curve_rule():
    rng = rng_for(SEED, "blocks")
    s = rng.integers(0, 3, 1000)
    for n in (1, 2, 5):
        c = block_counts(s, 3, n)
        assert c.sum() == s.size - n + 1
    assert block_counts(s, 3, 0).tolist() == [1]
    assert block_entropy(s, 3, 0) == 0.0
    # a large alphabet whose block codes would overflow int64 is still counted exactly
    big = rng.integers(0, 10 ** 6, 500)
    assert block_counts(big, 10 ** 6, 4).sum() == 497
    curve = entropy_rate_curve(s, 3, 6)
    # 3**4 = 81 <= 1000 / 5 = 200 < 3**5: orders 0..4 are reported, 5 and 6 are NaN
    assert np.all(np.isfinite(curve[:5])) and np.all(np.isnan(curve[5:]))
    assert curve[1] == pytest.approx(conditional_entropy(s, 3, 1))


# Measured mean bias (bits) over 50 i.i.d. uniform draws per cell, seed SEED, as estimate minus log2(k) for the
# symbol entropy and as the drop H(order 0) - H(order 1) (truth 0):
#   k   n        symbol entropy:  Grassberger  Miller-Madow   drop: plug-in  Grassberger  Miller-Madow
#   2   200                        +0.0006      -0.0012              +0.0034   +0.0004      -0.0002
#   2   2000/20000                  |bias| < 0.0005 for every estimator
#   4   200                        +0.0019      -0.0012              +0.0343   +0.0056      +0.0016
#   4   2000                       +0.0004      -0.0001              +0.0032   +0.0005      -0.0001
#   26  200                        +0.0034      -0.0054              +1.8528   +0.6826      +1.4124
#   26  2000                       +0.0004      +0.0003              +0.2508   +0.0036      +0.0380
#   26  20000                      +0.0001      +0.0000              +0.0225   +0.0001      -0.0000
# (measured with rng_for(1, "bias", k, n); the asserted values below use SEED and differ in the last digits.)
# At k = 26, n = 200 the 676 possible pairs are seen by only 199 overlapping pairs: Grassberger removes about
# 63 % of the plug-in bias and Miller-Madow only about 24 %; no estimator is unbiased there, which is why the
# structure tests compare against surrogates rather than against an absolute threshold.
@pytest.mark.parametrize("k", [2, 4, 26])
@pytest.mark.parametrize("n", [200, 2000, 20000])
def test_estimator_bias_on_iid_sequences(k, n):
    rng = rng_for(SEED, "bias", k, n)
    draws = 50
    h0 = {m: np.empty(draws) for m in ("plugin", "grassberger", "miller_madow")}
    drop = {m: np.empty(draws) for m in h0}
    for d in range(draws):
        s = rng.integers(0, k, n)
        for m in h0:
            a = conditional_entropy(s, k, 0, m)
            h0[m][d] = a
            drop[m][d] = a - conditional_entropy(s, k, 1, m)
    for m in ("grassberger", "miller_madow"):
        assert abs(h0[m].mean() - math.log2(k)) < 0.05, (m, h0[m].mean())
    if (k, n) == (26, 200):
        plug = drop["plugin"].mean()
        assert plug > 1.6                                  # the SeaQuest-style plug-in drop (about 1.8 bits)
        assert drop["grassberger"].mean() < 0.5 * plug     # Grassberger removes most of the bias
        assert drop["miller_madow"].mean() < 0.9 * plug    # Miller-Madow removes some of it
    else:
        for m in ("grassberger", "miller_madow"):
            assert abs(drop[m].mean()) < 0.05, (m, drop[m].mean())


def test_known_entropy_rate_of_markov_chain():
    T = np.array([[0.9, 0.1], [0.3, 0.7]])
    pi = np.array([0.75, 0.25])                            # stationary distribution of T
    rate = float(-(pi[:, None] * T * np.log2(T)).sum())    # about 0.572 bits per symbol
    s = markov_chain(T, 20000, rng_for(SEED, "markov-rate"))
    assert conditional_entropy(s, 2, 1) == pytest.approx(rate, abs=0.02)
    curve = entropy_rate_curve(s, 2, 4)
    assert np.all(np.abs(curve[1:] - rate) < 0.03)        # no further drop beyond order 1
    assert curve[0] - curve[1] > 0.2


def test_lagged_mi_corrections():
    rng = rng_for(SEED, "lagged")
    s = rng.integers(0, 8, 2000)
    lags = np.arange(1, 11)
    plug = lagged_mutual_information(s, 8, lags, correction="none")
    mm = lagged_mutual_information(s, 8, lags, correction="miller_madow")
    sh = lagged_mutual_information(s, 8, lags, correction="shuffle", rng=rng, n_shuffle=100)
    assert np.all(plug > 0.01)                             # plug-in bias about 49 / (2 * 2000 ln 2) = 0.018 bits
    assert abs(mm.mean()) < 0.005 and abs(sh.mean()) < 0.005
    with pytest.raises(ValueError):
        lagged_mutual_information(s, 8, lags, correction="shuffle")


def test_compression():
    rng = rng_for(SEED, "compress")
    iid = rng.integers(0, 2, 20000)
    periodic = np.tile(rng.integers(0, 2, 16), 1250)
    assert compression_ratio(iid, 2, "lzma") > 0.95
    assert compression_ratio(periodic, 2, "lzma") < 0.1
    for method in ("zlib", "bz2", "lzma"):
        assert compression_bits(periodic, 2, method) < compression_bits(iid, 2, method)
    text = rng.integers(0, 26, 5000)
    assert compression_ratio(text, 26, "zlib") > 0.8


# ---------------------------------------------------------------------------------------------------------------
# surrogates
# ---------------------------------------------------------------------------------------------------------------

def test_shuffles_keep_what_they_promise():
    rng = rng_for(SEED, "shuffles")
    s = rng.integers(0, 4, 1003)
    assert np.array_equal(np.sort(shuffle(s, rng)), np.sort(s))
    b = block_shuffle(s, 10, rng)
    assert b.size == s.size and np.array_equal(np.sort(b), np.sort(s))
    d = doublet_shuffle(s, 4, rng)
    pairs = lambda x: np.bincount(x[:-1] * 4 + x[1:], minlength=16)
    assert np.array_equal(pairs(d), pairs(s)) and d[0] == s[0] and d[-1] == s[-1]
    assert not np.array_equal(d, s)


def test_permutation_pvalue_formula():
    rng = rng_for(SEED, "pvalue")
    obs, null, p = permutation_pvalue(lambda x: float(x[0]), np.arange(10), shuffle, 99, rng)
    assert obs == 0.0 and null.size == 99
    assert p == pytest.approx((1 + np.count_nonzero(null >= 0)) / 100) and p == 1.0
    _, _, p_low = permutation_pvalue(lambda x: float(x[0]), np.arange(10), shuffle, 99, rng,
                                     larger_is_more_extreme=False)
    # lower tail: only surrogates that start with 0 are as extreme, about 1 in 10
    assert p_low < 0.25


def test_markov_surrogate_keeps_order1_and_destroys_period16():
    T = np.array([[0.8, 0.2], [0.4, 0.6]])
    s = markov_chain(T, 20000, rng_for(SEED, "msur-chain"))
    sur = markov_surrogate(s, 2, 1, rng_for(SEED, "msur-1"))
    assert sur.size == s.size
    assert conditional_entropy(sur, 2, 1) == pytest.approx(conditional_entropy(s, 2, 1), abs=0.03)

    periodic = noisy_word(rng_for(SEED, "msur-word"), repeats=1250, flip=0.02)
    sur16 = markov_surrogate(periodic, 2, 1, rng_for(SEED, "msur-16"))
    mi_orig = lagged_mutual_information(periodic, 2, [16], correction="miller_madow")[0]
    mi_sur = lagged_mutual_information(sur16, 2, [16], correction="miller_madow")[0]
    assert mi_orig > 0.5 and mi_sur < 0.02
    # a Markov surrogate of order 16 or more would keep it; order 3 still does not
    sur3 = markov_surrogate(periodic, 2, 3, rng_for(SEED, "msur-3"))
    assert lagged_mutual_information(sur3, 2, [16], correction="miller_madow")[0] < 0.05


# ---------------------------------------------------------------------------------------------------------------
# n-gram models
# ---------------------------------------------------------------------------------------------------------------

def test_ngram_model_matches_witten_bell_by_hand():
    # train "0 1 0 1": order 0 counts {0: 2, 1: 2}, n = 4, t = 2 -> lam = 2/3; P0(0) = 2/3 * 1/2 + 1/3 * 1/2 = 1/2
    m = NGramModel(1, 2).fit([np.array([0, 1, 0, 1])])
    lp = m.level_log2_probs(np.array([0, 1]))
    assert lp[0, 0] == pytest.approx(-1.0)
    # context "0" seen twice, always followed by 1: n = 2, t = 1, lam = 2/3 -> P1(1 | 0) = 2/3 + 1/3 * 1/2 = 5/6
    assert 2 ** lp[1, 1] == pytest.approx(5 / 6)
    assert m.context_counts(1) == {0: (2, 1), 1: (1, 1)}
    assert m.joint_counts(0) == {0: 2, 1: 2}


def test_held_out_bits_and_gain_on_markov_data():
    T = np.array([[0.9, 0.1], [0.3, 0.7]])
    s = markov_chain(T, 5000, rng_for(SEED, "gain-chain"))
    h1 = held_out_bits_per_symbol(s[:4000], s[4000:], 2, 1)
    h0 = held_out_bits_per_symbol(s[:4000], s[4000:], 2, 0)
    assert h1 < h0 - 0.2
    res = gain_in_bits(s, 2, 1, 100, rng_for(SEED, "gain"))
    assert res["split"] == "blocks" and res["n_units"] == 10
    # truth: H(0.75) - 0.572 = 0.239 bits per symbol
    assert res["gain_bits"] > 0.15 and res["p_value"] < 0.01
    assert res["ci_low"] < res["gain_bits"] < res["ci_high"] and res["ci_low"] > 0.05
    # several sequences are split by sequence
    many = [markov_chain(T, 300, rng_for(SEED, "gain-many", i)) for i in range(12)]
    res_many = gain_in_bits(many, 2, 2, 20, rng_for(SEED, "gain-many"))
    assert res_many["split"] == "sequences" and res_many["p_value"] < 0.01


def test_transition_matrix_pvalues():
    T = np.array([[0.1, 0.9, 0.0], [0.0, 0.1, 0.9], [0.9, 0.0, 0.1]])
    T[0, 2] = T[1, 0] = T[2, 1] = 0.0
    s = markov_chain(T / T.sum(1, keepdims=True), 600, rng_for(SEED, "tm-chain"))
    res = transition_matrix(s, 3, 200, rng_for(SEED, "tm"))
    assert res["counts"].sum() == s.size - 1
    assert np.allclose(res["probabilities"].sum(1), 1.0)
    assert np.all(res["p_adjusted"] >= res["p_values"])
    assert res["p_adjusted"][0, 1] < 0.05 and res["p_adjusted"][0, 2] < 0.05
    iid = rng_for(SEED, "tm-iid").integers(0, 3, 600)
    assert np.all(transition_matrix(iid, 3, 200, rng_for(SEED, "tm2"))["p_adjusted"] > 0.05)
    assert np.allclose(benjamini_hochberg(np.array([0.01, 0.04, 0.03])), [0.03, 0.04, 0.04])


# ---------------------------------------------------------------------------------------------------------------
# false-positive rates
# ---------------------------------------------------------------------------------------------------------------

def _rejection_rate(test, n_sequences: int = 300) -> float:
    rejections = 0
    for i in range(n_sequences):
        rng = rng_for(SEED, "fpr", test.__name__, i)
        if test(rng.integers(0, 4, 500), rng) < 0.05:
            rejections += 1
    return rejections / n_sequences


def _gain_p(seq: np.ndarray, rng: np.random.Generator) -> float:
    return gain_in_bits(seq, 4, 1, 0, rng)["p_value"]


def _grammar_p(seq: np.ndarray, rng: np.random.Generator) -> float:
    return grammar_gain(seq, 4, rng)["p_value"]


def test_gain_in_bits_false_positive_rate():
    rate = _rejection_rate(_gain_p)
    assert 0.02 <= rate <= 0.09, rate


@pytest.mark.slow
def test_grammar_gain_false_positive_rate():
    rate = _rejection_rate(_grammar_p)                    # about 60 s: 300 x 201 Re-Pair runs
    assert 0.02 <= rate <= 0.09, rate


# ---------------------------------------------------------------------------------------------------------------
# power: unary primes and a noisy repeated word
# ---------------------------------------------------------------------------------------------------------------
# Which statistic detects which signal (p-values with 200 surrogates; 1/201 = 0.005 is the smallest possible):
#   unary primes:    grammar_gain (exact Markov-1 null) p = 0.005; the max-lag MI test (lags 1..32, shuffle null)
#                    does NOT detect it (p = 0.22): the zeros are so rare that pairs at short lags carry almost no
#                    information, while the grammar compresses the runs of ones far better than their rearrangements
#                    (the run lengths are all prime, nearly all odd and increase slowly).
#   noisy 16-bit word repeated 250 times, 5 % flips: both detect it, p = 0.005 each (max plug-in MI 0.56 bits).
#   i.i.d. bits: neither (grammar p = 0.24, max-lag MI p = 0.14).
# Measured false-positive rates at p < 0.05 on the 300 i.i.d. sequences: gain_in_bits 6.0 %, grammar_gain 5.0 %.

def test_power_unary_primes():
    primes = unary_primes(600)
    res = grammar_gain(primes, 2, rng_for(SEED, "power", "primes"), n_surrogates=200)
    assert res["p_value"] < 0.01, res["p_value"]


def test_power_noisy_repeated_word():
    bits = noisy_word(rng_for(SEED, "power", "word"))
    res = grammar_gain(bits, 2, rng_for(SEED, "power", "word-grammar"), n_surrogates=200)
    assert res["p_value"] < 0.01
    _, _, p = permutation_pvalue(max_lag_mi, bits, shuffle, 200, rng_for(SEED, "power", "word-mi"))
    assert p < 0.01


def test_power_iid_bits_not_detected():
    bits = rng_for(SEED, "power", "iid").integers(0, 2, 4000)
    res = grammar_gain(bits, 2, rng_for(SEED, "power", "iid-grammar"), n_surrogates=200)
    assert res["p_value"] >= 0.01
    _, _, p = permutation_pvalue(max_lag_mi, bits, shuffle, 200, rng_for(SEED, "power", "iid-mi"))
    assert p >= 0.01


def test_entropy_rate_curve_sees_the_noisy_word_but_not_the_first_order_drop():
    # the SeaQuest entropy-drop control looked one symbol back only; at order 1 the noisy word shows almost no drop,
    # at order 8 a large one
    bits = noisy_word(rng_for(SEED, "curve", "word"))
    curve = entropy_rate_curve(bits, 2, 8)
    assert curve[0] - curve[1] < 0.1
    assert curve[0] - curve[8] > 0.4


# ---------------------------------------------------------------------------------------------------------------
# time-aware lags
# ---------------------------------------------------------------------------------------------------------------

def test_time_binned_mi_finds_the_time_separation():
    # a random binary pattern repeating every P = 20 time units with 5 % noise, observed at an irregular 30 % of
    # integer times; symbols 20 time units apart are dependent, symbols at any other separation below 40 are not
    rng = rng_for(SEED, "time-lag")
    period, n_grid = 20, 20000
    pattern = rng.integers(0, 2, period)
    grid = np.arange(n_grid)
    signal = pattern[grid % period] ^ (rng.random(n_grid) < 0.05)
    keep = rng.random(n_grid) < 0.3
    t, s = grid[keep].astype(float), signal[keep]
    t_jittered = t + rng.uniform(-0.1, 0.1, t.size)       # irregular, non-integer times
    edges = np.arange(1, 31) - 0.5                         # bins centred on separations 1, 2, ..., 29
    res = mutual_information_decay(s, 2, 0, 100, rng_for(SEED, "time-lag-null"), time=t_jittered, lag_edges=edges)
    peak = int(np.nanargmax(res["mi"]))
    assert res["lags"][peak] == pytest.approx(period)
    assert res["p_value"][peak] < 0.01 and res["mi"][peak] > 0.3
    # the index-lag version, on the same data, has no peak at lag 20 and only weak MI anywhere
    index_mi = lagged_mutual_information(s, 2, np.arange(1, 30), correction="shuffle", rng=rng, n_shuffle=100)
    assert int(np.argmax(index_mi)) + 1 != period
    assert index_mi.max() < 0.2 * res["mi"][peak]


# ---------------------------------------------------------------------------------------------------------------
# grammar and segmentation
# ---------------------------------------------------------------------------------------------------------------

def _expand(rules: np.ndarray, sequence: np.ndarray, k: int) -> list[int]:
    memo: dict[int, list[int]] = {}

    def ex(x: int) -> list[int]:
        if x < k:
            return [x]
        if x not in memo:
            a, b = rules[x - k]
            memo[x] = ex(int(a)) + ex(int(b))
        return memo[x]

    return [y for x in sequence.tolist() for y in ex(x)]


def test_repair_round_trip_and_compression():
    rng = rng_for(SEED, "repair")
    for k, n in [(2, 2000), (5, 3000), (40, 1000)]:
        s = rng.integers(0, k, n)
        g = repair(s, k)
        assert _expand(g["rules"], g["sequence"], k) == s.tolist()
        pairs = g["sequence"][:-1] * (k + g["n_rules"]) + g["sequence"][1:]
        assert np.unique(pairs, return_counts=True)[1].max() <= 2   # only overlapping repeats like "a a a" remain
    s = np.array([0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 0])
    g = repair(s, 3)
    assert _expand(g["rules"], g["sequence"], 3) == s.tolist()
    word = np.tile(rng.integers(0, 4, 50), 40)
    gw = repair(word, 4)
    assert gw["description_bits"] < 0.2 * gw["raw_bits"]
    assert kt_code_bits(np.zeros(10)) == 0.0             # one symbol used: nothing to code
    nine_one = np.array([0] * 9 + [1])
    expected = -math.log2(math.gamma(9.5) * math.gamma(1.5) / (math.gamma(0.5) ** 2 * math.gamma(11)))
    assert kt_code_bits(nine_one) == pytest.approx(expected)


def test_tp_segmentation_and_its_null():
    rng = rng_for(SEED, "tp")
    words = [rng.integers(0, 12, 4) for _ in range(6)]
    text = np.concatenate([words[i] for i in rng.integers(0, 6, 300)])
    b = tp_segmentation(text, 12)
    assert b[0] == 0 and b[-1] == text.size and np.all(np.diff(b) > 0)
    res = tp_segmentation_null(text, 12, rng, n=100)
    assert res["p_mean_chunk_frequency"] < 0.05           # recurring "words" in the text
    noise = rng.integers(0, 12, 1200)
    res_noise = tp_segmentation_null(noise, 12, rng, n=100)
    assert res_noise["n_boundaries"] > 100                 # the chunker still cuts noise into pieces ...
    assert res_noise["p_mean_chunk_frequency"] > 0.05      # ... but the null shows the chunks mean nothing


# ---------------------------------------------------------------------------------------------------------------
# ported behaviour and speed
# ---------------------------------------------------------------------------------------------------------------

def test_fit_decay_prefers_exponential_on_markov_chain():
    T = np.array([[0.9, 0.1], [0.3, 0.7]])
    s = markov_chain(T, 20000, rng_for(SEED, "decay-chain"))
    res = mutual_information_decay(s, 2, 20, 200, rng_for(SEED, "decay"))
    assert res["p_value"][0] < 0.01 and res["mi"][0] > 0.1
    fit = fit_decay(res["lags"], res["mi"])
    assert fit["preferred"] == "exponential", fit["aic"]
    assert "weak evidence" in fit["caveat"]
    # a clean power law is recognised as such
    x = np.arange(1, 31, dtype=float)
    assert fit_decay(x, 0.3 * x ** -0.8)["preferred"] in ("power_law", "composite")


@pytest.mark.slow
def test_speed_entropy_rate_curve_million_bits():
    s = rng_for(SEED, "speed").integers(0, 2, 10 ** 6)
    start = time.perf_counter()
    curve = entropy_rate_curve(s, 2, 8)
    assert time.perf_counter() - start < 10.0
    assert np.all(np.abs(curve - 1.0) < 0.01)


@pytest.mark.slow
def test_speed_repair_hundred_thousand_symbols():
    s = rng_for(SEED, "speed-repair").integers(0, 4, 10 ** 5)
    start = time.perf_counter()
    repair(s, 4)
    assert time.perf_counter() - start < 10.0
