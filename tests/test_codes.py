"""Tests for the blind code recogniser (stargazer.decode.gf2 and stargazer.decode.codes).

Test signals are built here from first principles (no stargazer.bench): Hamming(7,4) and (15,11), an even-parity
(8,7) code, a (3,1) repetition code, the PRBS-style LFSRs x^15 + x^14 + 1 and x^23 + x^18 + 1, a framer with the
CCSDS sync word 0x1ACFFC1D, and a binary symmetric channel.
"""
from __future__ import annotations

import itertools
import time

import numpy as np
import pytest

from stargazer.core.rng import rng_for
from stargazer.core.types import BitStream
from stargazer.decode import gf2
from stargazer.decode.codes import block_code_scan, describe, frame_scan, lfsr_test, sync_word

SEED = 20261122
SYNC = np.array([int(c) for c in f"{0x1ACFFC1D:032b}"], dtype=np.uint8)


# ------------------------------------------------------------------------------------------------ signal helpers

def _systematic(parity: np.ndarray) -> np.ndarray:
    k = parity.shape[0]
    return np.concatenate([np.eye(k, dtype=np.uint8), parity.astype(np.uint8)], axis=1)


G_HAMMING74 = _systematic(np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]]))
G_HAMMING1511 = _systematic(np.array([c for c in itertools.product([0, 1], repeat=4) if sum(c) >= 2]))
G_PARITY87 = _systematic(np.ones((7, 1), dtype=np.uint8))
G_REP31 = np.ones((1, 3), dtype=np.uint8)


def encode(generator: np.ndarray, n_words: int, rng: np.random.Generator) -> np.ndarray:
    """Random messages through a linear code, concatenated."""
    msgs = rng.integers(0, 2, (n_words, generator.shape[0]), dtype=np.uint8)
    return ((msgs.astype(int) @ generator) % 2).astype(np.uint8).ravel()


def lfsr(taps: tuple[int, int], n: int, rng: np.random.Generator) -> np.ndarray:
    """Fibonacci LFSR s[i] = s[i-a] xor s[i-b] (tap notation: x^b + x^a + 1 with b > a), random non-zero seed."""
    a, b = taps
    s = np.zeros(n, dtype=np.uint8)
    seed = rng.integers(0, 2, b, dtype=np.uint8)
    seed[0] = 1
    s[:b] = seed
    for i in range(b, n):
        s[i] = s[i - a] ^ s[i - b]
    return s


def lfsr_poly(taps: tuple[int, int]) -> np.ndarray:
    """Connection polynomial (1 + x^a + x^b) of the generator above, as Berlekamp-Massey returns it."""
    c = np.zeros(taps[1] + 1, dtype=np.uint8)
    c[[0, taps[0], taps[1]]] = 1
    return c


def framed(n_frames: int, rng: np.random.Generator, payload: int = 480) -> np.ndarray:
    """n_frames frames of SYNC followed by `payload` random bits (512-bit frames)."""
    frames = [np.concatenate([SYNC, rng.integers(0, 2, payload, dtype=np.uint8)]) for _ in range(n_frames)]
    return np.concatenate(frames)


def channel(bits: np.ndarray, ber: float, rng: np.random.Generator) -> np.ndarray:
    """Binary symmetric channel."""
    return (bits ^ (rng.random(bits.size) < ber)).astype(np.uint8)


def random_bits(n: int, rng: np.random.Generator) -> np.ndarray:
    return rng.integers(0, 2, n, dtype=np.uint8)


# -------------------------------------------------------------------------------------------------------- gf2

def _slow_rank(m: np.ndarray) -> int:
    """Textbook Gauss-Jordan elimination on a uint8 matrix, one row at a time."""
    a = m.copy() % 2
    rows, cols = a.shape
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if a[i, c]), None)
        if piv is None:
            continue
        a[[r, piv]] = a[[piv, r]]
        for i in range(rows):
            if i != r and a[i, c]:
                a[i] ^= a[r]
        r += 1
        if r == rows:
            break
    return r


def test_gf2_rank_null_space_against_reference() -> None:
    g = rng_for(SEED, "gf2", "matrices")
    for trial in range(200):
        rows = int(g.integers(1, 90))
        cols = int(g.integers(1, 140))
        density = float(g.choice([0.1, 0.5, 0.9]))
        m = (g.random((rows, cols)) < density).astype(np.uint8)
        if trial % 4 == 0 and rows > 2:                 # force dependencies
            m[-1] = m[0] ^ m[1]
        ref_rank = _slow_rank(m)
        assert gf2.rank(m) == ref_rank
        h = gf2.null_space(m)
        assert h.shape == (cols - ref_rank, cols)
        assert not ((m.astype(int) @ h.T.astype(int)) % 2).any()
        assert gf2.rank(h) == h.shape[0] if h.size else True
        reduced, pivots = gf2.row_reduce(m)
        assert len(pivots) == ref_rank
        assert np.array_equal(gf2.unpack_rows(gf2.pack_rows(m), cols), m)
        if cols <= 10:                                    # brute-force null-space size
            vecs = np.array(list(itertools.product([0, 1], repeat=cols)), dtype=int)
            assert int((((m.astype(int) @ vecs.T) % 2).sum(axis=0) == 0).sum()) == 2 ** (cols - ref_rank)


def test_batched_null_space_matches_single() -> None:
    g = rng_for(SEED, "gf2", "batched")
    for n in (3, 7, 20, 64):
        mats = g.integers(0, 2, (30, n - 1, n), dtype=np.uint8)
        packed = np.stack([gf2.pack_rows(mm)[:, 0] for mm in mats])
        basis, valid = gf2.batched_null_space(packed, n)
        for b in range(30):
            ns = gf2.null_space(mats[b])
            assert valid[b].sum() == ns.shape[0]
            got = gf2.unpack_rows(basis[b][valid[b]][:, None], n)
            assert not ((mats[b].astype(int) @ got.T.astype(int)) % 2).any()
            assert gf2.rank(got) == ns.shape[0] if ns.size else True


def test_berlekamp_massey() -> None:
    g = rng_for(SEED, "gf2", "bm")
    s = lfsr((14, 15), 200, g)
    big_l, poly = gf2.berlekamp_massey(s)
    assert big_l == 15
    assert np.array_equal(poly, lfsr_poly((14, 15)))
    # The target "L within 2 of n/2 for all 50 streams" is not a property random bits have: for n = 2000 exactly,
    # P(L <= 997) = (2 * 4^997 + 1) / (3 * 2^2000) = 0.0104 and P(L >= 1003) = sum_{d>=3} 4^-d = 0.0208, so a stream
    # falls 3 or more away with probability 0.031 and all 50 stay within 2 only 21% of the time. Measured with this
    # seed: deviations -3 x1, -2 x1, -1 x4, 0 x22, +1 x17, +2 x4, +3 x1 (48 of 50 within 2). We therefore assert the
    # distributional facts: at least 45 within 2 (expected 48.4), none beyond 6 (P = 1.2e-4 per stream), and a mean
    # deviation near the theoretical +2/9.
    deviations = []
    for _ in range(50):
        big_l, poly = gf2.berlekamp_massey(random_bits(2000, g))
        assert poly.size == big_l + 1 and poly[0] == 1
        deviations.append(big_l - 1000)
    dev = np.array(deviations)
    assert np.sum(np.abs(dev) <= 2) >= 45
    assert np.all(np.abs(dev) <= 6)
    assert abs(dev.mean() - 2 / 9) < 0.5
    prof = gf2.linear_complexity_profile(s)
    assert prof[-1] == 15 and prof.size == 200 and np.all(np.diff(prof) >= 0)
    assert np.array_equal(prof, [gf2.berlekamp_massey(s[: i + 1])[0] for i in range(200)])
    t0 = time.perf_counter()
    big_l, _ = gf2.berlekamp_massey(random_bits(100_000, g))
    assert time.perf_counter() - t0 < 5.0
    assert abs(big_l - 50_000) <= 6


# -------------------------------------------------------------------------------------------- frame and sync

CUT = 20                       # the stream starts 20 bits into a frame, so the sync word wraps the frame boundary
PHASE = 512 - CUT


def _framed_stream(ber: float, tag: str) -> np.ndarray:
    g = rng_for(SEED, "framed", tag)
    return channel(framed(40, g), ber, g)[CUT:]


@pytest.mark.parametrize("ber", [0.0, 0.01, 0.05])
def test_frame_scan_finds_512(ber: float) -> None:
    res = frame_scan(BitStream(_framed_stream(ber, f"frame{ber}")), rng=rng_for(SEED, "frame_scan", ber))
    assert res.frame_len == 512
    assert res.p_value < 0.01 and res.p_bonferroni < 1e-20
    assert np.all(res.harmonic_z > 3)


def test_frame_scan_respects_known_mask() -> None:
    bits = _framed_stream(0.01, "erased")
    g = rng_for(SEED, "erasures")
    known = g.random(bits.size) > 0.2
    garbage = np.where(known, bits, g.integers(0, 2, bits.size)).astype(np.uint8)
    res = frame_scan(BitStream(garbage, known), n_surrogates=0)
    assert res.frame_len == 512 and res.p_value < 1e-10
    assert res.n_pairs[0] < 0.7 * (bits.size - res.lags[0])     # erased pairs were not counted


def test_frame_scan_false_alarms_bonferroni() -> None:
    """100 random streams, analytic (Bonferroni) correction; the surrogate version is the slow test below."""
    g = rng_for(SEED, "frame_noise")
    p = np.array([frame_scan(BitStream(random_bits(20_000, g)), n_surrogates=0).p_value for _ in range(100)])
    assert np.sum(p > 0.01) >= 98


@pytest.mark.slow
def test_frame_scan_false_alarms_surrogates() -> None:
    g = rng_for(SEED, "frame_noise")
    p = np.array([frame_scan(BitStream(random_bits(20_000, g)), rng=rng_for(SEED, "frame_sur", i)).p_value
                  for i in range(100)])
    assert np.sum(p > 0.01) >= 98


@pytest.mark.parametrize("ber", [0.0, 0.01])
def test_sync_word_exact(ber: float) -> None:
    res = sync_word(BitStream(_framed_stream(ber, f"sync{ber}")), 512)
    assert res.phase == PHASE
    assert res.run_span == 32 and res.run_length == 32
    assert np.array_equal(res.pattern, SYNC)
    assert res.pattern_hex == "0x1ACFFC1D"
    assert res.p_value < 1e-20


def test_sync_word_five_percent() -> None:
    res = sync_word(BitStream(_framed_stream(0.05, "sync5")), 512)
    true_res = (PHASE + np.arange(32)) % 512
    got = dict(zip(res.pattern_positions.tolist(), res.pattern.tolist()))
    constant = set(res.positions.tolist())
    correct = sum(1 for r, b in zip(true_res.tolist(), SYNC.tolist()) if r in constant and got.get(r) == b)
    assert correct >= 28
    assert res.p_value < 1e-10


def test_sync_word_on_noise() -> None:
    res = sync_word(BitStream(random_bits(20_480, rng_for(SEED, "sync_noise"))), 512)
    assert res.p_value > 0.01


# ------------------------------------------------------------------------------------------------ block codes

def _coded(generator: np.ndarray, n_words: int, ber: float, tag: str, stray: int = 0) -> np.ndarray:
    g = rng_for(SEED, "coded", tag, ber)
    return np.concatenate([random_bits(stray, g), channel(encode(generator, n_words, g), ber, g)])


def _assert_top(res: list, n: int, k: int, offset: int, generator: np.ndarray) -> None:
    assert res, "nothing detected"
    top = res[0]
    assert (top.n, top.k, top.offset) == (n, k, offset)
    assert top.p_value < 1e-10
    # Every reported check is orthogonal to the code: G H^T = 0.
    assert not ((generator.astype(int) @ top.H.T.astype(int)) % 2).any()


@pytest.mark.parametrize("ber", [0.0, 0.001, 0.01])
def test_block_code_hamming74(ber: float) -> None:
    x = _coded(G_HAMMING74, 2000, ber, "h74", stray=3)
    res = block_code_scan(BitStream(x), rng=rng_for(SEED, "scan", "h74", ber))
    _assert_top(res, 7, 4, 3, G_HAMMING74)
    assert abs(res[0].ber_estimate - ber) < 0.004
    assert (14, 3) in res[0].collapsed                          # the multiple of 7 at the same alignment


@pytest.mark.parametrize("ber", [0.0, 0.001])
def test_block_code_hamming1511(ber: float) -> None:
    x = _coded(G_HAMMING1511, 934, ber, "h1511")
    res = block_code_scan(BitStream(x), rng=rng_for(SEED, "scan", "h1511", ber))
    _assert_top(res, 15, 11, 0, G_HAMMING1511)


def test_block_code_parity_and_repetition() -> None:
    res = block_code_scan(BitStream(_coded(G_PARITY87, 1750, 0.0, "par")), rng=rng_for(SEED, "scan", "par"))
    _assert_top(res, 8, 7, 0, G_PARITY87)
    assert res[0].weights.tolist() == [8]
    res = block_code_scan(BitStream(_coded(G_REP31, 4667, 0.0, "rep")), rng=rng_for(SEED, "scan", "rep"))
    _assert_top(res, 3, 1, 0, G_REP31)
    assert res[0].exact_rank == 1


def _noise_scan_pvalues(count: int) -> np.ndarray:
    g = rng_for(SEED, "block_noise")
    out = []
    for i in range(count):
        res = block_code_scan(BitStream(random_bits(14_000, g)), rng=rng_for(SEED, "block_noise_scan", i))
        out.append(min((r.p_value for r in res), default=1.0))
    return np.array(out)


def test_block_code_false_alarms_quick() -> None:
    """20 random streams in the default run (about 1 s each); the 100-stream version is marked slow."""
    assert np.sum(_noise_scan_pvalues(20) >= 0.01) == 20


@pytest.mark.slow
def test_block_code_false_alarms() -> None:
    assert np.sum(_noise_scan_pvalues(100) >= 0.01) >= 98


# ------------------------------------------------------------------------------------------------------ LFSR

@pytest.mark.parametrize("taps", [(14, 15), (18, 23)])
def test_lfsr_clean(taps: tuple[int, int]) -> None:
    s = lfsr(taps, 5000, rng_for(SEED, "lfsr", taps))
    res = lfsr_test(BitStream(s), rng=rng_for(SEED, "lfsr_test", taps))
    assert res.L == taps[1]
    assert np.array_equal(res.polynomial, lfsr_poly(taps))
    assert res.regen_mismatches == 0
    assert res.p_global < 1e-100 and res.p_value < 0.01


@pytest.mark.parametrize("taps", [(14, 15), (18, 23)])
def test_lfsr_noisy_windows(taps: tuple[int, int]) -> None:
    g = rng_for(SEED, "lfsr_noisy", taps)
    clean = lfsr(taps, 5000, g)
    s = channel(clean, 0.002, g)
    res = lfsr_test(BitStream(s), rng=rng_for(SEED, "lfsr_noisy_test", taps))
    assert res.p_global > 0.01                                 # a few flipped bits ruin the global estimate
    assert res.window_L_min == taps[1]
    assert np.array_equal(res.window_polynomial, lfsr_poly(taps))
    assert res.window_agree >= 5
    assert res.p_window < 0.01 and res.p_value < 0.01
    assert res.window_regen_mismatches == int(np.sum(s != clean))   # regeneration recovers the clean stream


def test_lfsr_false_alarms() -> None:
    g = rng_for(SEED, "lfsr_noise")
    p = np.array([lfsr_test(BitStream(random_bits(5000, g)), rng=rng_for(SEED, "lfsr_sur", i)).p_value
                  for i in range(100)])
    assert np.sum(p > 0.01) >= 98


# ------------------------------------------------------------------------------- the hard case: LFSR xor code

def test_lfsr_xor_hamming() -> None:
    """The 15-bit LFSR added (xor) to random Hamming(7,4) codewords, 14,000 bits.

    What each recogniser reports, and why (the same on every seed tried):
    * frame_scan: nothing (p ~ 0.5). There is no repeated pattern, and an m-sequence has a flat autocorrelation.
    * lfsr_test: nothing. L = n/2 = 7000 globally and about 125 in every 256-bit window, because the codeword bits
      are random; no two windows share a polynomial.
    * block_code_scan: finds linear redundancy, but not at n = 7. Inside any window the stream is a vector of the
      15-dimensional space of LFSR windows plus whole Hamming codewords, so aligned blocks of n = 7m bits lie in a
      space of dimension at most 15 + 4m, which is below n once m >= 6; for this trinomial the LFSR windows meet
      the Hamming check space earlier, and the exact ranks are 20 at n = 21, 26 at 28, 32 at 35 and 43 at 49. The
      scan reports these as error-free codes, at offsets that are multiples of 7 (and slightly weaker ones at the
      other offsets of n = 49, 35, 28, because the LFSR part does not care about alignment). At n = 7 itself the
      LFSR randomises every syndrome and nothing is found. n = 42, 56, 63 are folded in as multiples.
    * Not done by the recognisers, but it would peel both layers: the three Hamming syndrome sequences of the
      aligned 7-bit blocks are themselves LFSR sequences (Berlekamp-Massey gives linear complexity 15 for each).
    """
    g = rng_for(SEED, "xor")
    x = lfsr((14, 15), 14_000, g) ^ encode(G_HAMMING74, 2000, g)
    stream = BitStream(x)
    fr = frame_scan(stream, rng=rng_for(SEED, "xor_frame"))
    assert fr.p_value > 0.01
    lf = lfsr_test(stream, rng=rng_for(SEED, "xor_lfsr"))
    assert lf.p_value > 0.01 and lf.L >= 6990 and lf.window_agree <= 1
    res = block_code_scan(stream, rng=rng_for(SEED, "xor_scan"))
    assert res and res[0].p_value < 1e-10
    assert all(r.n % 7 == 0 and r.n >= 21 for r in res if r.offset % 7 == 0)
    assert not any(r.n < 21 for r in res)
    exact = {r.n: r.exact_rank for r in res if r.offset == 0}
    assert exact.get(21) == 20 and exact.get(28) == 26 and exact.get(35) == 32 and exact.get(49) == 43
    blocks = x.reshape(-1, 7)
    syndromes = (blocks.astype(int) @ gf2.null_space(G_HAMMING74).T.astype(int)) % 2
    assert [gf2.berlekamp_massey(syndromes[:, j])[0] for j in range(3)] == [15, 15, 15]


# ------------------------------------------------------------------------------------------------- describe

def test_describe_is_plain_and_finds_the_code() -> None:
    import json

    x = _coded(G_HAMMING74, 2000, 0.001, "describe", stray=3)
    out = describe(BitStream(x), rng=rng_for(SEED, "describe"))
    json.dumps(out)
    assert out["block_codes"][0]["n"] == 7 and out["block_codes"][0]["k"] == 4
    assert out["significant"] == ["block_code"]
