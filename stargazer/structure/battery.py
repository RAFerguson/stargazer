"""One question, several statistics, one probability: does this sequence have structure its null lacks?

Each statistic looks for a different kind of order:

* `entropy_drop`: how much knowing the last few symbols helps predict the next (short-range rules);
* `lzma_saving`, `bz2_saving`: what general-purpose compressors save (repeats and long matches);
* `lag_peak`: the strongest dependence between symbols any fixed distance apart (frames, rows, periods);
* `grammar_saving`: what a hierarchical grammar saves (nested, reused units).

Running five tests and reporting the best would inflate significance fivefold. Instead every statistic is
computed on the sequence and on each surrogate, each value is turned into its rank among all of them, and the
sequence's best rank is compared with every surrogate's best rank. The single probability that results is
exact whatever the dependence between the statistics.

The null matters more than the statistics. Against shuffled symbols almost every natural process shows
"structure", because natural processes have memory. The default null here is therefore a Markov chain fitted
to the sequence itself, with its order chosen by the Bayesian information criterion up to a small maximum:
a short-memory random process with exactly the sequence's own short-range habits. What this battery detects
is structure that such a process does not produce: framing, long-range order, hierarchy. It is not a test of
origin; a clock or a crystal passes it.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np

from . import estimators, grammar, surrogates


def markov_order_bic(seq: np.ndarray, k: int, max_order: int = 3) -> int:
    """The Markov order (0..max_order) that minimises the Bayesian information criterion."""
    seq = np.asarray(seq)
    n = seq.size
    limit = estimators.max_trustworthy_order(n, k, max_order)
    curve = estimators.entropy_rate_curve(seq, k, limit, method="plugin")
    best, best_bic = 0, np.inf
    for order in range(limit + 1):
        bic = 2.0 * n * curve[order] * np.log(2.0) + (k ** order) * (k - 1) * np.log(n)
        if bic < best_bic:
            best, best_bic = order, bic
    return best


def lag_peak(seq: np.ndarray, k: int, max_lag: int | None = None) -> tuple[float, int]:
    """The strongest dependence at any single lag, and the lag.

    For two symbols: the largest absolute autocorrelation of the +-1 sequence over lags 1..max_lag, in units
    of its standard error (computed for all lags at once with a Fourier transform). For larger alphabets: the
    largest plug-in mutual information between symbols one lag apart, in bits, over a smaller range of lags.
    """
    seq = np.asarray(seq)
    n = seq.size
    if k == 2:
        max_lag = min(max_lag or 4096, n // 4)
        x = 2.0 * seq - 1.0
        x = x - x.mean()
        size = int(2 ** np.ceil(np.log2(2 * n)))
        spectrum = np.fft.rfft(x, size)
        acf = np.fft.irfft(spectrum * np.conj(spectrum), size)[1:max_lag + 1]
        lags = np.arange(1, max_lag + 1)
        z = np.abs(acf) / (np.var(x) * np.sqrt(n - lags) + 1e-300)
        best = int(np.argmax(z))
        return float(z[best]), int(lags[best])
    max_lag = min(max_lag or 64, n // 8)
    lags = np.arange(1, max_lag + 1)
    pairings, _ = estimators.lag_pair_index(n, lags)
    mi, _, _ = estimators._mi_curve(seq, k, pairings)
    best = int(np.argmax(mi))
    return float(mi[best]), int(lags[best])


def statistics(seq: np.ndarray, k: int, *, include_grammar: bool = True) -> dict[str, float]:
    """Every statistic of the battery for one sequence. Larger always means more structure."""
    seq = np.asarray(seq)
    n = seq.size
    order = estimators.max_trustworthy_order(n, k, 8 if k == 2 else 3)
    curve = estimators.entropy_rate_curve(seq, k, order, method="grassberger")
    raw_bits = n * np.log2(k)
    out = {"entropy_drop": float(curve[0] - np.nanmin(curve)),
           "lzma_saving": float(raw_bits - estimators.compression_bits(seq, k, "lzma")),
           "bz2_saving": float(raw_bits - estimators.compression_bits(seq, k, "bz2")),
           "lag_peak": lag_peak(seq, k)[0]}
    if include_grammar:
        out["grammar_saving"] = float(grammar._grammar_saving(seq, k))
    return out


@dataclass
class BatteryResult:
    n: int
    k: int
    null: str                         # e.g. "markov order 1" or "shuffle"
    null_order: int
    n_surrogates: int
    p_value: float                    # the omnibus probability; at least 1 / (n_surrogates + 1)
    detected: bool
    alpha: float
    statistics: dict[str, dict] = field(default_factory=dict)   # name -> observed, null mean, null sd, p
    strongest: str = ""
    lag_of_peak: int = 0

    def to_dict(self) -> dict:
        return asdict(self)


def battery(seq: np.ndarray, k: int, rng: np.random.Generator, *, null: str = "markov", max_null_order: int = 3,
            n_surrogates: int = 199, alpha: float = 0.01, include_grammar: bool | None = None) -> BatteryResult:
    """Test `seq` (integers in 0..k-1) for structure beyond its null.

    null: "markov" (order chosen by `markov_order_bic` up to `max_null_order`) or "shuffle" (order 0 with the
    symbol counts kept exactly). The grammar statistic is skipped by default for sequences longer than 60,000
    symbols, where it dominates the cost.
    """
    seq = np.asarray(seq)
    if include_grammar is None:
        include_grammar = seq.size <= 60_000
    if null == "shuffle":
        order, name = 0, "shuffle"

        def draw() -> np.ndarray:
            return surrogates.shuffle(seq, rng)
    elif null == "markov":
        order = markov_order_bic(seq, k, max_null_order)
        name = f"markov order {order}"

        def draw() -> np.ndarray:
            return surrogates.markov_surrogate(seq, k, order, rng)
    else:
        raise ValueError('null must be "markov" or "shuffle"')
    observed = statistics(seq, k, include_grammar=include_grammar)
    names = list(observed)
    table = np.empty((n_surrogates + 1, len(names)))
    table[0] = [observed[name_] for name_ in names]
    for i in range(1, n_surrogates + 1):
        values = statistics(draw(), k, include_grammar=include_grammar)
        table[i] = [values[name_] for name_ in names]
    table = np.nan_to_num(table, nan=-np.inf)
    # Rank-based probability of every entry within its column, then the best per row, then the rank of row 0.
    upper = (table[None, :, :] >= table[:, None, :]).sum(axis=1) / (n_surrogates + 1.0)
    best = upper.min(axis=1)
    p_value = float(np.sum(best <= best[0]) / (n_surrogates + 1.0))
    detail = {}
    for j, name_ in enumerate(names):
        null_values = table[1:, j]
        finite = null_values[np.isfinite(null_values)]
        detail[name_] = {"observed": float(table[0, j]),
                         "null_mean": float(finite.mean()) if finite.size else float("nan"),
                         "null_sd": float(finite.std()) if finite.size else float("nan"),
                         "p": float(upper[0, j])}
    strongest = names[int(np.argmin(upper[0]))]
    return BatteryResult(int(seq.size), int(k), name, int(order), int(n_surrogates), p_value, p_value <= alpha,
                         alpha, detail, strongest, lag_peak(seq, k)[1])
