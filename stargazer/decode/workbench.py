"""The decoding workbench: what, if anything, can be said about the content of an unknown signal.

Give it a bit stream, a symbol sequence, a light curve or a list of event times. It climbs the content
ladder (`core.tiers.ContentTier`) one rung at a time and stops where the evidence stops:

1. **Symbols.** Does the signal sit on a few discrete states (levels of a light curve, slots of a clock)?
2. **Structure.** Does the sequence have order that a short-memory random process fitted to it does not?
3. **Layout.** Is there a geometry or a code: a raster width, a frame with a sync word, parity checks, a
   shift register?
4. **Anchor.** Does any of it match knowledge a sender and receiver could share without meeting: counting,
   primes, arithmetic, mathematical constants, or (for human signals) a known character code?
5. **Reading.** Is there a specific interpretation that predicts held-out parts of the signal?

Rules the workbench keeps:

* Every test's probability is corrected for the search inside the test, and the set of tests is fixed in
  advance, so the number of things tried is known and reported (`Reading.tried`).
* Layout claims must beat two nulls: independent bits (analytically) and a short-memory chain fitted to the
  stream (by simulation). Natural processes have memory; memory alone must not look like engineering.
* If a complete generative explanation is found (a shift register that reproduces every bit), nothing above
  layout is claimed: every other pattern in the stream is a consequence of that one.
* A strictly repeating sequence is reduced to one cycle. A clock carries no information beyond its period.
* The output says what was not found and what that does and does not mean. Compressed or encrypted content
  is indistinguishable from noise by design, so "nothing found" is never "nothing there".

The workbench decides how much structure has been recovered. It does not decide where a signal came from.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np

from ..core import BitStream, ContentTier, EventList, Series, SymbolSequence, content_phrase, rng_for
from ..structure import surrogates
from ..structure.battery import BatteryResult, battery, markov_order_bic
from . import anchors, codes, layout, symbolise


@dataclass(frozen=True)
class Settings:
    """Thresholds, fixed before looking at any signal."""

    alpha_structure: float = 0.01      # the structure battery (one omnibus probability)
    alpha_layout: float = 1e-4         # each layout test, after its own internal correction
    alpha_anchor: float = 1e-3         # anchors, after correction for the whole anchor library
    fold_threshold_bits: float = 20.0  # code-length gain a raster must reach
    markov_margin: float = 6.0         # a layout statistic must exceed short-memory surrogates by this many sd
    n_markov: int = 20
    n_surrogates: int = 199            # for the structure battery
    block_n_max: int = 32              # largest block length searched for parity checks
    min_bits: int = 128
    min_reading_accuracy: float = 0.9  # held-out accuracy a reading needs
    min_text_letters: float = 0.75     # share of letters and spaces a text reading needs


@dataclass
class Evidence:
    rung: int
    test: str
    p_value: float
    sentence: str
    detail: dict = field(default_factory=dict)


@dataclass
class Reading:
    name: str
    kind: str
    n: int
    tier: int = 0
    tier_phrase: str = ""
    evidence: list[Evidence] = field(default_factory=list)
    tried: list[dict] = field(default_factory=list)      # every test run, found or not
    layout: dict = field(default_factory=dict)
    picture: str = ""
    code_lengths: dict = field(default_factory=dict)     # bits needed to write the stream down under each model
    statement: str = ""
    limits: str = ""
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    def add(self, rung: int, test: str, p_value: float, sentence: str, **detail: Any) -> None:
        self.evidence.append(Evidence(int(rung), test, float(p_value), sentence, detail))

    def tested(self, test: str, p_value: float, found: bool, note: str = "") -> None:
        self.tried.append({"test": test, "p_value": float(p_value), "found": bool(found), "note": note})


# ------------------------------------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------------------------------------

def _markov_margin(statistic, bits: np.ndarray, order: int, rng: np.random.Generator, n: int) -> float:
    """How many standard deviations the stream's statistic lies above the same statistic on short-memory
    surrogates of the stream."""
    observed = float(statistic(bits))
    null = np.array([statistic(surrogates.markov_surrogate(bits, 2, order, rng)) for _ in range(n)])
    return float((observed - null.mean()) / max(null.std(), 1e-9))


def _p(value: Any) -> float:
    """A probability as a float, whether it arrived as a number or as a Monte Carlo result."""
    return float(getattr(value, "p", value))


def _lag_agreement(bits: np.ndarray, lag: int) -> float:
    return float(np.mean(bits[:-lag] == bits[lag:])) if 0 < lag < bits.size else 0.0


def _check_satisfaction(bits: np.ndarray, n: int, offset: int, check: np.ndarray) -> float:
    blocks = bits[offset:offset + ((bits.size - offset) // n) * n].reshape(-1, n)
    return float(np.mean((blocks[:, check.astype(bool)].sum(axis=1) % 2) == 0)) if blocks.size else 0.0


def _text_like(preview: str, space_most_common: bool, settings: Settings) -> bool:
    """A text reading must look like words: mostly letters and spaces, with the space the commonest character.
    Printable characters alone are not enough; several codes leave printable gibberish."""
    if not preview:
        return False
    letters = sum(ch.isalpha() or ch == " " for ch in preview) / len(preview)
    return bool(space_most_common and letters >= settings.min_text_letters)


def _periodic_unit(symbols: np.ndarray, k: int) -> tuple[int, float, np.ndarray] | None:
    """If the sequence repeats, the period, the share of symbols the repeating unit explains, and the unit
    (majority symbol at each phase). None if no period explains at least 98% of the sequence."""
    n = symbols.size
    if n < 16:
        return None
    best = None
    max_lag = n // 3
    agreement = np.array([np.mean(symbols[:-lag] == symbols[lag:]) for lag in range(1, max_lag + 1)])
    for lag in np.argsort(-agreement)[:5] + 1:
        phase = np.arange(n) % lag
        counts = np.zeros((lag, k))
        np.add.at(counts, (phase, symbols), 1)
        unit = counts.argmax(axis=1)
        explained = float(np.mean(unit[phase] == symbols))
        if explained >= 0.98 and (best is None or lag < best[0]):
            best = (int(lag), explained, unit)
    return best


# ------------------------------------------------------------------------------------------------------------
# Bit streams
# ------------------------------------------------------------------------------------------------------------

def _analyse_bits(stream: BitStream, reading: Reading, rng_seed: tuple, settings: Settings) -> None:
    bits = np.asarray(stream.bits, dtype=np.uint8)
    n = bits.size
    if n < settings.min_bits:
        reading.notes.append(f"Only {n} bits: too few for any of the tests.")
        return
    order = markov_order_bic(bits, 2, 3)

    # Rung 2: structure beyond a short-memory process.
    result: BatteryResult = battery(bits, 2, rng_for(*rng_seed, "battery"), n_surrogates=settings.n_surrogates,
                                    alpha=settings.alpha_structure)
    reading.tested("structure battery", result.p_value, result.detected, f"null: {result.null}")
    if result.detected:
        reading.add(ContentTier.STRUCTURE, "structure battery", result.p_value,
                    f"The stream has order that a {result.null} chain fitted to it does not produce (probability "
                    f"{result.p_value:.3f} over {result.n_surrogates} simulated chains; strongest sign: "
                    f"{result.strongest.replace('_', ' ')}).", statistics=result.statistics)

    # Rung 3a: a shift register. If it reproduces the whole stream, that is the whole story.
    lfsr = codes.lfsr_test(stream, rng=rng_for(*rng_seed, "lfsr"))
    exact = bool(lfsr.L >= 3 and lfsr.p_global <= settings.alpha_layout and lfsr.regen_mismatches == 0)
    reading.tested("shift register", lfsr.p_global, exact, f"linear complexity {lfsr.L} of {lfsr.n} bits")
    reading.code_lengths["raw"] = float(n)
    if exact:
        reading.add(ContentTier.LAYOUT, "shift register", lfsr.p_global,
                    f"Every one of the {n:,} bits follows from the first {lfsr.L} by one fixed rule (a linear "
                    f"feedback shift register of length {lfsr.L}). Random bits would need a register about "
                    f"{n // 2:,} long.", polynomial=codes._poly_text(lfsr.polynomial), length=int(lfsr.L))
        reading.layout = {"kind": "shift register", "length": int(lfsr.L),
                          "polynomial": codes._poly_text(lfsr.polynomial)}
        reading.code_lengths["shift register"] = float(2 * lfsr.L + 2 * np.log2(n))
        reading.notes.append("A shift register that reproduces the whole stream explains every other pattern in "
                             "it, so no further tests were made. Such sequences are used as test patterns, "
                             "spreading codes and scramblers; this one carries no content beyond its rule.")
        return

    # Rung 3b: a raster.
    fold = layout.fold_scan(stream, threshold_bits=settings.fold_threshold_bits)
    reading.tested("raster width", fold.p_bound, fold.detected, f"best width {fold.width}, gain {fold.gain_bits:.0f} bits")
    reading.code_lengths["short-memory model"] = float(fold.baseline_bits)
    raster = False
    if fold.detected:
        raster = True
        rows = (n - fold.header_bits) // fold.width
        sentence = (f"Written as rows of {fold.width} bits the stream takes {fold.gain_bits:,.0f} fewer bits to "
                    f"describe than as a line with short memory, after paying for the search over widths: "
                    f"neighbouring rows resemble each other. That is a raster {fold.width} wide and {rows} deep")
        sentence += f", after {fold.header_bits} leading bits that do not belong to it." if fold.header_bits else "."
        if fold.factor_note:
            sentence += f" ({fold.factor_note}.)"
        reading.add(ContentTier.LAYOUT, "raster width", fold.p_bound, sentence, width=int(fold.width),
                    gain_bits=float(fold.gain_bits), header_bits=int(fold.header_bits))
        reading.layout = {"kind": "raster", "width": int(fold.width), "rows": int(rows),
                          "header_bits": int(fold.header_bits)}
        reading.code_lengths[f"raster of width {fold.width}"] = float(fold.folded_bits)
        volume = layout.volume_scan(stream, fold.width, threshold_bits=settings.fold_threshold_bits)
        reading.tested("volume height", volume.p_bound, volume.detected, f"best height {volume.height}")
        if volume.detected:
            reading.add(ContentTier.LAYOUT, "volume height", volume.p_bound,
                        f"The rows in turn group into planes of {volume.height} rows: every plane resembles the "
                        f"one before it, saving a further {volume.gain_bits:,.0f} bits. That is a block "
                        f"{fold.width} x {volume.height} x {volume.depth_layers}.", height=int(volume.height),
                        layers=int(volume.depth_layers))
            reading.layout.update({"kind": "volume", "height": int(volume.height),
                                   "layers": int(volume.depth_layers)})
        image = layout.render(stream, fold.width, offset=fold.header_bits)
        if image.shape[0] <= 200 and fold.width <= 160:
            reading.picture = layout.to_text(image)

    # Rung 3c: a frame, by agreement of the stream with itself one frame later.
    frame = codes.frame_scan(stream, rng=rng_for(*rng_seed, "frame"), n_surrogates=0)
    frame_found = False
    if frame.p_bonferroni <= settings.alpha_layout:
        margin = _markov_margin(lambda b: _lag_agreement(b, frame.frame_len), bits, order,
                                rng_for(*rng_seed, "frame-markov"), settings.n_markov)
        frame_found = margin >= settings.markov_margin
        note = f"frame {frame.frame_len}, {margin:.1f} sd above short-memory surrogates"
    else:
        note = f"best lag {frame.frame_len}"
    reading.tested("frame length", frame.p_bonferroni, frame_found, note)
    if frame_found:
        sync = codes.sync_word(stream, frame.frame_len)
        sentence = (f"The stream agrees with itself {frame.frame_len} bits later far more often than chance or "
                    f"short memory allows: it is built of units {frame.frame_len} bits long.")
        detail: dict[str, Any] = {"frame_len": int(frame.frame_len)}
        if sync.p_value <= settings.alpha_layout and sync.run_length >= 2:
            pattern = "".join(str(int(b)) for b in sync.pattern)
            sentence += (f" {sync.run_length} positions of every unit hold the same value each time (pattern "
                         f"{pattern}, starting at bit {sync.phase} of the unit): a fixed marker.")
            detail.update({"marker_bits": int(sync.run_length), "marker_phase": int(sync.phase), "marker": pattern})
        reading.add(ContentTier.LAYOUT, "frame length", frame.p_bonferroni, sentence, **detail)
        if not reading.layout:
            reading.layout = {"kind": "frames", **detail}

    # Rung 3d: parity checks. Skipped when a raster already accounts for the redundancy (and slow there).
    if not raster:
        found_codes = codes.block_code_scan(stream, n_max=settings.block_n_max, rng=rng_for(*rng_seed, "block"))
        code_found = False
        note = "none"
        if found_codes and found_codes[0].p_value <= settings.alpha_layout and found_codes[0].k >= 1:
            top = found_codes[0]
            strongest = int(np.argmax(top.satisfaction))
            margin = _markov_margin(lambda b: _check_satisfaction(b, top.n, top.offset, top.H[strongest]), bits,
                                    order, rng_for(*rng_seed, "block-markov"), settings.n_markov)
            code_found = margin >= settings.markov_margin
            note = f"({top.n},{top.k}) at offset {top.offset}, {margin:.1f} sd above short-memory surrogates"
            if code_found:
                reading.add(ContentTier.LAYOUT, "parity checks", top.p_value,
                            f"Cut into blocks of {top.n} bits from bit {top.offset}, the stream obeys "
                            f"{top.n - top.k} independent parity rules in nearly every block (the best in "
                            f"{100 * float(top.satisfaction.max()):.1f}% of blocks). Each block therefore "
                            f"carries at most {top.k} free bits out of {top.n}: fixed positions, or an "
                            f"error-detecting code.", n=int(top.n), k=int(top.k), offset=int(top.offset),
                            bit_error_rate=float(top.ber_estimate))
                if not reading.layout:
                    reading.layout = {"kind": "block code", "n": int(top.n), "k": int(top.k),
                                      "offset": int(top.offset)}
        reading.tested("parity checks", found_codes[0].p_value if found_codes else 1.0, code_found, note)

    # Rungs 4 and 5: anchors, from the fixed library.
    report = anchors.search(stream, rng_for(*rng_seed, "anchors"), alpha=settings.alpha_anchor)
    for hit, corrected in zip(report.hits, report.p_corrected):
        found = corrected <= settings.alpha_anchor
        rung = None
        sentence = hit.sentence
        if found and hit.anchor == "repeated_block":
            rung = ContentTier.LAYOUT
            if hit.evidence.get("agreement_with_consensus", 0.0) < 0.95:
                sentence = (f"The stream agrees with itself shifted by {hit.parameters['best_lag']} bits at "
                            f"{100 * hit.evidence['agreement']:.0f}% of positions, where a short-memory model "
                            f"predicts {100 * hit.evidence['agreement_null']:.0f}%: a unit of that length recurs "
                            "with changing content.")
        elif found and hit.anchor == "text_code":
            if _text_like(hit.evidence.get("preview", ""), hit.evidence.get("space_most_common", False), settings):
                rung = ContentTier.ANCHORED
                sentence += f' It begins: "{hit.evidence["preview"][:72]}".'
            else:
                rung = ContentTier.LAYOUT
                sentence += (" The characters do not form words, so this shows byte-level structure, not a "
                             "reading of the content.")
        elif found and hit.anchor == "functional_relations":
            rung = ContentTier.ANCHORED
            held = (hit.reading or {}).get("held_out") or {}
            if held.get("n", 0) >= 20 and held.get("correct", 0) / held["n"] >= settings.min_reading_accuracy:
                rung = ContentTier.HYPOTHESIS
                operators = ", ".join(f"{value} as '{op}'" for value, op in hit.reading["operators"].items())
                sentence += (f" As a hypothesis: the stream states arithmetic, with {operators}. Worked out from "
                             f"the first half of the stream alone, this reading predicts the result field of "
                             f"{held['correct']} of {held['n']} records in the second half.")
        elif found:
            rung = ContentTier.ANCHORED
        reading.tested(f"anchor: {hit.anchor}", corrected, rung is not None,
                       f"library of {report.n_anchors}; {hit.n_internal:,} settings inside this anchor")
        if rung is not None:
            reading.add(rung, f"anchor: {hit.anchor}", corrected, sentence, parameters=hit.parameters,
                        reading=hit.reading)


# ------------------------------------------------------------------------------------------------------------
# Other payloads
# ------------------------------------------------------------------------------------------------------------

def _analyse_symbols(seq: SymbolSequence, reading: Reading, rng_seed: tuple, settings: Settings) -> None:
    symbols = np.asarray(seq.symbols)
    k = int(seq.alphabet_size)
    if symbols.size < settings.min_bits:
        reading.notes.append(f"Only {symbols.size} symbols: too few for the structure tests.")
        return
    unit = _periodic_unit(symbols, k)
    if unit is not None:
        period, explained, pattern = unit
        reading.notes.append(f"The sequence repeats every {period} symbols ({100 * explained:.1f}% of symbols are "
                             "explained by one repeating cycle). A repeating cycle is a clock: it carries no "
                             "information beyond the cycle itself.")
        reading.layout = {"kind": "repeating cycle", "period": period, "explained": explained}
        changes = int(np.sum(pattern != np.roll(pattern, 1)))
        if period < settings.min_bits or changes <= 4:
            reading.notes.append(f"The cycle has {changes} changes of state: nothing further to analyse.")
            return
        symbols = pattern
    if k == 2:
        _analyse_bits(BitStream(symbols.astype(np.uint8)), reading, rng_seed, settings)
        return
    result = battery(symbols, k, rng_for(*rng_seed, "battery"), max_null_order=1,
                     n_surrogates=settings.n_surrogates, alpha=settings.alpha_structure)
    reading.tested("structure battery", result.p_value, result.detected, f"null: {result.null}")
    if result.detected:
        reading.add(ContentTier.STRUCTURE, "structure battery", result.p_value,
                    f"The sequence of {k} symbols has order that a {result.null} chain fitted to it does not "
                    f"produce (probability {result.p_value:.3f}).", statistics=result.statistics)


def _analyse_series(series: Series, reading: Reading, rng_seed: tuple, settings: Settings) -> None:
    levels = symbolise.discreteness_test(series, rng_for(*rng_seed, "levels"))
    reading.tested("discrete levels", levels.p_value, levels.detected, levels.note)
    if not levels.detected:
        reading.notes.append("The values do not sit on separate levels: the series is continuous as far as its "
                             "error bars can tell. Periodicity in a continuous series is a question for the "
                             "periodicity search, not for this workbench.")
        return
    reading.add(ContentTier.SYMBOLS, "discrete levels", levels.p_value,
                f"All but {100 * levels.off_level_fraction:.1f}% of the points lie within their error bars of one "
                f"of {levels.k} levels, which are at least {levels.separation:.0f} error bars apart; the series "
                f"changes level {levels.transitions} times.", levels=[float(v) for v in levels.levels],
                occupancy=[float(v) for v in levels.occupancy])
    gaps = np.diff(series.time)
    if gaps.size and np.max(gaps) > 1.5 * np.median(gaps):
        reading.notes.append("The sampling is uneven, so the order of states was not analysed as a sequence.")
        return
    _analyse_symbols(levels.symbols, reading, rng_seed, settings)


def _analyse_events(events: EventList, reading: Reading, rng_seed: tuple, settings: Settings) -> None:
    from ..detect import intervals

    clock = intervals.clock_test(events, rng_for(*rng_seed, "clock"))
    reading.tested("clock", _p(clock.p_value), clock.detected, f"period {clock.period:.6g}")
    if clock.detected:
        reading.add(ContentTier.SYMBOLS, "clock", _p(clock.p_value),
                    f"The events fall on the ticks of a clock of period {clock.period:.6g} {events.time_unit}, "
                    f"scattering by {100 * clock.residual_rms:.1f}% of a period; {100 * clock.occupancy:.0f}% of "
                    "ticks carry an event. Each tick is therefore a slot that is either full or empty.",
                    period=float(clock.period), occupancy=float(clock.occupancy))
        slots = intervals.to_bits(events, clock.period, clock.phase)
        reading.layout = {"kind": "clocked slots", "period": float(clock.period), "slots": len(slots)}
        _analyse_symbols(SymbolSequence(np.asarray(slots.bits, dtype=int), 2), reading, rng_seed, settings)
        return
    spacing = intervals.spacing_test(events, rng_for(*rng_seed, "spacing"))
    reading.tested("interval quantum", _p(spacing.p_value), spacing.detected, f"quantum {spacing.quantum:.6g}")
    if spacing.detected:
        reading.add(ContentTier.SYMBOLS, "interval quantum", _p(spacing.p_value),
                    f"The gaps between events are whole multiples of {spacing.quantum:.6g} {events.time_unit}.",
                    quantum=float(spacing.quantum))
    else:
        reading.notes.append("The event times show no clock and no common unit of spacing.")


# ------------------------------------------------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------------------------------------------------

_STATEMENTS = {
    0: ("No structure was found. By these tests the signal cannot be told from its null model. That is not "
        "proof that it carries nothing: well-compressed or encrypted content looks exactly like this."),
    1: ("The signal takes discrete states, so it can be written down as a sequence of symbols. The sequence "
        "shows no order beyond what a short-memory or repeating natural process produces."),
    2: ("The sequence has order beyond a short-memory process, but no layout, code or anchor was recovered. "
        "What the order is remains open."),
    3: ("A layout was recovered: how the signal is organised is known. What it says is not. Nothing in it "
        "matched the fixed library of shared knowledge."),
    4: ("Part of the content matches shared knowledge from the fixed library, with the probability of a chance "
        "match corrected for everything in that library."),
    5: ("A specific reading is proposed and has predicted held-out parts of the signal. It is a hypothesis, "
        "open to refutation by more of the same signal."),
}
_LIMITS = ("Tests made: discrete levels and clocks; a battery of structure statistics against a fitted "
           "short-memory chain; raster and volume layouts; frame lengths and fixed markers; linear parity "
           "checks up to the block length set; linear shift registers; and a fixed library of six anchors. "
           "Not seen by any of them: non-linear codes, variable-length frames, compression, encryption, and "
           "any convention outside the library.")


def analyse(payload: BitStream | SymbolSequence | Series | EventList, *, name: str = "", seed: int = 0,
            settings: Settings = Settings()) -> Reading:
    """Climb the content ladder for one payload and report where the evidence stops."""
    kind = type(payload).__name__
    reading = Reading(name or getattr(payload, "source_id", "") or kind, kind, len(payload))
    rng_seed = (seed, "workbench", reading.name)
    if isinstance(payload, BitStream):
        _analyse_bits(payload, reading, rng_seed, settings)
    elif isinstance(payload, SymbolSequence):
        _analyse_symbols(payload, reading, rng_seed, settings)
    elif isinstance(payload, Series):
        _analyse_series(payload, reading, rng_seed, settings)
    elif isinstance(payload, EventList):
        _analyse_events(payload, reading, rng_seed, settings)
    else:
        raise TypeError(f"cannot analyse a {kind}")
    reading.tier = max((e.rung for e in reading.evidence), default=0)
    reading.tier_phrase = content_phrase(reading.tier)
    reading.statement = _STATEMENTS[reading.tier]
    reading.limits = _LIMITS
    return reading


def to_markdown(reading: Reading) -> str:
    """A reading as Markdown: conclusion, evidence, picture, what was tried, limits."""
    lines = [f"# Reading of {reading.name}", "",
             f"**Content tier {reading.tier} of 5: {reading.tier_phrase}.** {reading.statement}", "",
             f"Input: {reading.kind}, {reading.n:,} elements.", ""]
    if reading.evidence:
        lines += ["## Evidence", ""]
        for item in sorted(reading.evidence, key=lambda e: e.rung):
            lines.append(f"- (tier {item.rung}, {item.test}, p = {item.p_value:.2g}) {item.sentence}")
        lines.append("")
    if reading.notes:
        lines += ["## Notes", ""] + [f"- {note}" for note in reading.notes] + [""]
    if reading.picture:
        lines += ["## The raster", "", "```", reading.picture, "```", ""]
    if reading.code_lengths:
        lines += ["## Bits needed to write the stream down", ""]
        lines += [f"- {model}: {bits:,.0f}" for model, bits in reading.code_lengths.items()] + [""]
    lines += ["## Everything that was tried", "", "| Test | Probability | Found | Note |", "|---|---|---|---|"]
    lines += [f"| {t['test']} | {t['p_value']:.2g} | {'yes' if t['found'] else 'no'} | {t['note']} |"
              for t in reading.tried]
    lines += ["", "## Limits", "", reading.limits, ""]
    return "\n".join(lines)
