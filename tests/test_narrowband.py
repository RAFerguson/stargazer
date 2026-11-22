"""The narrowband drift search: normalisation, recovery against drift and S/N, calibration on noise, the cadence
filter, and the filterbank reader and writer. Run with -s to see the measured numbers."""
import time

import numpy as np
import pytest
from scipy import stats

from stargazer.bench import spectra
from stargazer.core import rng_for
from stargazer.detect import narrowband as nb
from stargazer.sources import filterbank

DEFAULT = nb.DriftConfig()
NARROW = nb.DriftConfig(max_drift_hz_per_s=1.0)


def _match(hits, signal, spec, drift_step):
    """The hit that recovers a planted tone: start within 2 channels and drift within 3 trial steps."""
    for h in hits:
        if abs(h.channel - signal["channel"]) <= 2.0 and abs(h.drift_hz_per_s - signal["drift_hz_per_s"]) <= \
                3 * drift_step + 1e-9:
            return h
    return None


def _recovery(drifts, snr, n_spec, config, seed, n_freq=4096, tones_per_spec=4):
    """Fraction of planted tones recovered above threshold at each |drift|, and the recovered S/N ratios."""
    out = {}
    for drift in drifts:
        found, ratios = 0, []
        for k in range(n_spec):
            rng = rng_for(seed, "recovery", drift, k)
            signs = rng.choice([-1.0, 1.0], size=tones_per_spec)
            spec, _ = spectra.noise(rng, 16, n_freq)
            sweep = abs(drift) * spec.duration / abs(spec.f_step)
            slots = np.linspace(sweep + 10, n_freq - sweep - 10, tones_per_spec)
            tones = [{"channel": float(c + rng.uniform(-0.5, 0.5)), "drift_hz_per_s": float(s * drift), "snr": snr}
                     for c, s in zip(slots, signs)]
            spec, truth = spectra.drifting_tones(rng, tones, n_freq=n_freq, n_time=16)
            result = nb.dedoppler(spec, config)
            hits = nb.find_hits(result, config)
            for signal in truth["signals"]:
                h = _match(hits, signal, spec, result.drift_step_hz_per_s)
                found += h is not None
                ratios.append(h.snr / snr if h else 0.0)
        out[drift] = (found / (n_spec * tones_per_spec), float(np.median(ratios)))
    return out


# ------------------------------------------------------------------------------------------------ basics

def test_spectrogram_checks_its_inputs():
    with pytest.raises(ValueError):
        nb.Spectrogram(np.ones(10), 1.0, 0.0, 1.0)
    with pytest.raises(ValueError):
        nb.Spectrogram(-np.ones((2, 10)), 1.0, 0.0, 1.0)
    with pytest.raises(ValueError):
        nb.Spectrogram(np.ones((2, 10)), 1.0, 0.0, 0.0)
    spec = nb.Spectrogram(np.ones((4, 10)), 2.0, 100.0, -3.0)
    assert spec.duration == 8.0 and spec.frequencies()[-1] == pytest.approx(73.0)


def test_drift_grid_includes_zero_and_spans_the_range():
    spec, _ = spectra.noise(rng_for(1, "grid"), 16, 256)
    hz, ch = nb.drift_grid(spec, DEFAULT)
    assert np.any(hz == 0) and hz[0] >= -4.0 - 1e-9 and hz[-1] <= 4.0 + 1e-9 and hz[-1] > 3.99
    assert np.all(np.diff(hz) > 0)
    step = abs(spec.f_step) / (16 * spec.t_samp)
    assert np.allclose(np.diff(hz), step)
    assert np.allclose(ch * spec.f_step / spec.t_samp, hz)


def test_boxcar_widens_with_drift_and_naive_takes_one_channel():
    lo, hi = nb.track_ranges(0.0, 0.0, 16)
    assert np.all(lo == 0) and np.all(hi == 0)
    lo, hi = nb.track_ranges(5.0, 0.0, 16)
    assert np.all(hi - lo + 1 >= 5) and np.all(hi - lo + 1 <= 6)
    assert lo[1] - lo[0] == 5                          # the run moves five channels per integration
    lo, hi = nb.track_ranges(-5.0, 0.5, 16, boxcar=False)
    assert np.all(lo == hi) and hi[0] == -2


def test_estimated_n_avg_and_exact_gamma_tail():
    # few raw samples per cell: the de-drifted sum is visibly skewed, and the gamma tail must be used
    rng = rng_for(2, "tail")
    spec, _ = spectra.noise(rng, 16, 4096, n_avg=1, record_n_avg=False)
    # one start offset, so that every plane cell is one trial (with two, it is the larger of two correlated ones)
    config = nb.DriftConfig(max_drift_hz_per_s=0.2, keep_plane_cells=10 ** 8, start_offsets=1)
    result = nb.dedoppler(spec, config)
    assert result.noise.n_avg_source != "meta" and result.noise.n_avg == pytest.approx(1.0, rel=0.1)
    j = np.repeat(np.arange(result.plane.shape[0])[:, None], spec.n_freq, axis=1)
    finite = np.isfinite(result.plane)
    p = result.p_local(result.plane[finite].astype(float), j[finite])
    gauss = stats.norm.sf(result.plane[finite].astype(float))
    print(f"\nn_avg=1, {p.size} trials: P(p_local<=0.01)={np.mean(p <= 0.01):.4f}, "
          f"<=0.001: {np.mean(p <= 1e-3):.5f};"
          f" Gaussian tail would give {np.mean(gauss <= 0.01):.4f} and {np.mean(gauss <= 1e-3):.5f}")
    assert 0.007 <= np.mean(p <= 0.01) <= 0.013
    assert 0.0005 <= np.mean(p <= 1e-3) <= 0.0016
    assert np.mean(gauss <= 1e-3) > 2.0 * 1e-3                 # a Gaussian would overstate significance


def test_snr_spread_is_one_at_every_drift_with_a_bandpass():
    spec, _ = spectra.noise(rng_for(3, "spread"), 16, 16384)
    result = nb.dedoppler(spec, DEFAULT)
    print(f"\nS/N spread over channels: {result.snr_scale.min():.3f} to {result.snr_scale.max():.3f} "
          f"(inflation a = {result.inflation:.2e})")
    assert np.all(np.abs(result.snr_scale - 1.0) < 0.12)
    assert np.median(np.abs(result.snr_scale - 1.0)) < 0.05


# ---------------------------------------------------------------------------------------- noise only

def test_false_alarms_on_noise():
    worst, any_hit = [], 0
    for k in range(100):
        spec, _ = spectra.noise(rng_for(4, "false-alarm", k), 16, 4096)
        result = nb.dedoppler(spec, NARROW)
        any_hit += bool(nb.find_hits(result, NARROW))
        worst.append(float(result.snr_max.max()))
    print(f"\n100 noise spectrograms (16 x 4096, +-1 Hz/s, {result.n_trials} trials each): {any_hit} with a hit "
          f">= 10; largest S/N median {np.median(worst):.2f}, max {max(worst):.2f}")
    assert any_hit == 0 and max(worst) < 8.0


def test_decoy_probability_is_not_small_on_noise_and_tiny_for_a_tone():
    config = nb.DriftConfig(max_drift_hz_per_s=0.5, n_decoys=40)
    ps = []
    for k in range(8):
        spec, _ = spectra.noise(rng_for(5, "decoy-noise", k), 16, 512)
        ps.append(nb.search(spec, rng_for(5, "decoys", k), config).p)
    print(f"\n8 noise spectrograms: decoy probabilities {np.round(sorted(ps), 3).tolist()}")
    assert sum(p <= 0.05 for p in ps) <= 2 and np.median(ps) > 0.15
    spec, _ = spectra.drifting_tones(rng_for(5, "tone"), [{"channel": 200.3, "drift_hz_per_s": 0.2, "snr": 15}],
                                     n_freq=512)
    result = nb.search(spec, rng_for(5, "tone-decoys"), config)
    assert result.p == pytest.approx(1 / 41) and result.decoy.exceedances == 0
    again = nb.search(spec, rng_for(5, "tone-decoys"), config)
    assert again.p == result.p and again.statistic == result.statistic


# -------------------------------------------------------------------------------------------- recovery

def test_recovery_holds_at_high_drift_with_the_boxcar_and_falls_without():
    drifts = [0.0, 0.3, 1.0, 2.0, 3.9]
    box = _recovery(drifts, 15.0, 2, DEFAULT, 6)
    naive = _recovery([0.0, 1.0, 3.9], 15.0, 1, nb.DriftConfig(boxcar=False), 6)
    print("\n|drift| Hz/s: boxcar recovered (median S/N ratio) | naive")
    for d in drifts:
        other = f"{naive[d][0]:.2f} ({naive[d][1]:.2f})" if d in naive else "-"
        print(f"  {d:4.2f}: {box[d][0]:.2f} ({box[d][1]:.2f}) | {other}")
    assert np.mean([box[d][0] for d in drifts]) >= 0.9
    assert all(box[d][0] >= 0.75 for d in drifts)
    assert naive[1.0][0] <= 0.25 and naive[3.9][0] <= 0.25 and naive[0.0][0] >= 0.75


def test_frequency_and_drift_accuracy():
    rng = rng_for(7, "accuracy")
    spec, truth = spectra.drifting_tones(rng, n_tones=4, snr=30.0)
    result = nb.dedoppler(spec, DEFAULT)
    hits = nb.find_hits(result, DEFAULT)
    assert len(hits) == 4                                       # one hit per tone: no butterfly left over
    for signal in truth["signals"]:
        h = min(hits, key=lambda h: abs(h.channel - signal["channel"]))
        assert abs(h.frequency_hz - signal["f0_hz"]) <= 1.0 * abs(spec.f_step)
        assert abs(h.drift_hz_per_s - signal["drift_hz_per_s"]) <= 2 * result.drift_step_hz_per_s
        assert h.snr == pytest.approx(30.0, rel=0.25) and h.p_local < 1e-20


def test_wide_tone_and_band_edges():
    rng = rng_for(8, "edges")
    spec, _ = spectra.noise(rng, 16, 2048)
    near_edge = spec.f_start + 3.0 * spec.f_step
    spec = nb.inject_tone(spec, near_edge, 1.0 * spec.f_step / spec.t_samp, 20.0)   # leaves the band
    hits = nb.find_hits(spec, DEFAULT)
    assert hits and abs(hits[0].channel - 3.0) <= 2 and hits[0].snr > 10


# ---------------------------------------------------------------------------------------- interference

def test_broadband_pulse_is_divided_out():
    spec, truth = spectra.broadband_pulse(rng_for(9, "pulse"))
    result = nb.dedoppler(spec, DEFAULT)
    assert not nb.find_hits(result, DEFAULT)
    assert result.noise.level[truth["pulse_row"]].min() > 5 * np.median(result.noise.level)


def test_wandering_interferer_is_not_a_line():
    spec, truth = spectra.wandering_interferer(rng_for(10, "wander"), snr=40.0)
    hits = nb.find_hits(spec, DEFAULT)
    best = max((h.snr for h in hits), default=0.0)
    print(f"\nwandering interferer, ideal S/N 40: strongest linear hit {best:.1f} from {len(hits)} hits")
    assert best < 0.75 * truth["signals"][0]["snr"]


def test_cadence_keeps_the_candidate_and_rejects_interference():
    scans, truth = spectra.cadence(rng_for(11, "cadence"))
    result = nb.cadence_search(scans, truth["is_on"], NARROW)
    signals = {s["name"]: s for s in truth["signals"]}
    step = abs(scans[0].f_step)

    def group(name):
        return min(result.candidates, key=lambda c: abs(c.frequency_hz - signals[name]["f0_hz"]))

    # the constant interferer is a real narrowband signal: found in every ON scan, rejected by the OFF scans
    constant = group("constant")
    assert abs(constant.frequency_hz - signals["constant"]["f0_hz"]) < 3 * step
    assert all(h is not None for h in constant.on_hits) and not constant.accepted and "OFF" in constant.reason
    hidden = group("hidden")
    assert not hidden.accepted and all(h.snr < 10 for _, h in hidden.off_hits)
    accepted = result.accepted
    assert len(accepted) == 1 and abs(accepted[0].frequency_hz - signals["candidate"]["f0_hz"]) < 3 * step
    assert accepted[0].drift_hz_per_s == pytest.approx(signals["candidate"]["drift_hz_per_s"], abs=0.005)
    assert result.swapped_accepted == 0
    # searching the OFF scans at the ON threshold lets the hidden interferer through: the BLC1 lesson
    lax = nb.cadence_search(scans, truth["is_on"], nb.DriftConfig(max_drift_hz_per_s=1.0, off_snr_threshold=10.0))
    assert len(lax.accepted) == 2


def test_cadence_filter_reasons():
    def hit(f, d, snr=20.0):
        return nb.Hit(f, d, snr, 0.0, 0, 1.0, 1e-30)

    on = [[hit(1000.0, 0.5), hit(5000.0, 0.0)], [hit(1300.0, 0.5), hit(5000.0, 0.0)], [hit(5000.0, 0.0)]]
    off = [[], [hit(5000.0, 0.0, 6.0)]]
    out = nb.cadence_filter(on, off, [0.0, 600.0, 1200.0], [300.0, 900.0], 300.0, freq_tolerance_hz=10.0,
                            drift_tolerance_hz_per_s=0.05)
    by_freq = {round(c.frequency_hz): c for c in out}
    assert not by_freq[1000].accepted and "absent from ON scan 3" in by_freq[1000].reason
    assert not by_freq[5000].accepted and "OFF scan 2" in by_freq[5000].reason
    clean = nb.cadence_filter(on, [[], []], [0.0, 600.0, 1200.0], [300.0, 900.0], 300.0, freq_tolerance_hz=10.0,
                              drift_tolerance_hz_per_s=0.05)
    zero = [c for c in clean if c.accepted]
    assert len(zero) == 1 and zero[0].flags and "zero drift" in zero[0].flags[0]


# ------------------------------------------------------------------------------------------- filterbank

@pytest.mark.parametrize("nbits", [8, 16, 32])
def test_filterbank_round_trip(tmp_path, nbits):
    rng = rng_for(12, "fil", nbits)
    spec, _ = spectra.noise(rng, 8, 300, level=100.0 if nbits == 8 else 3000.0)
    if nbits < 32:
        spec = spec.with_power(np.round(spec.power).clip(0, 2 ** nbits - 1))
    spec.meta.update({"src_raj": 123456.7, "telescope_id": 6})
    spec.source_id = "HIP 12345"
    path = filterbank.write(spec, tmp_path / f"x{nbits}.fil", nbits=nbits)
    back = filterbank.read(path)
    assert back.power.shape == spec.power.shape
    if nbits == 32:
        assert np.allclose(back.power, spec.power, rtol=1e-6)
    else:
        assert np.array_equal(back.power, spec.power)
    assert back.t_samp == spec.t_samp and back.f_start == pytest.approx(spec.f_start, abs=1e-6)
    assert back.f_step == pytest.approx(spec.f_step, rel=1e-12)
    assert back.source_id == "HIP 12345" and back.meta["src_raj"] == 123456.7 and back.meta["telescope_id"] == 6
    assert back.meta["tstart"] == spec.meta["tstart"] and back.meta["nbits"] == nbits


def test_filterbank_slicing_and_refusals(tmp_path):
    spec, _ = spectra.noise(rng_for(13, "slice"), 10, 500)
    path = filterbank.write(spec, tmp_path / "s.fil")
    freqs = spec.frequencies()
    lo, hi = sorted((freqs[100], freqs[199]))
    part = filterbank.read(path, f_range=(lo, hi), t_range=(3 * spec.t_samp, 6.5 * spec.t_samp))
    assert part.power.shape == (4, 100)
    assert np.allclose(part.power, spec.power[3:7, 100:200], rtol=1e-6)
    assert part.f_start == pytest.approx(freqs[100], abs=1e-3) and part.meta["channel_offset"] == 100
    assert part.meta["tstart"] == pytest.approx(spec.meta["tstart"] + 3 * spec.t_samp / 86400.0)
    with pytest.raises(ValueError):
        filterbank.write(spec, tmp_path / "bad.fil", nbits=8)            # not integers
    with pytest.raises(ValueError):
        filterbank.write(spec, tmp_path / "bad.fil", nbits=4)
    (tmp_path / "junk.fil").write_bytes(b"\x00" * 64)
    with pytest.raises(filterbank.FilterbankError):
        filterbank.read(tmp_path / "junk.fil")


def test_search_reads_a_written_file(tmp_path):
    spec, truth = spectra.drifting_tones(rng_for(14, "file-search"), [{"channel": 700.2, "drift_hz_per_s": -0.6,
                                                                        "snr": 25.0}], n_freq=2048)
    back = filterbank.read(filterbank.write(spec, tmp_path / "t.fil"))
    hits = nb.find_hits(back, DEFAULT)            # meta has no n_avg after the round trip: it is estimated
    assert hits and abs(hits[0].channel - 700.2) < 2 and hits[0].drift_hz_per_s == pytest.approx(-0.6, abs=0.03)


# ------------------------------------------------------------------------------------------------ slow

@pytest.mark.slow
def test_slow_false_alarms_full_drift_range():
    any_hit, worst = 0, []
    for k in range(100):
        spec, _ = spectra.noise(rng_for(20, "fa-full", k), 16, 4096)
        result = nb.dedoppler(spec, DEFAULT)
        any_hit += bool(nb.find_hits(result, DEFAULT))
        worst.append(float(result.snr_max.max()))
    print(f"\n100 noise spectrograms (16 x 4096, +-4 Hz/s, {result.n_trials} trials each): {any_hit} with a hit; "
          f"largest S/N median {np.median(worst):.2f}, max {max(worst):.2f}")
    assert any_hit == 0


@pytest.mark.slow
def test_slow_decoy_calibration():
    config = nb.DriftConfig(max_drift_hz_per_s=1.0, n_decoys=100)
    ps = []
    for k in range(100):
        spec, _ = spectra.noise(rng_for(21, "calibration", k), 16, 4096)
        ps.append(nb.search(spec, rng_for(21, "calibration-decoys", k), config).p)
    ps = np.array(ps)
    print(f"\n100 noise spectrograms (16 x 4096, +-1 Hz/s, <= 100 decoys): fraction p <= 0.05 = "
          f"{np.mean(ps <= 0.05):.2f}, p <= 0.1 = {np.mean(ps <= 0.1):.2f}, p <= 0.5 = {np.mean(ps <= 0.5):.2f}")
    assert np.mean(ps <= 0.05) <= 0.12 and np.mean(ps <= 0.5) <= 0.65


@pytest.mark.slow
def test_slow_recovery_curves():
    drifts = [0.0, 0.08, 0.15, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 3.9]
    print("\nrecovery at S/N 15, 16 x 4096, 40 tones per point: boxcar (median S/N ratio) | naive")
    box = _recovery(drifts, 15.0, 10, DEFAULT, 22)
    naive = _recovery(drifts, 15.0, 10, nb.DriftConfig(boxcar=False), 22)
    for d in drifts:
        print(f"  {d:4.2f} Hz/s: {box[d][0]:.3f} ({box[d][1]:.2f}) | {naive[d][0]:.3f} ({naive[d][1]:.2f})")
    print("recovery against S/N with the boxcar, 40 tones per point, at 0.3 and 3.0 Hz/s:")
    for snr in (8.0, 10.0, 12.0, 15.0, 20.0):
        curve = _recovery([0.3, 3.0], snr, 10, DEFAULT, 23)
        print(f"  S/N {snr:4.1f}: {curve[0.3][0]:.3f} | {curve[3.0][0]:.3f}")
    assert all(box[d][0] >= 0.9 for d in drifts)


@pytest.mark.slow
def test_slow_timing_65536_channels():
    spec, _ = spectra.noise(rng_for(24, "timing"), 16, 65536)
    start = time.perf_counter()
    result = nb.dedoppler(spec, DEFAULT)
    elapsed = time.perf_counter() - start
    print(f"\n16 x 65536, +-4 Hz/s ({result.drift_hz_per_s.size} drifts, {result.n_trials} trials): "
          f"{elapsed:.1f} s")
    assert elapsed < 30
