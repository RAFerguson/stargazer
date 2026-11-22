"""Error control: sequential probabilities, online false-discovery control, sealing, empirical null, injections."""
import numpy as np
import pytest
from scipy.stats import norm

from stargazer.core import Series
from stargazer.trust import blind, empirical_null, inject
from stargazer.trust.decoys import benjamini_hochberg, sequential_p, uniformity
from stargazer.trust.onlinefdr import OnlineFDR, gamma


def _sequential(p0, rng, **kw):
    return sequential_p(1.0, lambda n: 2.0 * (rng.random(n) < p0), **kw)


def test_sequential_probabilities_are_calibrated():
    rng = np.random.default_rng(0)
    p = np.array([_sequential(rng.random(), rng, max_draws=400).p for _ in range(3000)])
    check = uniformity(p)
    assert abs(check["z_below_0.05"]) < 3.5 and abs(check["z_below_0.01"]) < 3.5
    assert p.min() >= 1 / 401


def test_sequential_stops_early_for_unremarkable_objects():
    rng = np.random.default_rng(1)
    result = _sequential(0.5, rng, max_draws=2000)
    assert result.stopped_early and result.draws < 60
    strong = _sequential(0.0, rng, max_draws=300)
    assert not strong.stopped_early and strong.p == pytest.approx(1 / 301)


def test_benjamini_hochberg():
    p = np.array([1e-5, 0.01, 0.02, 0.5, 0.9])
    rejected, adjusted = benjamini_hochberg(p, q=0.05)
    assert rejected.tolist() == [True, True, True, False, False]
    assert np.all(adjusted >= p) and np.all(adjusted <= 1)


def test_gamma_sequence_spends_no_more_than_the_budget():
    assert 0.5 < gamma(np.arange(1, 2_000_001)).sum() < 1.0


def test_online_fdr_controls_false_discoveries_and_can_be_rebuilt():
    rng = np.random.default_rng(2)
    false, total = 0, 0
    for _ in range(150):
        state, decisions = OnlineFDR(0.05), []
        for _ in range(200):
            real = rng.random() < 0.1
            p = float(np.exp(-rng.exponential(12))) if real else rng.random()
            reject, level = state.test(p)
            decisions.append(reject)
            total += reject
            false += reject and not real
        rebuilt = OnlineFDR.from_history(decisions, 0.05)
        assert rebuilt.next_level() == pytest.approx(state.next_level())
    assert total > 500 and false / total < 0.05


def test_online_level_never_depends_on_the_next_p():
    state = OnlineFDR(0.05)
    level = state.next_level()
    assert state.test(0.9) == (False, level)
    with pytest.raises(ValueError):
        state.test(1.5)


def _subject(name, kind="real", parameters=None):
    time = np.arange(50.0)
    series = Series(time, np.sin(time), np.full(50, 0.1), "g", f"survey:{name}", meta={"objectid": 12345})
    return blind.Subject(name, [series], kind, parameters or {})


def test_sealing_hides_kinds_and_can_be_verified():
    subjects = [_subject("a"), _subject("a", "decoy"), _subject("b"), _subject("b", "plant", {"amplitude": 0.1})]
    sealed = blind.seal(subjects, np.random.default_rng(3))
    assert len(sealed.subjects) == 4 and {s.kind for s in sealed.subjects} == {"sealed"}
    for subject in sealed.subjects:
        assert subject.parameters == {} and len(subject.name) == 16
        assert all(band.source_id == subject.name and band.meta == {} and band.provenance is None
                   for band in subject.bands)
    key = sealed.unseal()
    assert sorted(v["kind"] for v in key.values()) == ["decoy", "plant", "real", "real"]
    proof = sealed.proof()
    assert blind.SealedSet.verify(proof)
    proof["key"][next(iter(proof["key"]))]["kind"] = "real!"
    assert not blind.SealedSet.verify(proof)


def test_decoy_estimate():
    out = blind.decoy_estimate(3, 1, 1000, 1000)
    assert out["expected_false_among_real"] == 1.0 and out["estimated_fdr"] == pytest.approx(2 / 3)
    assert 0 < out["false_alarm_probability_upper95"] < 0.01


def test_empirical_null_recovers_centre_and_width():
    rng = np.random.default_rng(4)
    theoretical = empirical_null.fit(rng.random(4000))
    assert abs(theoretical.centre) < 0.08 and theoretical.width == pytest.approx(1.0, abs=0.06)
    wide = empirical_null.fit(norm.sf(0.2 + 1.3 * rng.standard_normal(4000)))
    assert wide.width == pytest.approx(1.3, abs=0.1) and wide.centre == pytest.approx(0.2, abs=0.12)
    p = np.array([1e-6, 1e-3, 0.3])
    assert np.all(wide.recalibrate(p) > p[:3] * 0.99) and np.all(np.diff(wide.recalibrate(p)) > 0)
    narrow = empirical_null.fit(norm.sf(0.7 * rng.standard_normal(4000)))
    assert narrow.width_used == 1.0                       # never allowed to make results more significant
    assert empirical_null.fit(rng.random(50)).width == 1.0


def test_injections():
    time = np.linspace(0, 1000, 500)
    base = Series(time, np.zeros(500), np.full(500, 0.01))
    assert np.ptp(inject.sinusoid(base, 100.0, 0.2).value) == pytest.approx(0.4, rel=0.01)
    assert np.ptp(inject.sawtooth(base, 100.0, 0.2).value) == pytest.approx(0.4, rel=0.05)
    flares = inject.periodic_flares(base, 100.0, 0.3, 0.1).value
    assert flares.max() == pytest.approx(0.3, rel=0.05) and np.mean(flares > 0.15) == pytest.approx(0.1, abs=0.03)
    keyed = inject.on_off_keying(base, np.array([1, 0, 1, 1]), 50.0, 0.1).value
    assert set(np.round(keyed, 6)) == {0.0, 0.1} and base.value.max() == 0.0
    cells = inject.recovery_table([base] * 4, inject.sinusoid, [{"period_days": 50.0, "amplitude": 0.1}],
                                  lambda series, parameters: np.ptp(series.value) > 0.1)
    assert cells[0].completeness == 1.0 and cells[0].interval()[0] > 0.3
