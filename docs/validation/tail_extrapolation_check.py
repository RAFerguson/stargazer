"""How good is the fitted-tail extrapolation? Compare it with brute-force simulation.
For one realistic sampling pattern: draw many null statistics, fit the tail to 1,000 of them (as the pipeline
does), and compare the predicted tail probability with the empirical one at rarer levels."""
import sys, time, json
import numpy as np
from stargazer.core import rng_for
from stargazer.bench.natural import survey_times
from stargazer.detect import periodicity as P
from stargazer.nulls import drw

n_total = int(sys.argv[1])
rng = rng_for(3, "tail-check")
t = np.unique(np.floor(np.concatenate([survey_times(rng), survey_times(rng)]))) + 0.3
e = np.full(t.size, 0.03)
cfg = P.SearchConfig()
freq = P.frequency_grid(t[-1] - t[0], cfg)
sa, ta = drw.default_grids(t, 0.15 * np.ones(t.size) * rng.standard_normal(t.size), e, n_sigma=cfg.n_sigma, n_tau=cfg.n_tau)
stats = []
t0 = time.time()
while len(stats) * 50 < n_total:
    sigma = np.full(50, 0.15); tau = np.full(50, 300.0)
    fake = drw.simulate(t, e, sigma, tau, rng)
    stats.append(P.evaluate(t, e, fake, freq, sa, ta, cfg).statistic)
s = np.concatenate(stats)
print(f"{s.size} null statistics in {time.time()-t0:.0f}s; N={t.size}, F={freq.size}")
out = {"n": int(s.size)}
srt = np.sort(s)[::-1]
for level in (1e-2, 1e-3, 3e-4, 1e-4):
    k = int(round(level * s.size))
    if k < 3: continue
    threshold = srt[k - 1]
    preds = []
    for j in range(20):
        sub = rng.choice(s, 1000, replace=False)
        preds.append(P.tail_probability(float(threshold), sub)[0])
    preds = np.array(preds)
    out[str(level)] = {"threshold": float(threshold), "predicted_median": float(np.median(preds)), "predicted_min": float(preds.min()), "predicted_max": float(preds.max())}
    print(f"empirical {level:.0e} (statistic {threshold:.1f}): fitted tail predicts median {np.median(preds):.1e}, range {preds.min():.1e} to {preds.max():.1e}; ratio {np.median(preds)/level:.2f}")
json.dump(out, open("tail_check.json", "w"))
