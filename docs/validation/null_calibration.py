"""Null calibration for two-band objects through merge_bands + search. usage: calib3.py n [mode]"""
import sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from stargazer.core import Series, rng_for
from stargazer.bench.natural import survey_times, simulate_drw
from stargazer.detect.periodicity import SearchConfig, search, merge_bands
from stargazer.trust.decoys import uniformity

CFG = SearchConfig(max_decoys=400)
MODE = sys.argv[2] if len(sys.argv) > 2 else "two"

def make(i):
    rng = rng_for(77, "obj", i)
    sigma, tau = rng.uniform(0.08, 0.25), 10 ** rng.uniform(1.7, 3.0)
    tg, tr = survey_times(rng, nights=(25, 40)), survey_times(rng, nights=(30, 45))
    t_all = np.sort(np.concatenate([tg, tr]))
    x = simulate_drw(t_all, sigma, tau, rng)
    bands = []
    for name, t, err, scale in (("g", tg, 0.04, 1.0), ("r", tr, 0.03, 0.85)):
        idx = np.searchsorted(t_all, t)
        e = err * rng.uniform(0.8, 1.3, t.size)
        bands.append(Series(t, 18.5 + scale * x[idx] + e * rng.standard_normal(t.size), e, name, f"sim:{i}:{name}"))
    if MODE == "one":
        return bands[0]
    return merge_bands(bands)

def run(i):
    r = search(make(i), rng_for(9, "search3", i), CFG)
    return r.p, r.period_days, r.statistic, r.noise["tau_null_days"]

if __name__ == "__main__":
    n = int(sys.argv[1])
    t0 = time.time()
    with ProcessPoolExecutor(4) as pool:
        out = np.array(list(pool.map(run, range(n), chunksize=4)))
    print(MODE, f"{n} objects in {time.time()-t0:.0f}s", uniformity(out[:, 0]))
    hits = out[out[:, 0] <= 0.05]
    print("periods of hits:", np.round(np.sort(hits[:, 1]), 0))
