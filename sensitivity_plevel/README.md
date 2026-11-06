# Substitute for preregistered sensitivity #6: pressure-level RH (deviation D3)

The registered version used ERA5 925/850 hPa RH via CDS, which is not accessible here. This folder holds a substitute that uses GFS pressure-level RH from the Open-Meteo historical-forecast API, covering 2021-03-23 to 2024-12-31. It is exploratory. The design, coverage and results are in [`../deviations/DEVIATIONS.md`](../deviations/DEVIATIONS.md), entry D3.

## Files

- `fetch_plevel.py` — download script. A failed site raises an error instead of reusing another site's data.
- `fetch_log.txt`, `fetch_start_utc.txt`, `fetch_end_utc.txt` — download record.
- `plevel_analysis.py` — refits the frozen model with four RH predictors on identical hours.
- `results/` — outputs:
  - `plevel_site_results.csv`: site-level slopes;
  - `plevel_meta.csv`: pooled slopes by class;
  - `plevel_hypotheses.csv`: H1a–H3 for each predictor;
  - `plevel_coverage.csv`: which sites were kept or dropped.
- `gfs_rh/` — the downloaded GFS series for 47 sites (110 MB). It is not committed, consistent with the rest of the repository, which keeps no raw downloads; `fetch_plevel.py` re-creates it.
