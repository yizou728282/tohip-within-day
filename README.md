# Within-day humidity response of AERONET optics

Code, the frozen 72-site list and the derived results for the preregistered study at https://osf.io/5vfyq. The question is whether within-day changes in relative humidity produce mode-specific responses in AERONET intensive optics: fine-mode growth at pollution sites (alpha_f down, FMF up, total AE up) and coarse sea-salt growth at marine sites (total AE down, FMF down), with dust and biomass-burning sites near zero.

These files are the `within_day` directory of commit 7a4f326dfd6960e0304f6894d129a905b167a0ec (2026-11-06). The manuscript draft is not included. Raw AERONET and Open-Meteo downloads are not included.

## Layout

| Path | Role |
|------|------|
| `discovery/` | 18 discovery sites and every result seen before preregistration |
| `prereg/confirmation_sites.csv` | Frozen confirmation site list (72 sites; SHA-256 of this LF file: bb07b254d3920a716e9171f219e09a835bbc0d9782c8c52b2fbcd459a97a8d04) |
| `scripts/select_sites.py` | Site selection from metadata and daily Level 2.0 climatology only |
| `scripts/fetch_daily_all.py` | All-network daily Level 2.0 AOD/SDA, input to site selection |
| `scripts/fetch_hourly.py` | Hourly-averaged Level 2.0 all-points AOD/AE/PW plus SDA FMF/alpha_f, matched to Open-Meteo hourly RH |
| `scripts/fetch_rh.py` | Full hourly RH series, used for the lag-7-day placebo and the rising/falling split |
| `scripts/confirm_analysis.py` | Preregistered confirmatory analysis: site fits, meta-analysis, Holm, TOST, sensitivity |
| `scripts/within_day2.py`, `scripts/robust.py` | Discovery-stage scripts, not used for confirmation |
| `confirmation/` | Confirmatory run on the 72 sites: results, run log, summary, post-hoc joint model |
| `deviations/DEVIATIONS.md` | Deviation log (D1 to D6), including the pre-fix re-run in `deviations/prefix_rerun/` |
| `exploratory/` | Meta-regression, f(RH) comparison, and H1a to H3 under each sensitivity specification |
| `sensitivity_plevel/` | Substitute for sensitivity 6: GFS pressure-level RH (deviation D3) |
| `mie_check/` | Post-hoc Mie and kappa-Koehler check |
| `manuscript/figures/` | Site map and within-day slope figure |
| `figures/site_map.py` | Script that draws the site map |

## Tables

In `confirmation/site_results.csv`, `beta` is the within-day slope per +10 % RH, except `sens_xk`, where the predictor is RH/(100-RH) and the slope is not rescaled. `lo` and `hi` are the 2.5 and 97.5 percentiles of the day-clustered bootstrap. `confirmation/meta_results.csv` pools those slopes. The run that wrote the confirmation files is in `confirmation/run_log.txt`. What changed after the freeze is in `deviations/DEVIATIONS.md`.

## Reproduce

```bash
python scripts/fetch_hourly.py prereg/confirmation_sites.csv data/hourly
python scripts/fetch_rh.py prereg/confirmation_sites.csv data/rh_full
python scripts/confirm_analysis.py prereg/confirmation_sites.csv data/hourly data/rh_full results/confirmation
```

Python 3.10 or newer. Packages are listed in `requirements.txt`. AERONET Version 3 Level 2.0 is public at https://aeronet.gsfc.nasa.gov. Hourly RH is taken from the Open-Meteo archive API (https://open-meteo.com).

`select_sites.py` also needs local AERONET location and daily files, and Natural Earth coastline and populated-place files. Those inputs are not in this repository. The confirmation list is the CSV above. `figures/site_map.py` needs a Natural Earth 1:10 m coastline GeoJSON, which is not included; the map files already in `manuscript/figures/` were drawn with it.

## Preregistration

https://osf.io/5vfyq was registered at 2026-11-05T14:51:32Z. The frozen code is tag `prereg-v2`, commit b824d746576fb91788514d669b4ea35ef19f5d5a. That tag adds a descriptive dust-core sensitivity (median AE < 0.6) to S3. `prereg-v1` is an earlier freeze of the same site list. The header of `scripts/confirm_analysis.py` still says `prereg-v1`.

Hourly files for the confirmation sites were downloaded after that timestamp. The times are in `confirmation/run_log.txt`. Later changes are logged deviations.

## Licence

The code and the text files here are under the MIT licence. AERONET and Open-Meteo data are not redistributed in this repository and stay under the terms of those sources.
