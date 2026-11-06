# Pre-fix re-run of the confirmatory analysis (deviation D2)

The protocol requires that when a bug is fixed, results before and after the fix are both reported.

The bug was in the frozen `fetch_rh.py`. For 9 marine sites the Open-Meteo archive returned HTTP 429 on all retries, and the script then wrote the previous site's (Cape_San_Juan) full RH series into those sites' files:

> El_Arenosillo, Gozo, Guadeloup, La_Parguera, Lampedusa, Learmonth, Maldives_Gan, Medenine-IRA, Monterey

The reported results use correct files obtained by re-downloading those sites one at a time. This folder holds the same frozen code (`prereg-v2`, `git archive`) run on a reconstruction of the corrupted state, made by copying Cape_San_Juan.csv over the 9 files. Hourly AERONET files and the 2 m RH in them were never affected.

## What the corruption affected

The full RH series feeds only two specifications:
- the lag-7-day placebo (`RH_lag7`);
- the rising/falling split (`RH_prev1h`).

So only marine-class rows of `placebo_lag7`, `expl_fall` and `expl_rise_minus_fall` differ between the pre-fix and post-fix runs.

## Result

| Quantity | Pre-fix | Reported (post-fix) |
|---|---|---|
| H1a, H1b | supported, Holm p = 2 × 10⁻⁵ | same |
| H2, H3 | not supported, Holm p = 0.17 | same |
| S6, share of sites with a significant placebo: AE | 11.1 % (fail) | 13.9 % (fail) |
| S6, share: α_f | 15.3 % (fail) | 15.3 % (fail) |
| S6, share: FMF | 13.9 % (fail) | 16.7 % (fail) |
| Marine pooled placebo, AE | −0.0005 | −0.0034 |

No verdict changes. The marine AE placebo signal reported in the paper appears only with the correct data. With the corrupted series, nine marine sites were regressed on another island's humidity, which is close to a pure null.
