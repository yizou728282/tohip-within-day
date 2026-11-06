# Deviation log: preregistered within-day study (OSF osf.io/5vfyq)

**Frozen materials:** Git tag `prereg-v2`, commit b824d746576fb91788514d669b4ea35ef19f5d5a; site list SHA-256 bb07b254…a8d04.
**Registered:** 2026-11-05T14:51:32Z.

The protocol requires every deviation to be logged with its date, reason and impact, and reported in the paper. This file covers that, and also lists:

- analyses that the protocol pre-specified but that were run differently or not run;
- choices the protocol left open, recorded before results were computed;
- analyses added after the confirmation results were seen.

Times are UTC.

## Summary

| ID | Date | Type | What | Effect on any preregistered verdict |
|---|---|---|---|---|
| D1 | 2026-11-05 | Operational | `fetch_rh.py` does not create its output folder; folder created by hand and script re-run | None |
| D2 | 2026-11-05 | Bug (download script) | HTTP 429 made the frozen `fetch_rh.py` write Cape_San_Juan's RH into 9 marine files; files re-downloaded | None (pre- and post-fix results both reported) |
| D3 | 2026-11-06 | Planned analysis not run as specified | Sensitivity #6 (ERA5 925/850 hPa RH via CDS) could not run: no CDS access. A substitute with GFS pressure-level RH was run instead | Not a verdict test; substitute reported as exploratory |
| D4 | 2026-11-06 | Reporting completed after the run | Sensitivity analyses evaluated for H1a–H3 from the frozen outputs (the frozen script pooled them but did not test them) | None; H1b fails only under the κ-Köhler predictor |
| D5 | 2026-11-06 | Open choices in pre-specified exploratory analyses | Meta-regression and f(RH) validation details not fixed by the protocol; fixed before computing | None (exploratory) |
| D6 | 2026-11-05 to 06 | Analyses added after unblinding | Joint current + lag-7 RH model; placebo decomposition by class; Mie + κ-Köhler forward check | None; all labelled post hoc |

Applied as registered, with nothing to log:

- No site was dropped for insufficient data. The protocol requires dropped site–variable combinations to be recorded in advance, and there were none.
- No outlier removal.
- No change to the site list, classes or model after the data were downloaded.

## D1 — Missing output folder (operational)

- **When:** 2026-11-05T15:24:38Z.
- **What:** the frozen `fetch_rh.py` writes into a folder it does not create, so the first run failed.
- **Action:** the folder was created by hand and the unchanged script was re-run.
- **Effect:** none on data or results. Logged in `confirmation/run_log.txt`.

## D2 — Rate-limit failure wrote another site's RH (bug in a frozen download script)

- **When:** found and fixed 2026-11-05T16:04:45Z to 16:10:42Z, before the analysis started at 16:10:42Z.
- **What happened:** for 9 marine sites, Open-Meteo answered all six retries with HTTP 429. The frozen `fetch_rh.py` then saved the previous site's series (Cape_San_Juan) under the 9 sites' names: El_Arenosillo, Gozo, Guadeloup, La_Parguera, Lampedusa, Learmonth, Maldives_Gan, Medenine-IRA, Monterey.
- **How it was found:** a content-hash check before the analysis.
- **Fix:** no code was changed. The 9 files were deleted and downloaded again with one site per process, so a failed download raises an error instead of reusing stale data. Before the analysis started, the full RH series was checked against the 2 m RH in the hourly files at all 72 sites.
- **Scope:** the full RH series is used only for the lag-7 placebo and for the rising/falling split. The primary model, H1–H3 and S1–S5 use the 2 m RH in the hourly files, which was never affected.
- **Before vs after:** the protocol requires reporting results before and after a bug fix. The frozen code was re-run on a reconstruction of the corrupted state (`prefix_rerun/`).
  - No hypothesis verdict changes, and S6 fails in both runs.
  - Share of sites with a significant placebo, pre-fix vs post-fix: AE 11.1 % vs 13.9 %; α_f 15.3 % vs 15.3 %; FMF 13.9 % vs 16.7 %.
  - Marine pooled AE placebo: −0.0005 pre-fix vs −0.0034 post-fix.

## D3 — Sensitivity #6 (pressure-level RH) not runnable as registered; substitute run

- **Registered:** "若可访问 CDS：改用 ERA5 925 hPa 与 850 hPa 的逐小时 RH" (if CDS is accessible, use hourly ERA5 RH at 925 and 850 hPa).
- **Why it could not run (2026-11-06):** this machine has no Copernicus CDS account or API key. The Open-Meteo ERA5 archive returns no pressure-level RH; its `relative_humidity_925hPa` is undefined. Creating a CDS account and accepting its licence is the author's decision.
- **Substitute (exploratory, not registered):** hourly RH from GFS at 2 m, 925 hPa and 850 hPa, from the Open-Meteo historical-forecast API (model `gfs_global`). This is a different model and a shorter period: 2021-03-23 to 2024-12-31, the earliest date these fields exist there.
- **Design, fixed before running (`sensitivity_plevel/plevel_analysis.py`):**
  - Sites: pollution and marine only, the classes H1–H3 need.
  - The frozen `fit`/`meta` functions are refitted with four RH predictors on one common set of hours.
  - A pressure level must lie at least 100 m above the site.
  - MIN_DAYS = 200 and MIN_H = 3 are unchanged; sites that fall below MIN_DAYS in the shorter period are listed in `plevel_coverage.csv`.
- **Download script:** `fetch_plevel.py` raises an error on a failed site instead of reusing data.
- **Download:** 2026-11-06T00:55:55Z to 02:41:15Z. All 47 files are complete (33,120 hours each) and have distinct content. Three empty responses were retried successfully.
- **Coverage:** 18 of 25 pollution sites and 19 of 22 marine sites entered the analysis. The rest dropped out:
  - no AERONET Level 2.0 data after 2021-03: Kyiv, Moscow_MSU_MO, Frioul, Guadeloup, Trinidad_Head;
  - 925 hPa less than 100 m above the site: Bandung, UNC-Palmira;
  - fewer than 200 qualifying days: Berlin_FUB, Cabauw, Punta_Arenas_UMAG.
- **Result:** pooled slopes per +10 % RH on identical hours (`sensitivity_plevel/results/`):

  | RH predictor | Pollution α_f (H1a) | Pollution FMF (H1b) | Marine AE (H2) | Median within-day SD of RH |
  |---|---|---|---|---|
  | 2 m, Open-Meteo archive (registered) | −0.0153 [−0.0261, −0.0045] | +0.0057 [+0.0022, +0.0093] | −0.0007 [−0.0070, +0.0057] | 11.9 % |
  | 2 m, GFS | −0.0221 [−0.0320, −0.0122] | +0.0050 [+0.0019, +0.0081] | −0.0065 [−0.0140, +0.0010] | 10.5 % |
  | 925 hPa, GFS | −0.0099 [−0.0166, −0.0032] | +0.0023 [−0.0001, +0.0048] | −0.0026 [−0.0082, +0.0030] | 7.6 % |
  | 850 hPa, GFS | −0.0042 [−0.0079, −0.0006] | +0.0037 [+0.0024, +0.0051] | +0.0035 [−0.0000, +0.0070] | 8.8 % |

- **Reading:**
  - The H1a sign holds with every predictor.
  - Pressure-level RH gives smaller slopes than 2 m RH, not larger. So the idea that the response is weak because 2 m RH misses the humidity of the aerosol layer is not supported by this substitute.
  - The two 2 m sources differ by about 45 % in the α_f slope (−0.015 vs −0.022). Error in the RH predictor itself therefore attenuates the slope substantially.
  - H2 is still not significant, and no marine result here is a test.
- **If the author sets up CDS access**, the registered ERA5 version can still be run, and it will be reported alongside this substitute.

## D4 — Sensitivity analyses tested for H1a–H3 (reporting completed)

- **Registered:** "敏感性分析（对 H1a–H3 逐一重算）" (recompute H1a–H3 under each sensitivity analysis).
- **Gap:** the frozen script writes pooled estimates for every specification (`meta_results.csv`) but evaluates the hypotheses only for the primary one.
- **What was done:** H1a–H3 were recomputed from the frozen outputs with the same rules (one-sided HKSJ tests; H3 SE as registered; Holm within each specification), without changing any frozen code (`exploratory/sensitivity_hypotheses.py`, `sensitivity_hypotheses.csv`).
- **Result:**
  - H1a is supported in all 6 specifications (Holm p ≤ 6 × 10⁻⁴).
  - H1b is supported in 5 of 6. It fails only with the κ-Köhler predictor RH/(100 − RH): Holm p = 0.33.
  - H2 and H3 are not supported in any specification. The lowest is H3 with the ln(AOD) control: raw p = 0.045, Holm p = 0.091.

## D5 — Choices left open in pre-specified exploratory analyses

Both analyses below are in the protocol, but parts of their specification were not. The choices were fixed in the script headers before any slope was computed for the new sites.

**Meta-regression of heterogeneity** (`exploratory/meta_regression.py`):
- Moderators standardised as z-scores.
- Coast distance entered as log10(1 + km).
- Mean RH taken over the analysed hours.
- τ² estimated by REML, with Knapp–Hartung standard errors.
- Two scopes: all 72 sites, and the pollution class alone.

Result:
- Across all sites, sites with higher climatological FMF had a more negative α_f slope (−0.0070 per SD, p = 0.020) and a more positive FMF slope (+0.0019 per SD, p = 0.033).
- Mean RH and coast distance explained nothing.
- Together the moderators explained 5–6 % of τ².
- Within the pollution class, α_f heterogeneity was unexplained.

**f(RH) external validation** (`exploratory/frh/`):
- The protocol says the site list "按 Burgos et al. 2019 的数据集确定" (is determined by the Burgos et al. 2019 dataset).
- f(RH) values: Burgos et al. 2019, Table 5, f(RH=85%/RHdry), PM10 or whole air, median. Taken from the article via Europe PMC (PMC6706437), not from the raw EBAS files.
- Site matching: the nearest AERONET site within 50 km with at least 300 days of Level 2.0 SDA data in 2015–2024. This gave 13 pairs.
- Test: one-sided Spearman correlation.
- Nephelometer periods do not overlap the AERONET years at most sites.

Result:
- α_f slope versus f(RH=85%): ρ = −0.22, in the predicted direction but not significant (one-sided p = 0.23, n = 13).
- Descriptively, the α_f slope was negative at 11 of the 13 sites. These include fine-dominated sites outside the confirmation set: ARM_SGP −0.021, Granada −0.021, Hyytiälä −0.021, Melpitz −0.032.

## D6 — Analyses added after unblinding (post hoc)

None of these were in the protocol. They are reported as post hoc and do not change any preregistered verdict.

- **Joint current + lag-7 RH model** (`confirmation/posthoc_not_preregistered/`): the pollution-site current-RH slopes are unchanged and the lag-7 terms are ≈ 0. The marine AE signal splits between the two terms.
- **S6 failures broken down by class:** pollution 7/75, marine 13/66, dust 11/39, burning 2/36 site–variable pairs.
- **Mie + κ-Köhler forward check** (`mie_check/`): for κ = 0.2–0.4, the observed α_f slope is 26–42 % of the model's. The model's FMF/α_f slope ratio matches the observed one when the coarse mode is slightly hygroscopic. The model does not reproduce the differences between sites.

## Effect on what the paper may claim

Under the registered rule, the failed placebo (S6) means H1a and H1b are reported as "可能受日变化伪影影响" (possibly affected by diurnal artefacts) and not counted as valid evidence. The paper must say so wherever the H1 result is stated.

Everything listed under D3, D5 and D6 is exploratory and cannot lift that flag. The paper must also report:
- the H1b failure under the κ-Köhler predictor (D4);
- the null f(RH) correlation (D5);
- that pressure-level RH weakens the response rather than strengthening it (D3).

D3 bears directly on the manuscript's "attenuated" framing. The response is still smaller than the Mie model predicts (D6). But the explanation that 2 m RH misses the humidity of the aerosol layer is not supported. The evidence points instead to error in the RH predictor (D3: the two 2 m sources differ by 45 %), or to lower hygroscopicity or efflorescence. The manuscript should say this.
