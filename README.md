# What de-identified timestamps support in MIMIC-IV and eICU-CRD

Analysis code and aggregate outputs for:

> Mikkelsen Y. What De-identified Timestamps Support in MIMIC-IV and eICU-CRD:
> An Empirical Audit of Calendar Position, Concurrency, Clock Time and Precision.
> Manuscript, 2026.

De-identification changes what the timestamps of MIMIC-IV and eICU-CRD can be
used for. This repository tests, for each of four temporal properties, whether
the timestamps preserve it: **calendar position** (weekday, season),
**cross-patient concurrency** (census), **time of day** (including the clock
reference of eICU-CRD), and **timestamp precision** (date-only placeholders and
heaping). Each check compares an observed distribution with the distribution
expected if the property were preserved and the distribution expected if
de-identification had destroyed it. A worked example repeats a published
eICU-CRD clock-time analysis on the recorded clock and in local time.

---

## Checks

| Check | Database | What is tested |
|---|---|---|
| M1 | MIMIC-IV | Weekday of each patient's first scheduled admission (ELECTIVE, SURGICAL SAME DAY ADMISSION) against 1/7 per day |
| M2 | MIMIC-IV | Month of each patient's first influenza admission (ICD-9 487–488, ICD-10 J09–J11) against month length |
| M3 | MIMIC-IV | Hour of laboratory specimen times (`labevents.charttime`), per result row and per specimen (distinct `specimen_id`), and of scheduled admission times |
| M4 | MIMIC-IV | ICU census rebuilt from shifted `icustays` times against the census implied by the real span (`anchor_year_group`) |
| M5 | MIMIC-IV | Times at exactly 00:00:00 and exactly on the hour, and times with non-zero seconds: `admittime` by admission type, `icustays` intime/outtime |
| E1 | eICU-CRD | Every calendar or clock column in all 31 tables |
| E2 | eICU-CRD | Hour of unit admission and of laboratory results (`unitadmittime24` + `labresultoffset`) |
| E3 | eICU-CRD | A census built from offsets (shown to be the length-of-stay survival curve) |
| E4 | eICU-CRD | Placeholder, heaping and seconds shares of the four `*time24` fields |
| E5 | eICU-CRD | Clock reference, laboratory anchor: lag of each hospital's 10-minute laboratory-result profile against the pooled Northeast profile, by region |
| E6 | both | Clock reference, physiological anchor: lag of eICU-CRD heart-rate and respiratory-rate day–night profiles (recorded clock) against MIMIC-IV (Boston local time), per hospital and region. Two estimates: whole-hour circular cross-correlation, and first-harmonic cosinor acrophase difference corrected for within-hour timing |
| E7 | both | Clock reference, administrative anchor: lag of eICU-CRD live hospital-discharge and ICU-discharge time profiles (10-minute bins, exact midnight excluded) against MIMIC-IV, by cross-correlation and cosinor, per hospital and region |
| E8 | eICU-CRD | Worked example: the design of a published analysis of heart-rate peak time and in-hospital mortality (Yang et al., *Heart Rhythm* 2022, doi:10.1016/j.hrthm.2022.03.1230), repeated on the recorded clock and in local time. Per-stay cosinor peak time in each patient's first unit stay; rhythmic stays (zero-amplitude F test P < .05); peak 12:00–17:59 against all other times; logistic regression unadjusted and adjusted for the APACHE IV score, hospital-clustered 95% CIs. Local time: regional midpoint of the expected UTC offset, ±1 hour, and the E7 hospital-discharge offsets. Stroke subgroup: APACHE admission diagnosis "CVA, cerebrovascular accident/stroke" |

## Main results (run of 20 Oct 2026, script v11; `results/paper26_summary.txt`)

- **Weekday (MIMIC-IV).** 29.07% (95% CI 28.64–29.51%) of 41,683 first scheduled admissions fall on a weekend, slightly above the 28.57% expected if weekday is scrambled.
- **Season (MIMIC-IV).** 35.1% of 2,187 first influenza admissions fall in December–March, against 33.2% if season is scrambled (χ² P = .47).
- **Concurrency.** The census rebuilt from shifted MIMIC-IV times averages 9.0 ICU stays and is 0 for 7.8% of the time. The census implied by the real 2008–2022 span averages 62.6. The 7.0-fold ratio of the means is the ratio of the shifted to the real span. An eICU-CRD census built from offsets is the length-of-stay survival curve.
- **Time of day.** MIMIC-IV laboratory specimen times peak at 06:00. eICU-CRD laboratory result times peak at 10:00 on the recorded clock.
- **eICU-CRD clock reference.** Median hospital cosinor lags against MIMIC-IV local time, by region (Northeast, Midwest, South, West):

  | Anchor | Northeast | Midwest | South | West |
  |---|---|---|---|---|
  | Heart rate | +2.8 h | +3.9 h | +3.6 h | +6.3 h |
  | Respiratory rate | +4.0 h | +4.6 h | +4.5 h | +6.6 h |
  | Hospital discharge | +3.3 h | +3.9 h | +4.1 h | +6.4 h |
  | Expected if UTC | +4 to +5 h | +4 to +6 h | +4 to +6 h | +6 to +8 h |

  These exclude local time (expected 0) and a single US Eastern clock (expected 0, 0 to +1, 0 to +1, +2 to +3 h). They lie within the UTC ranges or up to 1.2 h below them. ICU discharge times matched MIMIC-IV less closely (pooled r 0.86–0.90) and fitted no single clock reference (+0.5, +2.6, +2.4, +4.3 h).
- **Placeholders and heaping (MIMIC-IV).**
  - 32.5% of SURGICAL SAME DAY ADMISSION and 7.1% of ELECTIVE admission times are exactly 00:00:00, against 0.07% expected.
  - Laboratory times are exactly 00:00:00 in 2.2% of result rows (3,482,197) and 1.25% of specimens (219,018 of 17,503,481). They fall exactly on the hour in 9.4% and 9.3% of other hours (1.7% expected).
  - ICU intime and outtime are mostly recorded to the second and carry no midnight placeholders.
- **Worked example (eICU-CRD).** Of 73,291 rhythmic first stays in 169 hospitals, 16.9% peaked between 12:00 and 17:59 on the recorded clock and had higher in-hospital mortality than other stays (11.1% vs 8.2%; APACHE IV-adjusted OR 1.28, 95% CI 1.19–1.38), in the published direction. In local time (regional midpoint offset), 86.1% of stays fell in a different 6-hour window and the association reversed (adjusted OR 0.79, 95% CI 0.74–0.85). It reversed under every offset tested (adjusted OR 0.76–0.86) and among 3,229 stroke admissions (adjusted OR 1.69 on the recorded clock, 0.59 in local time).

## Data

Neither database is redistributed. Both require credentialed access through
PhysioNet under their data use agreements:

- **MIMIC-IV v3.1:** https://physionet.org/content/mimiciv/3.1/ (doi:10.13026/kpb9-mt58)
- **eICU-CRD v2.0:** https://physionet.org/content/eicu-crd/2.0/ (doi:10.13026/C2WM1R)

Everything in `results/` and `figures/` is aggregate. No record-level data are
included. eICU-CRD `hospitalid` is not written to any output by default.
- **Hospital-level files** (`e5_hospital_lags.csv`, `e6_hospital_lags.csv`, `e7_hospital_lags.csv`) carry only region, volume, lag and fit statistics.
- **Row order** is sorted, so the original hospital order cannot be recovered.
- **E8 files** (`e8_peak_time_windows.csv`, `e8_models.csv`) carry only counts, proportions and odds ratios by cohort, clock and window.
- **Local runs:** `--keep-hospital-ids` restores the identifier for local use only.

## Running

Requirements: Python ≥ 3.11 and [uv](https://docs.astral.sh/uv/). Dependencies
(pandas, numpy, scipy, statsmodels, matplotlib) are declared inline with
[PEP 723](https://peps.python.org/pep-0723/). No GPU is needed.

```bash
uv run scripts/paper26_timestamp_checks_v11.py \
    --mimic-root ~/physionet.org/files/mimiciv/3.1 \
    --eicu-root  ~/physionet.org/files/eicu-crd/2.0 \
    --out-dir    ~/paper26_results
```

| Argument | Default | Meaning |
|---|---|---|
| `--mimic-root` | `~/physionet.org/files/mimiciv/3.1` | MIMIC-IV root (with `hosp/` and `icu/`) |
| `--eicu-root` | `~/physionet.org/files/eicu-crd/2.0` | eICU-CRD root |
| `--out-dir` | `~/paper26_results` | Output folder |
| `--chunksize` | 5,000,000 | Rows per chunk for the large tables |
| `--max-lab-rows` | 0 (all) | Stop each laboratory pass after N rows (testing) |
| `--min-hospital-labs` | 2,000 | E5: minimum laboratory results per hospital |
| `--reference-region` | Northeast | E5: reference region |
| `--skip-mimic`, `--skip-eicu` | off | Run one database only |
| `--skip-e6`, `--only-e6` | off | Skip, or run only, the physiological anchor (heavy passes) |
| `--skip-e7`, `--only-e7` | off | Skip, or run only, the discharge-time anchor (light) |
| `--skip-e8` | off | Skip the worked example (E8 uses the E6 heart-rate profiles, so it also does not run when E6 is skipped) |
| `--keep-hospital-ids` | off | Write eICU-CRD `hospitalid` (local use only) |

The full run took about 9 minutes (`results/run_manifest.txt` records the
versions, arguments and input files). Most of that time goes on single passes
over `labevents`, `chartevents`, eICU-CRD `lab` and `vitalPeriodic`. E8 also
reads eICU-CRD `apachePatientResult`.

## Testing on synthetic data

`synthetic/make_mock_v4.py` writes MIMIC-IV- and eICU-CRD-shaped tables with
planted effects (its docstring names script v10; v11 reads the same tables):

- weekday-clustered scheduled admissions and winter influenza;
- morning laboratory peaks, with results grouped into specimens (`specimen_id`);
- midnight placeholders (30% of same-day surgical admissions, 3% of laboratory times);
- heart-rate and respiratory-rate day–night rhythms;
- discharge times with a midday peak in local time, with discharge status and `patienthealthsystemstayid`;
- an alarm-limit item (220046) that must be excluded;
- for E8: a per-stay heart-rate peak in local time (25% of stays at 09:00, the rest at 16:00), higher in-hospital mortality with the morning peak, `uniquepid`, `unitvisitnumber`, `apacheadmissiondx` and an `apachePatientResult` table;
- the missing values found in the real eICU-CRD tables: hospitals with no region and stays with no `hospitaldischargestatus`;
- optionally, region-specific UTC offsets for the eICU-CRD clock.

```bash
uv run synthetic/make_mock_v4.py scrambled mock_utc utc     # or: preserved <dir> local
uv run scripts/paper26_timestamp_checks_v11.py \
    --mimic-root mock_utc/mimiciv/3.1 --eicu-root mock_utc/eicu-crd/2.0 \
    --out-dir mock_out
```

`scrambled` shifts each patient by a random number of days, so weekday and
season are destroyed. `preserved` shifts by whole weeks close to whole years.
The checks should recover what was planted, and return the destroyed-property
expectation when an effect is absent. With `utc`, the cosinor and discharge
lags recover the planted regional offsets; with `local`, they are about 0.
With `utc`, E8 shows the planted morning-peak mortality excess under 12:00–17:59
on the recorded clock (OR about 4) and under 06:00–11:59 after conversion, where
12:00–17:59 has the lower mortality (OR about 0.24).

## Outputs

| File | Content |
|---|---|
| `paper26_summary.txt` | All results in text form |
| `run_manifest.txt` | Script, versions, arguments, input files with sizes and dates |
| `m1_weekday.csv` | Weekday counts (first scheduled admission per patient; all scheduled; first admission of any type) |
| `m2_month.csv` | Monthly influenza counts and month-length expectation |
| `m3_hour.csv` | Hourly laboratory specimen times (per result row), on-the-hour counts, scheduled admission times |
| `m3_specimen_hour.csv` | Hourly laboratory specimen times and on-the-hour counts, one per specimen |
| `m4_census_shifted.csv`, `m4_census_summary.csv` | Distribution and summary of the census rebuilt from shifted times |
| `m5_placeholders.csv`, `e4_placeholders.csv` | Exactly-midnight, on-the-hour and non-zero-seconds shares by field |
| `e1_calendar_fields.csv` | Calendar and clock columns in eICU-CRD |
| `e2_hour.csv` | Hourly eICU-CRD unit admissions and laboratory results (recorded clock) |
| `e3_offset_census.csv` | Stays present at each hour after unit admission |
| `e5_hospital_lags.csv`, `e5_region_summary.csv` | Laboratory-anchor lags per hospital (no identifier) and by region |
| `e6_hospital_lags.csv`, `e6_region_summary.csv`, `e6_profiles.csv` | Physiological-anchor lags (cross-correlation and cosinor) and hourly profiles |
| `e7_hospital_lags.csv`, `e7_region_summary.csv`, `e7_profiles.csv` | Discharge-time anchor lags (cross-correlation and cosinor) and hourly profiles |
| `e8_peak_time_windows.csv` | Worked example: stays, deaths and mortality by 6-hour peak-time window, per cohort and clock |
| `e8_models.csv` | Worked example: share and mortality with a 12:00–17:59 peak, share of stays moved to another window, unadjusted and APACHE IV-adjusted odds ratios with 95% CIs, per cohort and clock |
| `figures/paper26_fig1.png` – `paper26_fig4.png` | Manuscript Figures 1–4 (300 dpi) |

## Version

This release (1.1.0) corresponds to script version 11. It is the version used
for the manuscript. All outputs of checks M1–M5 and E1–E7, and all four
figures, are byte-identical to release 1.0.0 (script v8); E8 is new.
- **v11 change:** E8 stroke subgroup defined by the exact APACHE admission diagnosis "CVA, cerebrovascular accident/stroke". The v10 text search for "stroke" also matched 11 stays with "Heat exhaustion/stroke".
- **v10 fix:** E8 failed on missing (NA) `hospitaldischargestatus` or region values ("boolean value of NA is ambiguous"). Missing values are now treated as False, and the summary reports how many stays lack an outcome or region.
- **v9 addition:** E8 worked example, with outputs `e8_peak_time_windows.csv` and `e8_models.csv`, and `--skip-e8`. Adds scipy as a dependency.
- **v8 change:** Figure 1C y-axis label corrected to "Share of result rows" (the panel counts `labevents` rows).
- **v7 additions:** per-specimen laboratory counts (M3); cosinor phase lags with a correction for within-hour timing (E6); discharge-time anchor (E7) and Figure 4.
- **v6 additions:** non-zero-seconds shares; expected UTC ranges corrected for the time zones each region spans; figures at 300 dpi; no hospital identifiers in outputs.
- **v1–v5:** development versions, not included. Their computations are unchanged in v11.

## Use of generative AI

Claude (Anthropic) was used to help write and test the code. No credentialed
data were shared with the tool. The analysis was run locally by the author,
who reviewed and verified all code and outputs.

## Licence

MIT; see `LICENSE`.
