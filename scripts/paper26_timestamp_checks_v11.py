#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["pandas", "numpy", "scipy", "statsmodels", "matplotlib"]
# ///
"""
Paper 26 - what MIMIC-IV and eICU-CRD timestamps can and cannot support.

Each check compares an observed distribution with what de-identification would
produce if the calendar or clock property had been scrambled.

MIMIC-IV
  M1  Day of week   Scheduled admissions (ELECTIVE, SURGICAL SAME DAY
                    ADMISSION) by weekday of admittime. Real scheduled
                    admissions cluster on weekdays; a scrambled weekday puts
                    about 2/7 of them on Saturday or Sunday.
  M2  Season        Admissions carrying an influenza diagnosis (ICD-9 487-488,
                    ICD-10 J09-J11) by month of admittime. Real influenza peaks
                    in winter; a scrambled season spreads admissions across
                    months in proportion to days per month.
  M3  Time of day   Laboratory specimen times (hosp/labevents.charttime) and
                    scheduled-admission times by hour. A preserved clock shows
                    routine peaks; a scrambled clock is flat.
  M4  Concurrency   ICU census reconstructed from shifted icustays times,
                    against the mean census implied by total ICU stay-time over
                    the real calendar span given by anchor_year_group.
  M5  Placeholder   Share of timestamps at exactly 00:00:00 and exactly on the
      times         hour: admissions by admission_type, icustays intime and
                    outtime. (The labevents equivalent is computed in the M3
                    pass.) With genuine minute-resolution times, 1/60 of the
                    times in any hour fall at minute 00.

eICU-CRD
  E1  Calendar fields  Every column, in every table, whose name suggests a
                       calendar or clock field.
  E2  Time of day      Unit admission clock time (unitadmittime24) and
                       laboratory result clock time (unitadmittime24 +
                       labresultoffset) by hour.
  E3  Concurrency      A "census" built from offsets: stays present at each
                       hour after unit admission. Every stay starts at offset
                       0, so the curve is the length-of-stay survival curve,
                       not a census.
  E4  Placeholder      Share of exactly 00:00:00 and on-the-hour values in the
      times            patient-table clock fields (*time24).

  E5  Clock reference  Are the eICU *time24 clock times local time or UTC?
                       Each hospital's laboratory-result clock-time profile
                       (10-minute bins) is aligned by circular
                       cross-correlation with the pooled profile of a
                       reference region (Northeast, leaving the hospital
                       itself out). If times are UTC, hospitals in later US
                       time zones sit later on the recorded clock (Midwest
                       and South 0 to +1 h, West +2 to +3 h; approximate,
                       as each region spans more than one zone); if times
                       are local, every region sits near 0. Secondary: eICU Northeast
                       result times against MIMIC-IV specimen times (Boston,
                       local), excluding exact-midnight placeholders: about
                       +5 to +6 h if eICU is UTC, about +1 h (laboratory
                       turnaround) if local.

v2 (20 Oct 2026): adds M5, E4, the midnight and on-the-hour counts in the
M3 labevents pass, and per-weekday adjusted residuals in M1.
  E6  Physiological    Absolute anchor for the eICU clock. Heart-rate (and
      anchor           respiratory-rate) day-night profiles, each stay's own
                       mean removed, by clock hour: MIMIC-IV chartevents
                       (itemids verified against d_items; Boston local time)
                       against eICU vitalPeriodic on its recorded clock
                       (unitadmittime24 + observationoffset), per hospital and
                       per region. If eICU *time24 is UTC: Northeast +4 to
                       +5 h, Midwest and South +4 to +6 h, West +6 to +8 h
                       (contiguous states; daylight saving time cannot be
                       resolved); a common US Eastern clock: Northeast 0,
                       West +2 to +3; local time: 0 everywhere.

v2 (20 Oct 2026): adds M5, E4, the midnight and on-the-hour counts in the
M3 labevents pass, and per-weekday adjusted residuals in M1.
v3 (20 Oct 2026): adds E5 and paper26_fig2.png; the M3 and E2 laboratory
passes now also accumulate 10-minute profiles (MIMIC-IV pooled, eICU per
hospital). All v2 outputs are unchanged.
v4 (20 Oct 2026): adds E6 and paper26_fig3.png. E6 adds one pass over
icu/chartevents and one over vitalPeriodic (the slowest steps of the run);
--only-e6 runs E6 alone, --skip-e6 omits it. All v3 outputs are unchanged.
v5 (20 Oct 2026): figure fixes only (figure 2B fixed axis with outlying
hospitals drawn as edge triangles; median labels clear of points in 2B and
3B) and an outlier count in the E5 summary. No computation changed.
  E7  Administrative  Discharge clock times as an anchor independent of
      anchor           physiology: live hospital discharges and live ICU
                       discharges, MIMIC-IV (dischtime, icustays.outtime;
                       Boston local time) against eICU-CRD
                       (hospitaldischargetime24, unitdischargetime24;
                       recorded clock), per hospital and per region.
                       Exact-midnight values are excluded on both sides.
  E8  Worked example  Consequence of the clock reference for a published
                       analysis design (Yang et al., Heart Rhythm 2022):
                       per-stay heart-rate cosinor peak time in 6-hour
                       windows and in-hospital mortality, on the recorded
                       clock and after conversion to local time with a
                       regional UTC offset (midpoint, -1 h, +1 h, and the E7
                       hospital-discharge anchor). One stay per patient,
                       rhythmic stays (zero-amplitude F test p < 0.05),
                       logistic regression unadjusted and adjusted for the
                       APACHE IV score, hospital-clustered 95% CIs.

v6 (20 Oct 2026): (1) M5 and E4 report the share of times with non-zero
seconds per field, because the 1/60 and 1/1440 expectations hold only for
minute-resolution times (second-resolution times: 1/3600 and 1/86400);
(2) E5 and E6 expected UTC ranges corrected for the time zones each region
spans (summary text and figure 3B key; no computation changed); (3) figures
saved at 300 dpi; (4) eICU-CRD hospitalid is no longer written to any output
(e5/e6_hospital_lags.csv, summary), so all outputs can be posted publicly;
--keep-hospital-ids restores it for local use. All other computations are
unchanged from v5.
v7 (20 Oct 2026): (1) M3 also counts laboratory specimens (distinct
specimen_id, time of the first row) as well as result rows; (2) E6 adds
sub-hour lags from first-harmonic cosinor acrophases (alongside the 1-hour
cross-correlation lags), corrected for the difference between the databases
in mean minute within the clock hour; (3) new check E7, an administrative anchor from
discharge clock times, with its own cosinor lags, and figure 4. All v6
computations and outputs are unchanged; v7 adds columns and files.
v8 (20 Oct 2026): figure 1C y-axis label corrected from 'Share of specimens'
to 'Share of result rows' (the panel counts labevents rows, not distinct
specimens). No computation or other output changed.
v9 (20 Oct 2026): adds E8 (worked example, outputs e8_peak_time_windows.csv
and e8_models.csv), computed from the E6 heart-rate stay profiles plus the
eICU-CRD apachePatientResult table; --skip-e8 omits it. Requires scipy
(already a statsmodels dependency). All v8 computations and outputs are
unchanged.
v10 (20 Oct 2026): E8 fix - missing hospitaldischargestatus or region values
(NA) raised "boolean value of NA is ambiguous"; masks now treat NA as False,
and the E8 summary reports the number of stays with a missing outcome or
region. No other change.
v11 (20 Oct 2026): E8 stroke subgroup defined by the exact APACHE admission
diagnosis 'CVA, cerebrovascular accident/stroke'; the v10 text search for
'stroke' also matched 11 stays with 'Heat exhaustion/stroke'. No other change.

Outputs in --out-dir: paper26_summary.txt, one CSV per check,
paper26_fig1.png to paper26_fig4.png (300 dpi) and run_manifest.txt.

Usage
  uv run paper26_timestamp_checks_v11.py \
      --mimic-root ~/physionet.org/files/mimiciv/3.1 \
      --eicu-root  ~/physionet.org/files/eicu-crd/2.0 \
      --out-dir    ~/paper26_results

Runtime is dominated by one pass over hosp/labevents and one over eICU lab.
--max-lab-rows N stops each of those passes after N rows (0 = read all).
A check that fails is reported in the summary; the others still run, and the
script exits with status 1.
"""

from __future__ import annotations

import argparse
import platform
import re
import sys
import time
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels  # noqa: E402
from statsmodels.stats.gof import chisquare  # noqa: E402
from statsmodels.stats.proportion import proportion_confint  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SCHEDULED_TYPES = ("ELECTIVE", "SURGICAL SAME DAY ADMISSION")
FLU_ICD9 = ("487", "488")
FLU_ICD10 = ("J09", "J10", "J11")
WEEKDAYS = ("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
DAYS_PER_MONTH = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
WINTER = (12, 1, 2, 3)
CALENDAR_PATTERN = re.compile(r"(year|date|month|week|time24|dayofweek)",
                              re.IGNORECASE)
CENSUS_HOURS_MAX = 720  # eICU offset "census" grid, hours after unit admission
BIN_MIN = 10            # clock-profile bin width, minutes
BINS = 1440 // BIN_MIN  # 144 bins per day
BOOT_SEED = 20261020
# E6 physiological anchor
VITALS = {  # name: (MIMIC itemid, expected d_items label, eICU column, range)
    "heart_rate": (220045, "Heart Rate", "heartrate", (20, 250)),
    "respiratory_rate": (220210, "Respiratory Rate", "respiration", (2, 80)),
}
MIN_STAY_HOURS = 18      # clock hours (of 24) a stay must cover
MIN_HOSPITAL_STAYS = 30  # stays per hospital for a per-hospital lag
AXIS_CLIP_HOURS = (3.5, 6.5)  # figure 2B axis: -3.5 h to +6.5 h
OUTLIER_HOURS = 6.0           # E5: |lag| beyond this is counted as outlying
KEEP_HOSPITAL_IDS = False     # set by --keep-hospital-ids; default omits IDs
FIG_DPI = 300
MIN_HOSPITAL_DISCHARGES = 50  # E7: discharges per hospital for a lag
REGION_COLOURS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
# E8 worked example
PEAK_WINDOWS = ((0, 6), (6, 12), (12, 18), (18, 24))  # clock-hour windows
RISK_WINDOW = (12, 18)  # Yang et al. 2022: highest mortality, peak 12-18 h
UTC_OFFSET_MID = {"Northeast": 4.5, "Midwest": 5.0, "South": 5.0,
                  "West": 7.0}  # midpoints of the expected UTC lags (h)
STROKE_DX = "CVA, cerebrovascular accident/stroke"  # exact APACHE label

INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
SERIES = "#2a78d6"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def find_table(root: Path, *candidates: str) -> Path:
    for c in candidates:
        p = root / c
        if p.exists():
            return p
        alt = p.with_suffix("")
        if alt.exists():
            return alt
    raise FileNotFoundError(f"none of {candidates} under {root}")


def read_table(path: Path, wanted: list[str], **kwargs):
    """Read only the wanted columns; fail loudly if any is absent."""
    header = pd.read_csv(path, nrows=0).columns.tolist()
    missing = [c for c in wanted if c not in header]
    if missing:
        raise KeyError(f"{path.name}: missing columns {missing}; "
                       f"header has {header}")
    return pd.read_csv(path, usecols=wanted, low_memory=False, **kwargs)


def share_ci(k: int, n: int) -> tuple[float, float, float]:
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return k / n, float(lo), float(hi)


def gof(observed: np.ndarray, expected_share: np.ndarray) -> tuple[float, float]:
    observed = np.asarray(observed, dtype=float)
    expected = expected_share / expected_share.sum() * observed.sum()
    stat, p = chisquare(observed, f_exp=expected)
    return float(stat), float(p)


def hour_profile(counts: np.ndarray) -> dict:
    """Peak hour, busiest circular 3-hour window and peak/trough ratio."""
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    share = counts / total
    window = np.array([share[[h, (h + 1) % 24, (h + 2) % 24]].sum()
                       for h in range(24)])
    start = int(window.argmax())
    trough = counts.min()
    return {
        "n": int(total),
        "peak_hour": int(counts.argmax()),
        "peak_window": f"{start:02d}:00-{(start + 2) % 24:02d}:59",
        "peak_window_share": float(window.max()),
        "peak_trough_ratio": float(counts.max() / trough) if trough > 0
        else float("inf"),
    }


def first_per_subject(frame: pd.DataFrame) -> pd.DataFrame:
    """Each patient's earliest admission in the frame (one date shift per
    patient, so later admissions of the same patient are not independent)."""
    return frame.sort_values(["subject_id", "admittime"]) \
        .drop_duplicates("subject_id", keep="first")


def hours_from_clock_strings(s: pd.Series) -> pd.Series:
    """'YYYY-MM-DD HH:MM:SS' -> integer hour; NaN if malformed.

    Fixed-position slice (fast on ~10^8 rows); rows that do not have the
    expected separators at positions 10 and 13 are returned as NaN."""
    s = s.astype("string")
    ok = (s.str.slice(10, 11) == " ") & (s.str.slice(13, 14) == ":")
    h = pd.to_numeric(s.str.slice(11, 13), errors="coerce")
    return h.where(ok.fillna(False) & h.between(0, 23))


def clock_parts_datetime(s: pd.Series):
    """'YYYY-MM-DD HH:MM:SS' -> (hour, minute, second); NaN where malformed."""
    s = s.astype("string")
    ok = ((s.str.slice(10, 11) == " ") & (s.str.slice(13, 14) == ":")
          & (s.str.slice(16, 17) == ":")).fillna(False)
    h = pd.to_numeric(s.str.slice(11, 13), errors="coerce")
    m = pd.to_numeric(s.str.slice(14, 16), errors="coerce")
    sec = pd.to_numeric(s.str.slice(17, 19), errors="coerce")
    ok = ok & h.between(0, 23) & m.between(0, 59) & sec.between(0, 59)
    return h.where(ok), m.where(ok), sec.where(ok)


def clock_parts_time24(s: pd.Series):
    """'HH:MM:SS' -> (hour, minute, second); NaN where malformed."""
    s = s.astype("string").str.strip()
    ok = ((s.str.slice(2, 3) == ":") & (s.str.slice(5, 6) == ":")).fillna(False)
    h = pd.to_numeric(s.str.slice(0, 2), errors="coerce")
    m = pd.to_numeric(s.str.slice(3, 5), errors="coerce")
    sec = pd.to_numeric(s.str.slice(6, 8), errors="coerce")
    ok = ok & h.between(0, 23) & m.between(0, 59) & sec.between(0, 59)
    return h.where(ok), m.where(ok), sec.where(ok)


def placeholder_row(label: str, h: pd.Series, m: pd.Series,
                    s: pd.Series) -> dict:
    good = h.notna()
    n = int(good.sum())
    exact_mid = int((good & (h == 0) & (m == 0) & (s == 0)).sum())
    hour0 = int((good & (h == 0)).sum())
    on_hour = int((good & (m == 0) & (s == 0)).sum())
    on_hour_not0 = int((good & (h != 0) & (m == 0) & (s == 0)).sum())
    not0 = int((good & (h != 0)).sum())
    nz_sec = int((good & (s != 0)).sum())
    return {
        "field": label,
        "n": n,
        "n_exact_midnight": exact_mid,
        "share_exact_midnight": exact_mid / n if n else np.nan,
        "n_hour_00": hour0,
        "share_of_hour_00_exact_midnight": exact_mid / hour0 if hour0
        else np.nan,
        "share_on_the_hour_all": on_hour / n if n else np.nan,
        "share_on_the_hour_hours_01_23": on_hour_not0 / not0 if not0
        else np.nan,
        "n_unparseable": int((~good).sum()),
        "n_nonzero_seconds": nz_sec,
        "share_nonzero_seconds": nz_sec / n if n else np.nan,
    }


def report_placeholder_rows(rep: Report, rows: list[dict]) -> None:
    for r in rows:
        rep.add(f"  {r['field']}: n = {r['n']:,}; exactly 00:00:00 = "
                f"{r['n_exact_midnight']:,} ({r['share_exact_midnight']:.4f}); "
                f"of hour-00 times exactly 00:00:00 = "
                f"{r['share_of_hour_00_exact_midnight']:.3f}; on the hour, "
                f"hours 01-23 = {r['share_on_the_hour_hours_01_23']:.3f}; "
                f"non-zero seconds = {r['share_nonzero_seconds']:.4f}"
                + (f"; unparseable {r['n_unparseable']:,}"
                   if r["n_unparseable"] else ""))


def minutes_from_clock_strings(s: pd.Series) -> pd.Series:
    s = s.astype("string")
    parts = s.str.extract(r"^\s*(\d{1,2}):(\d{2})", expand=True)
    h = pd.to_numeric(parts[0], errors="coerce")
    m = pd.to_numeric(parts[1], errors="coerce")
    mins = h * 60 + m
    return mins.where(h.between(0, 23) & m.between(0, 59))


class Report:
    def __init__(self) -> None:
        self.lines: list[str] = []
        self.failed: list[str] = []

    def add(self, text: str = "") -> None:
        log(text) if text else print(flush=True)
        self.lines.append(text)

    def section(self, title: str) -> None:
        self.add("")
        self.add(title)
        self.add("-" * len(title))

    def write(self, path: Path) -> None:
        path.write_text("\n".join(self.lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# MIMIC-IV
# ---------------------------------------------------------------------------

def load_admissions(root: Path) -> pd.DataFrame:
    p = find_table(root, "hosp/admissions.csv.gz")
    log(f"reading {p}")
    adm = read_table(p, ["subject_id", "hadm_id", "admittime",
                         "admission_type"])
    adm["admittime"] = pd.to_datetime(adm["admittime"], format="ISO8601",
                                      errors="coerce")
    return adm


def check_m1_weekday(adm: pd.DataFrame, rep: Report, out: Path) -> dict:
    rep.section("M1  MIMIC-IV day of week (admissions.admittime)")
    types = adm["admission_type"].value_counts(dropna=False)
    rep.add("admission_type counts: "
            + "; ".join(f"{k}={v:,}" for k, v in types.items()))
    present = [t for t in SCHEDULED_TYPES if t in types.index]
    absent = [t for t in SCHEDULED_TYPES if t not in types.index]
    if absent:
        rep.add(f"NOTE: scheduled types absent from this release: {absent}")
    if not present:
        raise ValueError("no scheduled admission types found")

    valid = adm.dropna(subset=["admittime"])
    rep.add(f"admissions with parseable admittime: {len(valid):,} "
            f"of {len(adm):,}")
    sched_all = valid[valid["admission_type"].isin(present)]
    # One shift per patient: admissions of the same patient are not
    # independent, so the primary statistic uses each patient's first
    # scheduled admission; all scheduled admissions are reported alongside.
    sched = first_per_subject(sched_all)

    def weekday_counts(frame: pd.DataFrame) -> pd.Series:
        return frame["admittime"].dt.dayofweek.value_counts() \
            .reindex(range(7), fill_value=0)

    wd_sched = weekday_counts(sched)
    wd_sched_all = weekday_counts(sched_all)
    wd_all = weekday_counts(first_per_subject(valid))

    table = pd.DataFrame({"weekday": WEEKDAYS,
                          "n_scheduled_first_per_patient": wd_sched.values,
                          "n_scheduled_all": wd_sched_all.values,
                          "n_all_first_per_patient": wd_all.values})
    for t in present:
        table[f"n_{t.lower().replace(' ', '_')}_first_per_patient"] = \
            weekday_counts(first_per_subject(
                sched_all[sched_all["admission_type"] == t])).values
    table.to_csv(out / "m1_weekday.csv", index=False)

    n = int(wd_sched.sum())
    k = int(wd_sched[5] + wd_sched[6])
    s, lo, hi = share_ci(k, n)
    chi, p = gof(wd_sched.values, np.ones(7))
    scrambled = 2 / 7
    rep.add(f"scheduled admissions ({', '.join(present)}): "
            f"{len(sched_all):,} admissions in {n:,} patients")
    rep.add(f"PRIMARY (first scheduled admission per patient): weekend share "
            f"= {s:.4f} (95% CI {lo:.4f}-{hi:.4f}); expectation if weekday is "
            f"scrambled = {scrambled:.4f}; expectation is "
            f"{'inside' if lo <= scrambled <= hi else 'outside'} the CI")
    rep.add(f"chi-square vs uniform weekdays: {chi:.2f} on 6 df, p = {p:.3g}")
    e = n / 7
    resid = (wd_sched.values - e) / np.sqrt(e * (1 - 1 / 7))
    rep.add("adjusted residuals by weekday (|z| > 1.96 flags a single day): "
            + "; ".join(f"{d} {z:+.2f}" for d, z in zip(WEEKDAYS, resid)))
    for t in present:
        col = table[f"n_{t.lower().replace(' ', '_')}_first_per_patient"]
        st, slo, shi = share_ci(int(col.iloc[5] + col.iloc[6]), int(col.sum()))
        rep.add(f"  {t}: n = {int(col.sum()):,} patients, weekend share "
                f"{st:.4f} (95% CI {slo:.4f}-{shi:.4f})")
    sa_all = (wd_sched_all[5] + wd_sched_all[6]) / wd_sched_all.sum()
    rep.add(f"secondary, all scheduled admissions (patients counted more than "
            f"once; no CI): weekend share {sa_all:.4f}")
    ka = int(wd_all[5] + wd_all[6])
    sa, alo, ahi = share_ci(ka, int(wd_all.sum()))
    rep.add(f"first admission of any type per patient: weekend share "
            f"{sa:.4f} (95% CI {alo:.4f}-{ahi:.4f})")
    return {"weekday_share": (wd_sched / n).values, "n": n}


def check_m2_season(root: Path, adm: pd.DataFrame, chunksize: int,
                    rep: Report, out: Path) -> dict:
    rep.section("M2  MIMIC-IV season (influenza admissions by month)")
    p = find_table(root, "hosp/diagnoses_icd.csv.gz")
    log(f"reading {p}")
    flu_hadm: set = set()
    n_rows = 0
    for chunk in read_table(p, ["hadm_id", "icd_code", "icd_version"],
                            dtype={"icd_code": "string"},
                            chunksize=chunksize):
        n_rows += len(chunk)
        code = chunk["icd_code"].str.strip().str.upper()
        ver = pd.to_numeric(chunk["icd_version"], errors="coerce")
        is_flu = (((ver == 9) & code.str.startswith(FLU_ICD9))
                  | ((ver == 10) & code.str.startswith(FLU_ICD10)))
        flu_hadm.update(chunk.loc[is_flu.fillna(False), "hadm_id"].tolist())
    rep.add(f"diagnoses_icd rows scanned: {n_rows:,}; admissions with an "
            f"influenza code: {len(flu_hadm):,}")
    if not flu_hadm:
        raise ValueError("no influenza-coded admissions found")

    valid = adm.dropna(subset=["admittime"])
    flu_all = valid[valid["hadm_id"].isin(flu_hadm)]
    flu = first_per_subject(flu_all)

    def month_counts(frame: pd.DataFrame) -> pd.Series:
        return frame["admittime"].dt.month.value_counts() \
            .reindex(range(1, 13), fill_value=0)

    m_flu = month_counts(flu)
    m_flu_all = month_counts(flu_all)
    m_all = month_counts(first_per_subject(valid))
    expected = DAYS_PER_MONTH / DAYS_PER_MONTH.sum()
    pd.DataFrame({"month": MONTHS,
                  "n_influenza_first_per_patient": m_flu.values,
                  "n_influenza_all": m_flu_all.values,
                  "n_all_first_per_patient": m_all.values,
                  "expected_share_if_scrambled": expected}) \
        .to_csv(out / "m2_month.csv", index=False)

    n = int(m_flu.sum())
    k = int(m_flu[list(WINTER)].sum())
    s, lo, hi = share_ci(k, n)
    exp_w = float(expected[[m - 1 for m in WINTER]].sum())
    chi, pv = gof(m_flu.values, expected)
    rep.add(f"influenza admissions with parseable admittime: "
            f"{len(flu_all):,} admissions in {n:,} patients")
    rep.add(f"PRIMARY (first influenza admission per patient): Dec-Mar share "
            f"= {s:.4f} (95% CI {lo:.4f}-{hi:.4f}); expectation if season is "
            f"scrambled = {exp_w:.4f}; expectation is "
            f"{'inside' if lo <= exp_w <= hi else 'outside'} the CI")
    rep.add(f"chi-square vs days-per-month expectation: {chi:.2f} on 11 df, "
            f"p = {pv:.3g}")
    sw_all = m_flu_all[list(WINTER)].sum() / m_flu_all.sum()
    rep.add(f"secondary, all influenza admissions (patients counted more than "
            f"once; no CI): Dec-Mar share {sw_all:.4f}")
    return {"month_share": (m_flu / n).values, "expected": expected, "n": n}


def check_m3_clock(root: Path, adm: pd.DataFrame, chunksize: int,
                   max_rows: int, rep: Report, out: Path) -> dict:
    rep.section("M3  MIMIC-IV time of day")
    p = find_table(root, "hosp/labevents.csv.gz")
    header = pd.read_csv(p, nrows=0).columns.tolist()
    has_spec = "specimen_id" in header
    cols = ["charttime"] + (["specimen_id"] if has_spec else [])
    log(f"reading {p} ({', '.join(cols)})")
    counts = np.zeros(24)
    on_hour = np.zeros(24)
    bins10 = np.zeros(BINS)
    n_rows = n_bad = n_sec = 0
    spec_ids, spec_min, spec_sec0 = [], [], []
    for i, chunk in enumerate(read_table(p, cols,
                                         dtype={"charttime": "string"},
                                         chunksize=chunksize)):
        if max_rows and n_rows >= max_rows:
            break
        if max_rows:
            chunk = chunk.iloc[: max_rows - n_rows]
        n_rows += len(chunk)
        h, m, s = clock_parts_datetime(chunk["charttime"])
        n_bad += int(h.isna().sum())
        good = h.notna()
        counts += np.bincount(h[good].astype(int), minlength=24)[:24]
        top = good & (m == 0) & (s == 0)
        on_hour += np.bincount(h[top].astype(int), minlength=24)[:24]
        n_sec += int((good & (s != 0)).sum())
        real = good & ~(top & (h == 0))  # exclude exact-midnight placeholders
        idx = (h[real] * 60 + m[real]).astype(int) // BIN_MIN
        bins10 += np.bincount(idx, minlength=BINS)[:BINS]
        if has_spec:
            sp = pd.to_numeric(chunk["specimen_id"], errors="coerce")
            ok = good & sp.notna()
            d = pd.DataFrame({
                "s": sp[ok].astype("int64").to_numpy(),
                "m": (h[ok] * 60 + m[ok]).astype("int16").to_numpy(),
                "z": (s[ok] == 0).to_numpy()}).drop_duplicates("s")
            spec_ids.append(d["s"].to_numpy())
            spec_min.append(d["m"].to_numpy())
            spec_sec0.append(d["z"].to_numpy())
        if (i + 1) % 10 == 0:
            log(f"  ...{n_rows:,} labevents rows")
    lab = hour_profile(counts)
    rep.add(f"labevents rows read: {n_rows:,}"
            + (" (capped by --max-lab-rows)" if max_rows else "")
            + f"; unparseable charttime: {n_bad:,}; with non-zero seconds: "
              f"{n_sec:,}")
    rep.add(f"lab specimen times: peak hour {lab['peak_hour']:02d}:00, "
            f"busiest 3-hour window {lab['peak_window']} holds "
            f"{lab['peak_window_share']:.3f} of draws (0.125 if flat); "
            f"peak/trough ratio {lab['peak_trough_ratio']:.2f}")
    n_lab = counts.sum()
    rep.add(f"lab charttime exactly 00:00:00: {int(on_hour[0]):,} "
            f"({on_hour[0] / n_lab:.4f} of all; 1/1440 = 0.0007 if minutes "
            f"were uniform); of the {int(counts[0]):,} hour-00 times, "
            f"{on_hour[0] / counts[0]:.3f} are exactly 00:00:00 "
            f"(1/60 = 0.017 if uniform within the hour)")
    other = np.delete(np.arange(24), 0)
    rep.add(f"lab charttime exactly on the hour, hours 01-23: "
            f"{on_hour[other].sum() / counts[other].sum():.3f} "
            f"(1/60 = 0.017 if uniform)")
    if has_spec and spec_ids:
        sid = np.concatenate(spec_ids)
        smin = np.concatenate(spec_min)
        sz = np.concatenate(spec_sec0)
        _, first = np.unique(sid, return_index=True)
        smin, sz = smin[first].astype(int), sz[first]
        sh = smin // 60
        n_sp = len(smin)
        sp_counts = np.bincount(sh, minlength=24)[:24]
        sp_on = np.bincount(sh[sz & (smin % 60 == 0)], minlength=24)[:24]
        rep.add(f"laboratory specimens (distinct specimen_id; time of the "
                f"first result row): {n_sp:,}; result rows per specimen "
                f"{counts.sum() / n_sp:.2f}")
        rep.add(f"specimens exactly 00:00:00: {int(sp_on[0]):,} "
                f"({sp_on[0] / n_sp:.4f}); of the {int(sp_counts[0]):,} "
                f"hour-00 specimens, {sp_on[0] / sp_counts[0]:.3f} are "
                f"exactly 00:00:00; specimens exactly on the hour, hours "
                f"01-23: {sp_on[other].sum() / sp_counts[other].sum():.3f}")
        sph = hour_profile(sp_counts)
        rep.add(f"specimen times: peak hour {sph['peak_hour']:02d}:00, "
                f"busiest 3-hour window {sph['peak_window']} holds "
                f"{sph['peak_window_share']:.3f}")
        pd.DataFrame({"hour": range(24), "n_specimens": sp_counts,
                      "n_specimens_on_the_hour": sp_on}) \
            .to_csv(out / "m3_specimen_hour.csv", index=False)
    elif not has_spec:
        rep.add("NOTE: labevents has no specimen_id column; specimen-level "
                "counts not computed")

    valid = adm.dropna(subset=["admittime"])
    sched = valid[valid["admission_type"].isin(SCHEDULED_TYPES)]
    adm_counts = sched["admittime"].dt.hour.value_counts() \
        .reindex(range(24), fill_value=0).values
    sch = hour_profile(adm_counts)
    rep.add(f"scheduled admission times: n = {sch['n']:,}, peak hour "
            f"{sch['peak_hour']:02d}:00, busiest 3-hour window "
            f"{sch['peak_window']} holds {sch['peak_window_share']:.3f} "
            f"(0.125 if flat)")
    pd.DataFrame({"hour": range(24), "n_lab_charttime": counts.astype(int),
                  "n_lab_on_the_hour": on_hour.astype(int),
                  "n_scheduled_admittime": adm_counts}) \
        .to_csv(out / "m3_hour.csv", index=False)
    return {"lab_share": counts / counts.sum(), "lab_bins10": bins10}


def check_m4_concurrency(root: Path, rep: Report, out: Path) -> dict:
    rep.section("M4  MIMIC-IV concurrency (ICU census from shifted times)")
    p = find_table(root, "icu/icustays.csv.gz")
    log(f"reading {p}")
    icu = read_table(p, ["subject_id", "stay_id", "intime", "outtime"])
    for c in ("intime", "outtime"):
        icu[c] = pd.to_datetime(icu[c], format="ISO8601", errors="coerce")
    n0 = len(icu)
    icu = icu.dropna(subset=["intime", "outtime"])
    icu = icu[icu["outtime"] > icu["intime"]]
    rep.add(f"ICU stays: {n0:,}; with valid intime < outtime: {len(icu):,}")

    epoch = pd.Timestamp("1970-01-01")
    t_in = ((icu["intime"] - epoch) / pd.Timedelta(hours=1)).to_numpy()
    t_out = ((icu["outtime"] - epoch) / pd.Timedelta(hours=1)).to_numpy()
    times = np.concatenate([t_in, t_out])
    delta = np.concatenate([np.ones(len(t_in)), -np.ones(len(t_out))])
    order = np.lexsort((delta, times))  # exits before entries at equal times
    times, delta = times[order], delta[order]
    level = np.cumsum(delta)[:-1]
    dur = np.diff(times)
    span_h = float(dur.sum())
    stay_h = float((t_out - t_in).sum())
    mean_shifted = float((level * dur).sum() / span_h)
    zero_share = float(dur[level == 0].sum() / span_h)
    srt = np.argsort(level, kind="mergesort")
    cumw = np.cumsum(dur[srt]) / span_h
    q = {qq: float(level[srt][np.searchsorted(cumw, qq)])
         for qq in (0.5, 0.9, 0.99)}
    first = icu["intime"].min()
    last = icu["outtime"].max()
    rep.add(f"shifted calendar span: {first:%Y-%m-%d} to {last:%Y-%m-%d} "
            f"({span_h / 24 / 365.25:.1f} years)")
    rep.add(f"census from shifted times (time-weighted): mean {mean_shifted:.2f}, "
            f"median {q[0.5]:.0f}, 90th pct {q[0.9]:.0f}, 99th pct {q[0.99]:.0f},"
            f" max {int(level.max())}; share of time with census 0: "
            f"{zero_share:.3f}")

    pp = find_table(root, "hosp/patients.csv.gz")
    pts = read_table(pp, ["subject_id", "anchor_year_group"])
    groups = pts.loc[pts["subject_id"].isin(icu["subject_id"]),
                     "anchor_year_group"].dropna().astype(str)
    years = groups.str.findall(r"\d{4}").explode().dropna().astype(int)
    real_mean = float("nan")
    if len(years):
        y0, y1 = int(years.min()), int(years.max())
        n_years = y1 - y0 + 1
        real_mean = stay_h / (n_years * 365.25 * 24)
        rep.add(f"anchor_year_group of ICU patients spans {y0}-{y1} "
                f"({n_years} years): {sorted(groups.unique())}")
        rep.add(f"implied mean real census, if total ICU stay-time "
                f"({stay_h / 24:,.0f} stay-days) is spread over {n_years} "
                f"calendar years: {real_mean:.1f} (approximate)")
        rep.add(f"ratio implied real / shifted mean census: "
                f"{real_mean / mean_shifted:.1f}")
    else:
        rep.add("NOTE: anchor_year_group not parseable; implied real census "
                "not computed")

    max_level = int(level.max())
    dist = np.bincount(level.astype(int), weights=dur,
                       minlength=max_level + 1) / span_h
    pd.DataFrame({"census_level": np.arange(len(dist)),
                  "time_share": dist}).to_csv(out / "m4_census_shifted.csv",
                                              index=False)
    pd.DataFrame([{"n_stays": len(icu), "span_years": span_h / 24 / 365.25,
                   "mean_census_shifted": mean_shifted,
                   "median_census_shifted": q[0.5],
                   "p90_census_shifted": q[0.9],
                   "p99_census_shifted": q[0.99],
                   "max_census_shifted": max_level,
                   "time_share_census_zero": zero_share,
                   "implied_mean_real_census": real_mean}]) \
        .to_csv(out / "m4_census_summary.csv", index=False)
    return {"mean_shifted": mean_shifted, "real_mean": real_mean}


def check_m5_placeholders(root: Path, adm: pd.DataFrame, rep: Report,
                          out: Path) -> None:
    rep.section("M5  MIMIC-IV placeholder times (exactly 00:00:00, on the hour)")
    valid = adm.dropna(subset=["admittime"])
    rows = []

    def from_dt(label: str, t: pd.Series) -> dict:
        return placeholder_row(label, t.dt.hour, t.dt.minute, t.dt.second)

    rows.append(from_dt("admittime, all admission types", valid["admittime"]))
    for t, grp in valid.groupby("admission_type", sort=True):
        rows.append(from_dt(f"admittime, {t}", grp["admittime"]))

    p = find_table(root, "icu/icustays.csv.gz")
    icu = read_table(p, ["intime", "outtime"],
                     dtype={"intime": "string", "outtime": "string"})
    for c in ("intime", "outtime"):
        rows.append(placeholder_row(f"icustays.{c}",
                                    *clock_parts_datetime(icu[c])))
    pd.DataFrame(rows).to_csv(out / "m5_placeholders.csv", index=False)
    rep.add("with genuine minute-resolution times: exactly 00:00:00 = 0.0007 "
            "of all times; of hour-00 times 0.017; on the hour 0.017")
    report_placeholder_rows(rep, rows)


# ---------------------------------------------------------------------------
# eICU-CRD
# ---------------------------------------------------------------------------

def check_e1_fields(root: Path, rep: Report, out: Path) -> None:
    rep.section("E1  eICU-CRD calendar and clock fields (all tables)")
    files = sorted(list(root.glob("*.csv.gz")) + list(root.glob("*.csv")))
    if not files:
        raise FileNotFoundError(f"no .csv or .csv.gz tables under {root}")
    rows = []
    for f in files:
        cols = pd.read_csv(f, nrows=0).columns.tolist()
        for c in cols:
            if CALENDAR_PATTERN.search(c):
                rows.append({"table": f.name, "column": c})
    pd.DataFrame(rows, columns=["table", "column"]) \
        .to_csv(out / "e1_calendar_fields.csv", index=False)
    rep.add(f"tables scanned: {len(files)}")
    if rows:
        for r in rows:
            rep.add(f"  {r['table']}: {r['column']}")
    else:
        rep.add("  no column names matched")


def load_eicu_patient(root: Path) -> pd.DataFrame:
    p = find_table(root, "patient.csv.gz")
    log(f"reading {p}")
    required = ["patientunitstayid", "hospitalid", "unitadmittime24",
                "unitdischargeoffset", "hospitaldischargeyear"]
    header = pd.read_csv(p, nrows=0).columns.tolist()
    optional = [c for c in ("unitdischargetime24", "hospitaladmittime24",
                            "hospitaldischargetime24",
                            "hospitaldischargestatus", "unitdischargestatus",
                            "patienthealthsystemstayid", "uniquepid",
                            "unitvisitnumber", "apacheadmissiondx")
                if c in header]
    pat = read_table(p, required + optional,
                     dtype={c: "string" for c in ["unitadmittime24"]
                            + optional
                            if c != "patienthealthsystemstayid"})
    pat["admit_min"] = minutes_from_clock_strings(pat["unitadmittime24"])
    return pat


def check_e4_placeholders(pat: pd.DataFrame, rep: Report, out: Path) -> None:
    rep.section("E4  eICU-CRD placeholder times (exactly 00:00:00, on the hour)")
    cols = [c for c in pat.columns if c.endswith("time24")]
    rows = [placeholder_row(f"patient.{c}", *clock_parts_time24(pat[c]))
            for c in cols]
    pd.DataFrame(rows).to_csv(out / "e4_placeholders.csv", index=False)
    rep.add("with genuine minute-resolution times: exactly 00:00:00 = 0.0007 "
            "of all times; of hour-00 times 0.017; on the hour 0.017")
    report_placeholder_rows(rep, rows)


def circular_lag(profile: np.ndarray, reference: np.ndarray,
                 bin_min: int = BIN_MIN) -> tuple:
    """Shift (hours) at which the reference, moved later by that shift, best
    matches the profile (Pearson r over circular shifts). Positive = the
    profile sits later on the recorded clock than the reference."""
    profile = np.asarray(profile, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if np.isnan(profile).any() or np.isnan(reference).any():
        return float("nan"), float("nan")
    p = profile - profile.mean()
    r = reference - reference.mean()
    denom = np.linalg.norm(p) * np.linalg.norm(r)
    if denom == 0:
        return float("nan"), float("nan")
    n = len(p)
    best_l, best_c = 0, -np.inf
    for lag in range(-(n // 2), n // 2):
        c = float(np.dot(p, np.roll(r, lag)) / denom)
        if c > best_c:
            best_l, best_c = lag, c
    return best_l * bin_min / 60, best_c


def cosinor_phase(profile: np.ndarray, t_hours: np.ndarray) -> tuple:
    """First-harmonic (24-hour) cosinor fit by least squares. Returns the
    acrophase (clock hour of the fitted peak, 0-24) and the amplitude."""
    y = np.asarray(profile, dtype=float)
    t = np.asarray(t_hours, dtype=float)
    ok = ~np.isnan(y)
    if ok.sum() < 6:
        return float("nan"), float("nan")
    w = 2 * np.pi * t[ok] / 24
    X = np.column_stack([np.ones(ok.sum()), np.cos(w), np.sin(w)])
    beta, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    a, b = beta[1], beta[2]
    return (float(np.mod(np.arctan2(b, a) * 24 / (2 * np.pi), 24)),
            float(np.hypot(a, b)))


def phase_lag(profile: np.ndarray, reference: np.ndarray,
              t_hours: np.ndarray) -> float:
    """Acrophase of the profile minus acrophase of the reference, wrapped
    to (-12, 12] hours. Positive = the profile peaks later."""
    pa, _ = cosinor_phase(profile, t_hours)
    pr, _ = cosinor_phase(reference, t_hours)
    if np.isnan(pa) or np.isnan(pr):
        return float("nan")
    d = (pa - pr) % 24
    return float(d - 24 if d > 12 else d)


HOUR_CENTRES = np.arange(24) + 0.5
BIN_CENTRES = (np.arange(BINS) + 0.5) * BIN_MIN / 60


def boot_median_ci(x: np.ndarray, rng, n_boot: int = 2000) -> tuple:
    x = np.asarray(x, dtype=float)
    if len(x) < 2:
        return float("nan"), float("nan")
    meds = np.median(rng.choice(x, size=(n_boot, len(x)), replace=True),
                     axis=1)
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def public_hospital_table(df: pd.DataFrame) -> pd.DataFrame:
    """Hospital-level rows for output. eICU-CRD hospitalid is dropped unless
    --keep-hospital-ids is given, and rows are sorted so that the original
    hospital order cannot be recovered."""
    if KEEP_HOSPITAL_IDS or "hospitalid" not in df.columns:
        return df
    keys = [c for c in ("signal", "region", "lag_hours", "r") if c in df.columns]
    return df.drop(columns=["hospitalid"]).sort_values(keys) \
        .reset_index(drop=True)


def check_e5_timezone(root: Path, e2: dict, m3: dict | None, min_labs: int,
                      ref_region: str, rep: Report, out: Path) -> dict:
    rep.section("E5  eICU-CRD clock reference: local time or UTC?")
    p = find_table(root, "hospital.csv.gz")
    hosp = read_table(p, ["hospitalid", "region"])
    rc = hosp["region"].value_counts(dropna=False)
    rep.add("hospital.region values: "
            + "; ".join(f"{k}={v}" for k, v in rc.items()))
    region_of = hosp.set_index("hospitalid")["region"]
    ids = e2["hosp_ids"]
    lab, adm = e2["lab_mat"], e2["adm_mat"]
    labels = pd.Series(ids).map(region_of).astype("string") \
        .fillna("(missing)").to_numpy()
    n_lab = lab.sum(axis=1)
    keep = n_lab >= min_labs
    rep.add(f"hospitals with laboratory results: {int((n_lab > 0).sum())}; "
            f"with >= {min_labs:,} results (analysed): {int(keep.sum())}")

    present = pd.Series(labels[keep]).value_counts()
    if ref_region not in present.index:
        new_ref = str(present.index[0])
        rep.add(f"NOTE: reference region '{ref_region}' not found among "
                f"analysed hospitals; using '{new_ref}' (most hospitals)")
        ref_region = new_ref
    ref_mask = keep & (labels == ref_region)
    ref_total = lab[ref_mask].sum(axis=0)

    rows = []
    for i in np.where(keep)[0]:
        ref_i = ref_total - lab[i] if ref_mask[i] else ref_total
        lag, r = circular_lag(lab[i] / lab[i].sum(), ref_i / ref_i.sum())
        rows.append({"hospitalid": ids[i], "region": labels[i],
                     "n_lab_results": int(n_lab[i]), "lag_hours": lag,
                     "r": r})
    df = pd.DataFrame(rows)
    public_hospital_table(df).to_csv(out / "e5_hospital_lags.csv",
                                     index=False)

    order = [ref_region] + sorted(r for r in present.index
                                  if r != ref_region)
    rng = np.random.default_rng(BOOT_SEED)
    rep.add(f"lag of each hospital's laboratory-result clock profile "
            f"(10-min bins) against the pooled {ref_region} profile "
            f"(leave-one-out for {ref_region} hospitals); positive = later "
            f"on the recorded clock")
    rep.add("expected if *time24 is UTC (relative to the Northeast, Eastern "
            "zone): Midwest 0 to +1 h, South 0 to +1 h, West +2 to +3 h "
            "(approximate: each region spans more than one zone); if local "
            "time: about 0 in every region")
    summary = []
    for reg in order:
        x = df.loc[df["region"] == reg, "lag_hours"].dropna().to_numpy()
        rr = df.loc[df["region"] == reg, "r"].dropna().to_numpy()
        lo, hi = boot_median_ci(x, rng)
        pooled = lab[keep & (labels == reg)].sum(axis=0)
        p_lag, p_r = (circular_lag(pooled / pooled.sum(),
                                   ref_total / ref_total.sum())
                      if reg != ref_region else (0.0, 1.0))
        adm_reg = adm[labels == reg].sum(axis=0)
        adm_ref = adm[labels == ref_region].sum(axis=0)
        a_lag, _ = (circular_lag(adm_reg / adm_reg.sum(),
                                 adm_ref / adm_ref.sum())
                    if reg != ref_region and adm_reg.sum() > 0
                    else (0.0, 1.0))
        summary.append({"region": reg, "n_hospitals": len(x),
                        "median_lag_hours": float(np.median(x)) if len(x)
                        else np.nan,
                        "q1": float(np.percentile(x, 25)) if len(x)
                        else np.nan,
                        "q3": float(np.percentile(x, 75)) if len(x)
                        else np.nan,
                        "median_ci_lo": lo, "median_ci_hi": hi,
                        "median_r": float(np.median(rr)) if len(rr)
                        else np.nan,
                        "pooled_lab_lag_hours": p_lag,
                        "pooled_unit_admission_lag_hours": a_lag})
        rep.add(f"  {reg}: {len(x)} hospitals; median lag "
                f"{summary[-1]['median_lag_hours']:+.2f} h "
                f"(IQR {summary[-1]['q1']:+.2f} to {summary[-1]['q3']:+.2f}; "
                f"bootstrap 95% CI of median {lo:+.2f} to {hi:+.2f}); "
                f"median r {summary[-1]['median_r']:.2f}; pooled lab lag "
                f"{p_lag:+.2f} h; pooled unit-admission lag {a_lag:+.2f} h")
    pd.DataFrame(summary).to_csv(out / "e5_region_summary.csv", index=False)
    outl = df[df["lag_hours"].abs() > OUTLIER_HOURS]
    rep.add(f"hospitals with |lag| > {OUTLIER_HOURS:.0f} h (outlying "
            f"profiles): {len(outl)}"
            + (" - " + "; ".join(
                (f"{r.hospitalid} " if KEEP_HOSPITAL_IDS else "")
                + f"({r.region}, {r.lag_hours:+.2f} h, r {r.r:.2f}, "
                f"{r.n_lab_results:,} results)"
                for r in outl.itertuples())
               if len(outl) else ""))

    mimic_lag = float("nan")
    if m3 is not None and m3.get("lab_bins10") is not None \
            and m3["lab_bins10"].sum() > 0:
        mb = m3["lab_bins10"]
        mimic_lag, m_r = circular_lag(ref_total / ref_total.sum(),
                                      mb / mb.sum())
        rep.add(f"eICU {ref_region} laboratory-result profile against "
                f"MIMIC-IV specimen-time profile (Boston local time, exact "
                f"midnight excluded): lag {mimic_lag:+.2f} h (r {m_r:.2f}); "
                f"expected about +5 to +6 h if eICU is UTC (offset plus "
                f"turnaround), about +1 h (turnaround) if local")
    else:
        rep.add("MIMIC-IV comparison not run (M3 profile unavailable)")

    hourly = {reg: lab[keep & (labels == reg)].sum(axis=0)
              .reshape(24, BINS // 24).sum(axis=1) for reg in order}
    return {"order": order, "hourly": hourly, "lags": df,
            "ref": ref_region, "mimic_lag": mimic_lag}


# ---------------------------------------------------------------------------
# E6 - physiological anchor for the eICU clock
# ---------------------------------------------------------------------------

def stay_deviation_profiles(sums: np.ndarray, cnts: np.ndarray):
    """Per-stay hourly means minus the stay's own mean (mean of its hourly
    means), for stays covering >= MIN_STAY_HOURS of the 24 clock hours.
    Returns (deviation matrix with NaN for uncovered hours, kept mask)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        m = sums / cnts
    covered = (cnts > 0).sum(axis=1)
    keep = covered >= MIN_STAY_HOURS
    mk = m[keep]
    with np.errstate(invalid="ignore"):
        dev = mk - np.nanmean(mk, axis=1, keepdims=True)
    return dev, keep


def mean_profile(dev: np.ndarray) -> np.ndarray:
    if len(dev) == 0:
        return np.full(24, np.nan)
    with np.errstate(invalid="ignore"):
        return np.nanmean(dev, axis=0)


def verify_items(mroot: Path, rep: Report) -> dict:
    p = find_table(mroot, "icu/d_items.csv.gz")
    di = read_table(p, ["itemid", "label"])
    labels = di.set_index("itemid")["label"].astype(str)
    usable = {}
    for name, (iid, expected, _, _) in VITALS.items():
        got = labels.get(iid, "(absent)")
        ok = got.strip().lower() == expected.lower()
        rep.add(f"d_items {iid}: '{got}' (expected '{expected}') -> "
                + ("used" if ok else "NOT USED"))
        if ok:
            usable[name] = iid
    return usable


def mimic_vital_profiles(mroot: Path, usable: dict, chunksize: int,
                         rep: Report) -> dict:
    icu = read_table(find_table(mroot, "icu/icustays.csv.gz"), ["stay_id"])
    stay_index = pd.Index(icu["stay_id"].unique())
    n = len(stay_index)
    acc = {v: [np.zeros(n * 24), np.zeros(n * 24)] for v in usable}
    minute_sum = {v: 0.0 for v in usable}
    item_to_name = {iid: name for name, iid in usable.items()}
    p = find_table(mroot, "icu/chartevents.csv.gz")
    log(f"reading {p} (stay_id, itemid, charttime, valuenum)")
    n_rows = 0
    for i, chunk in enumerate(read_table(
            p, ["stay_id", "itemid", "charttime", "valuenum"],
            dtype={"charttime": "string"}, chunksize=chunksize)):
        n_rows += len(chunk)
        chunk = chunk[chunk["itemid"].isin(list(item_to_name))]
        if len(chunk):
            h, mi, se = clock_parts_datetime(chunk["charttime"])
            frac = (mi + se / 60).to_numpy()
            val = pd.to_numeric(chunk["valuenum"], errors="coerce")
            si = stay_index.get_indexer(chunk["stay_id"])
            for iid, name in item_to_name.items():
                lo, hi = VITALS[name][3]
                ok = ((chunk["itemid"] == iid) & h.notna() & val.between(lo, hi)
                      & (si >= 0)).to_numpy()
                idx = si[ok] * 24 + h[ok].astype(int).to_numpy()
                acc[name][0] += np.bincount(idx, weights=val[ok].to_numpy(),
                                            minlength=n * 24)
                acc[name][1] += np.bincount(idx, minlength=n * 24)
                minute_sum[name] += float(np.nansum(frac[ok]))
        if (i + 1) % 10 == 0:
            log(f"  ...{n_rows:,} chartevents rows")
    rep.add(f"MIMIC-IV chartevents rows scanned: {n_rows:,}")
    out = {}
    for name, (s, c) in acc.items():
        dev, keep = stay_deviation_profiles(s.reshape(n, 24), c.reshape(n, 24))
        out[name] = mean_profile(dev)
        out[f"{name}__minute"] = minute_sum[name] / max(c.sum(), 1)
        rep.add(f"MIMIC-IV {name}: {int(c.sum()):,} values; stays covering "
                f">= {MIN_STAY_HOURS} clock hours: {int(keep.sum()):,}; "
                f"profile nadir at {int(np.nanargmin(out[name])):02d}:00, "
                f"peak at {int(np.nanargmax(out[name])):02d}:00 "
                f"(Boston local time); mean minute within the clock hour "
                f"{out[f'{name}__minute']:.1f}")
    return out


def eicu_vital_profiles(eroot: Path, pat: pd.DataFrame, chunksize: int,
                        rep: Report) -> dict:
    stay_index = pd.Index(pat["patientunitstayid"])
    admit = pat["admit_min"].to_numpy()
    n = len(stay_index)
    acc = {v: [np.zeros(n * 24), np.zeros(n * 24)] for v in VITALS}
    minute_sum = {v: 0.0 for v in VITALS}
    cols = ["patientunitstayid", "observationoffset"] + \
        [VITALS[v][2] for v in VITALS]
    p = find_table(eroot, "vitalPeriodic.csv.gz", "vitalperiodic.csv.gz")
    log(f"reading {p} ({', '.join(cols)})")
    n_rows = 0
    for i, chunk in enumerate(read_table(p, cols, chunksize=chunksize)):
        n_rows += len(chunk)
        si = stay_index.get_indexer(chunk["patientunitstayid"])
        off = pd.to_numeric(chunk["observationoffset"], errors="coerce") \
            .to_numpy()
        a = np.where(si >= 0, admit[np.clip(si, 0, None)], np.nan)
        clock = np.mod(a + off, 1440)
        hour = np.floor(clock / 60)
        base_ok = (si >= 0) & ~np.isnan(hour)
        for name in VITALS:
            col, (lo, hi) = VITALS[name][2], VITALS[name][3]
            val = pd.to_numeric(chunk[col], errors="coerce").to_numpy()
            ok = base_ok & (val >= lo) & (val <= hi)
            idx = si[ok] * 24 + hour[ok].astype(int)
            acc[name][0] += np.bincount(idx, weights=val[ok],
                                        minlength=n * 24)
            acc[name][1] += np.bincount(idx, minlength=n * 24)
            minute_sum[name] += float(np.sum(np.mod(clock[ok], 60)))
        if (i + 1) % 10 == 0:
            log(f"  ...{n_rows:,} vitalPeriodic rows")
    rep.add(f"eICU-CRD vitalPeriodic rows scanned: {n_rows:,}")
    out = {}
    for name, (s, c) in acc.items():
        dev, keep = stay_deviation_profiles(s.reshape(n, 24), c.reshape(n, 24))
        out[name] = (dev, keep)
        out[f"{name}__minute"] = minute_sum[name] / max(c.sum(), 1)
        rep.add(f"eICU-CRD {name}: {int(c.sum()):,} values; stays covering "
                f">= {MIN_STAY_HOURS} clock hours: {int(keep.sum()):,}; mean "
                f"minute within the clock hour "
                f"{out[f'{name}__minute']:.1f}")
    return out


def check_e6_physiology(mroot: Path, eroot: Path, pat: pd.DataFrame,
                        chunksize: int, ref_region: str, rep: Report,
                        out: Path) -> dict:
    rep.section("E6  eICU-CRD clock reference: physiological anchor "
                "(MIMIC-IV, Boston local time)")
    usable = verify_items(mroot, rep)
    if "heart_rate" not in usable:
        raise ValueError("heart-rate itemid not verified in d_items")
    mimic = mimic_vital_profiles(mroot, usable, chunksize, rep)
    eicu = eicu_vital_profiles(eroot, pat, chunksize, rep)

    hosp = read_table(find_table(eroot, "hospital.csv.gz"),
                      ["hospitalid", "region"])
    region_of = hosp.set_index("hospitalid")["region"]
    stay_hosp = pat["hospitalid"].to_numpy()
    stay_region = pd.Series(stay_hosp).map(region_of).astype("string") \
        .fillna("(missing)").to_numpy()
    rep.add("lag of eICU profile (recorded clock) against MIMIC-IV profile "
            "(local); positive = eICU later on its recorded clock")
    rep.add("expected if eICU *time24 is UTC: Northeast +4 to +5 h, "
            "Midwest +4 to +6 h, South +4 to +6 h, West +6 to +8 h "
            "(contiguous states; ranges include standard and daylight saving "
            "time, which cannot be resolved without dates); if a common US "
            "Eastern clock: Northeast 0, Midwest and South 0 to +1, West +2 "
            "to +3; if local time: about 0 everywhere")

    rng = np.random.default_rng(BOOT_SEED)
    hosp_rows, reg_rows, prof_rows = [], [], []
    for h in range(24):
        row = {"hour": h}
        for name in [v for v in VITALS if v in mimic]:
            row[f"mimic_{name}"] = mimic[name][h]
        prof_rows.append(row)

    rep.add("cosinor phase lags are corrected for the difference in mean "
            "minute within the clock hour between the two databases (values "
            "binned by hour are labelled at the bin centre on both sides)")
    for name in [v for v in VITALS if v in mimic]:
        dev, keep = eicu[name]
        adj = (eicu[f"{name}__minute"] - mimic[f"{name}__minute"]) / 60
        rep.add(f"{name}: bin-timing correction {adj:+.2f} h (eICU mean "
                f"minute {eicu[f'{name}__minute']:.1f}, MIMIC-IV "
                f"{mimic[f'{name}__minute']:.1f})")
        regions_k = stay_region[keep]
        hosp_k = stay_hosp[keep]
        ref_prof = mimic[name]
        # per hospital
        for hid in pd.unique(hosp_k):
            sel = hosp_k == hid
            if sel.sum() < MIN_HOSPITAL_STAYS:
                continue
            hp = mean_profile(dev[sel])
            lag, r = circular_lag(hp, ref_prof, 60)
            hosp_rows.append({"signal": name, "hospitalid": hid,
                              "region": regions_k[sel][0],
                              "n_stays": int(sel.sum()),
                              "lag_hours": lag, "r": r,
                              "phase_lag_hours": phase_lag(
                                  hp, ref_prof, HOUR_CENTRES) + adj,
                              "amplitude": cosinor_phase(
                                  hp, HOUR_CENTRES)[1]})
        hdf = pd.DataFrame([r for r in hosp_rows if r["signal"] == name])
        regs = [ref_region] + sorted(r for r in pd.unique(regions_k)
                                     if r != ref_region)
        rep.add(f"{name}:")
        for reg in regs:
            sel = regions_k == reg
            if not sel.any():
                continue
            prof = mean_profile(dev[sel])
            for h in range(24):
                prof_rows[h][f"eicu_{name}_{reg}"] = prof[h]
            p_lag, p_r = circular_lag(prof, ref_prof, 60)
            x = hdf.loc[hdf["region"] == reg, "lag_hours"].dropna() \
                .to_numpy() if len(hdf) else np.array([])
            lo, hi = boot_median_ci(x, rng)
            med = float(np.median(x)) if len(x) else float("nan")
            px = hdf.loc[hdf["region"] == reg, "phase_lag_hours"].dropna() \
                .to_numpy() if len(hdf) else np.array([])
            plo, phi = boot_median_ci(px, rng)
            pmed = float(np.median(px)) if len(px) else float("nan")
            pp_lag = phase_lag(prof, ref_prof, HOUR_CENTRES) + adj
            reg_rows.append({"signal": name, "region": reg,
                             "n_stays": int(sel.sum()), "n_hospitals": len(x),
                             "pooled_lag_hours": p_lag, "pooled_r": p_r,
                             "median_hospital_lag_hours": med,
                             "median_ci_lo": lo, "median_ci_hi": hi,
                             "pooled_phase_lag_hours": pp_lag,
                             "median_hospital_phase_lag_hours": pmed,
                             "phase_ci_lo": plo, "phase_ci_hi": phi})
            rep.add(f"  {reg}: {int(sel.sum()):,} stays; pooled lag "
                    f"{p_lag:+.0f} h (r {p_r:.2f}); {len(x)} hospitals with "
                    f">= {MIN_HOSPITAL_STAYS} stays: median lag {med:+.1f} h "
                    f"(bootstrap 95% CI {lo:+.1f} to {hi:+.1f}); cosinor "
                    f"phase lag: pooled {pp_lag:+.2f} h, hospital median "
                    f"{pmed:+.2f} h (bootstrap 95% CI {plo:+.2f} to "
                    f"{phi:+.2f})")
        all_prof = mean_profile(dev)
        a_lag, a_r = circular_lag(all_prof, ref_prof, 60)
        rep.add(f"  all eICU stays: cosinor phase lag "
                f"{phase_lag(all_prof, ref_prof, HOUR_CENTRES) + adj:+.2f} h; "
                f"MIMIC-IV acrophase "
                f"{cosinor_phase(ref_prof, HOUR_CENTRES)[0]:.2f} h")
        rep.add(f"  all eICU stays: pooled lag {a_lag:+.0f} h (r {a_r:.2f}); "
                f"eICU nadir at {int(np.nanargmin(all_prof)):02d}:00 and "
                f"MIMIC-IV nadir at {int(np.nanargmin(ref_prof)):02d}:00 on "
                f"their recorded clocks")

    public_hospital_table(pd.DataFrame(hosp_rows)).to_csv(
        out / "e6_hospital_lags.csv", index=False)
    pd.DataFrame(reg_rows).to_csv(out / "e6_region_summary.csv", index=False)
    pd.DataFrame(prof_rows).to_csv(out / "e6_profiles.csv", index=False)
    return {"profiles": pd.DataFrame(prof_rows),
            "hosp": pd.DataFrame(hosp_rows), "ref": ref_region,
            "eicu_hr": eicu["heart_rate"]}


def make_figure3(e6: dict, path: Path) -> None:
    prof = e6["profiles"]
    hosp = e6["hosp"]
    hr = hosp[hosp["signal"] == "heart_rate"]
    order = [r for r in ["Northeast", "Midwest", "South", "West"]
             if f"eicu_heart_rate_{r}" in prof.columns]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), facecolor="white",
                           gridspec_kw={"width_ratios": [1.5, 1]})
    a = ax[0]
    a.plot(prof["hour"], prof["mimic_heart_rate"], color=INK, linewidth=2,
           linestyle="--", label="MIMIC-IV (Boston, local time)")
    for colour, reg in zip(REGION_COLOURS, order):
        a.plot(prof["hour"], prof[f"eicu_heart_rate_{reg}"], color=colour,
               linewidth=2, label=f"eICU {reg} (recorded clock)")
    a.axhline(0, color=GRID, linewidth=1)
    a.set_xticks(range(0, 24, 3), [f"{h:02d}" for h in range(0, 24, 3)])
    a.set_xlim(-0.5, 23.5)
    a.set_xlabel("Clock hour", fontsize=8, color=INK_2)
    a.set_ylabel("Heart rate minus stay mean (beats/min)", fontsize=8,
                 color=INK_2)
    a.legend(frameon=True, facecolor="white", edgecolor=GRID, framealpha=0.9,
             fontsize=7.5, loc="best")
    style_axis(a, "A  Heart-rate day-night profile")

    b = ax[1]
    n_out = strip_panel(b, hr, order, -1.5, 10.5)
    b.set_ylabel("Shift against MIMIC-IV local time (hours)", fontsize=8,
                 color=INK_2)
    key = "0 if local; UTC: NE +4–5, MW/S +4–6, W +6–8"
    if n_out:
        key += f"; ▲▼ beyond axis: {n_out}"
    style_axis(b, "B  Per-hospital heart-rate profile shift", key)
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI, facecolor="white")
    plt.close(fig)


def clock_bins(h: pd.Series, m: pd.Series, s: pd.Series) -> np.ndarray:
    """10-minute bin counts of clock times, excluding exact-midnight values
    (potential date-only placeholders) and unparseable times."""
    good = h.notna() & ~((h == 0) & (m == 0) & (s == 0))
    idx = ((h[good] * 60 + m[good]).astype(int) // BIN_MIN).to_numpy()
    return np.bincount(idx, minlength=BINS)[:BINS].astype(float)


def check_e7_discharge(mroot: Path, eroot: Path, pat: pd.DataFrame,
                       ref_region: str, rep: Report, out: Path) -> dict:
    rep.section("E7  eICU-CRD clock reference: administrative anchor "
                "(discharge clock times)")
    need = ["hospitaldischargestatus", "unitdischargestatus",
            "patienthealthsystemstayid", "hospitaldischargetime24",
            "unitdischargetime24"]
    missing = [c for c in need if c not in pat.columns]
    if missing:
        raise KeyError(f"eICU patient columns missing: {missing}")
    adm = read_table(find_table(mroot, "hosp/admissions.csv.gz"),
                     ["hadm_id", "dischtime", "hospital_expire_flag"],
                     dtype={"dischtime": "string"})
    icu = read_table(find_table(mroot, "icu/icustays.csv.gz"),
                     ["hadm_id", "outtime"], dtype={"outtime": "string"})
    flag = pd.to_numeric(adm["hospital_expire_flag"], errors="coerce")
    alive = adm[(flag == 0) & adm["hadm_id"].isin(icu["hadm_id"])]
    icu_alive = icu[icu["hadm_id"].isin(alive["hadm_id"])]
    mimic = {"hospital_discharge": clock_bins(
                 *clock_parts_datetime(alive["dischtime"])),
             "unit_discharge": clock_bins(
                 *clock_parts_datetime(icu_alive["outtime"]))}
    rep.add(f"MIMIC-IV live hospital discharges after an ICU stay: "
            f"{int(mimic['hospital_discharge'].sum()):,}; their ICU "
            f"discharges: {int(mimic['unit_discharge'].sum()):,} "
            f"(exact-midnight values excluded)")
    for sig, v in mimic.items():
        hp = hour_profile(v.reshape(24, BINS // 24).sum(axis=1))
        rep.add(f"  MIMIC-IV {sig}: peak hour {hp['peak_hour']:02d}:00, "
                f"busiest 3-hour window {hp['peak_window']} holds "
                f"{hp['peak_window_share']:.3f}; cosinor acrophase "
                f"{cosinor_phase(v / v.sum(), BIN_CENTRES)[0]:.2f} h "
                f"(Boston local time)")

    hosp = read_table(find_table(eroot, "hospital.csv.gz"),
                      ["hospitalid", "region"])
    region_of = hosp.set_index("hospitalid")["region"]
    hstat = pat["hospitaldischargestatus"].astype("string").str.strip()
    ustat = pat["unitdischargestatus"].astype("string").str.strip()
    e_hosp = pat[hstat == "Alive"].drop_duplicates("patienthealthsystemstayid")
    e_unit = pat[(ustat == "Alive") & (hstat == "Alive")]
    frames = {"hospital_discharge": (e_hosp, "hospitaldischargetime24"),
              "unit_discharge": (e_unit, "unitdischargetime24")}
    rep.add(f"eICU-CRD live hospital discharges (one per "
            f"patienthealthsystemstayid): {len(e_hosp):,}; live unit "
            f"discharges of patients discharged alive: {len(e_unit):,}")
    rep.add("lag of eICU discharge-time profile (recorded clock, 10-min bins) "
            "against MIMIC-IV (local); positive = eICU later on its recorded "
            "clock")
    rep.add("expected if eICU *time24 is UTC: Northeast +4 to +5 h, Midwest "
            "+4 to +6 h, South +4 to +6 h, West +6 to +8 h; if a common US "
            "Eastern clock: Northeast 0, Midwest and South 0 to +1, West +2 "
            "to +3; if local time: about 0 everywhere")

    rng = np.random.default_rng(BOOT_SEED)
    hosp_rows, reg_rows = [], []
    prof_rows = [{"hour": h} for h in range(24)]
    for sig, (frame, col) in frames.items():
        ref = mimic[sig] / mimic[sig].sum()
        for h in range(24):
            prof_rows[h][f"mimic_{sig}"] = \
                ref.reshape(24, BINS // 24).sum(axis=1)[h]
        hh, mm, ss = clock_parts_time24(frame[col])
        regions = frame["hospitalid"].map(region_of).astype("string") \
            .fillna("(missing)").to_numpy()
        hids = frame["hospitalid"].to_numpy()
        for hid in pd.unique(hids):
            sel = hids == hid
            v = clock_bins(hh[sel], mm[sel], ss[sel])
            if v.sum() < MIN_HOSPITAL_DISCHARGES:
                continue
            lag, r = circular_lag(v / v.sum(), ref)
            hosp_rows.append({"signal": sig, "hospitalid": hid,
                              "region": regions[sel][0],
                              "n_discharges": int(v.sum()),
                              "lag_hours": lag, "r": r,
                              "phase_lag_hours": phase_lag(
                                  v / v.sum(), ref, BIN_CENTRES)})
        hdf = pd.DataFrame([r for r in hosp_rows if r["signal"] == sig])
        regs = [ref_region] + sorted(r for r in pd.unique(regions)
                                     if r != ref_region)
        rep.add(f"{sig}:")
        for reg in regs:
            sel = regions == reg
            if not sel.any():
                continue
            v = clock_bins(hh[sel], mm[sel], ss[sel])
            share = v / v.sum()
            for h in range(24):
                prof_rows[h][f"eicu_{sig}_{reg}"] = \
                    share.reshape(24, BINS // 24).sum(axis=1)[h]
            p_lag, p_r = circular_lag(share, ref)
            pp_lag = phase_lag(share, ref, BIN_CENTRES)
            x = hdf.loc[hdf["region"] == reg, "lag_hours"].dropna() \
                .to_numpy() if len(hdf) else np.array([])
            px = hdf.loc[hdf["region"] == reg, "phase_lag_hours"].dropna() \
                .to_numpy() if len(hdf) else np.array([])
            lo, hi = boot_median_ci(x, rng)
            plo, phi = boot_median_ci(px, rng)
            med = float(np.median(x)) if len(x) else float("nan")
            pmed = float(np.median(px)) if len(px) else float("nan")
            reg_rows.append({"signal": sig, "region": reg,
                             "n_discharges": int(v.sum()),
                             "n_hospitals": len(x),
                             "pooled_lag_hours": p_lag, "pooled_r": p_r,
                             "median_hospital_lag_hours": med,
                             "median_ci_lo": lo, "median_ci_hi": hi,
                             "pooled_phase_lag_hours": pp_lag,
                             "median_hospital_phase_lag_hours": pmed,
                             "phase_ci_lo": plo, "phase_ci_hi": phi})
            rep.add(f"  {reg}: {int(v.sum()):,} discharges; pooled lag "
                    f"{p_lag:+.2f} h (r {p_r:.2f}); {len(x)} hospitals with "
                    f">= {MIN_HOSPITAL_DISCHARGES} discharges: median lag "
                    f"{med:+.2f} h (bootstrap 95% CI {lo:+.2f} to "
                    f"{hi:+.2f}); cosinor phase lag: pooled {pp_lag:+.2f} h, "
                    f"hospital median {pmed:+.2f} h (bootstrap 95% CI "
                    f"{plo:+.2f} to {phi:+.2f})")

    public_hospital_table(pd.DataFrame(hosp_rows)).to_csv(
        out / "e7_hospital_lags.csv", index=False)
    pd.DataFrame(reg_rows).to_csv(out / "e7_region_summary.csv", index=False)
    pd.DataFrame(prof_rows).to_csv(out / "e7_profiles.csv", index=False)
    return {"profiles": pd.DataFrame(prof_rows),
            "hosp": pd.DataFrame(hosp_rows), "ref": ref_region}


def make_figure4(e7: dict, path: Path) -> None:
    prof = e7["profiles"]
    hosp = e7["hosp"]
    hd = hosp[hosp["signal"] == "hospital_discharge"]
    order = [r for r in ["Northeast", "Midwest", "South", "West"]
             if f"eicu_hospital_discharge_{r}" in prof.columns]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), facecolor="white",
                           gridspec_kw={"width_ratios": [1.5, 1]})
    a = ax[0]
    a.plot(prof["hour"], prof["mimic_hospital_discharge"], color=INK,
           linewidth=2, linestyle="--", label="MIMIC-IV (Boston, local time)")
    for colour, reg in zip(REGION_COLOURS, order):
        a.plot(prof["hour"], prof[f"eicu_hospital_discharge_{reg}"],
               color=colour, linewidth=2, label=f"eICU {reg} (recorded clock)")
    a.set_xticks(range(0, 24, 3), [f"{h:02d}" for h in range(0, 24, 3)])
    a.set_xlim(-0.5, 23.5)
    headroom(a, prof.drop(columns="hour").to_numpy())
    a.set_xlabel("Clock hour", fontsize=8, color=INK_2)
    a.set_ylabel("Share of live hospital discharges", fontsize=8,
                 color=INK_2)
    a.legend(frameon=True, facecolor="white", edgecolor=GRID, framealpha=0.9,
             fontsize=7.5, loc="best")
    style_axis(a, "A  Hospital discharge clock time")

    b = ax[1]
    n_out = strip_panel(b, hd, order, -3.5, 10.5)
    b.set_ylabel("Shift against MIMIC-IV local time (hours)", fontsize=8,
                 color=INK_2)
    key = "0 if local; UTC: NE +4\u20135, MW/S +4\u20136, W +6\u20138"
    if n_out:
        key += f"; \u25b2\u25bc beyond axis: {n_out}"
    style_axis(b, "B  Per-hospital discharge-time shift", key)
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI, facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------------
# E8 - worked example: a published clock-time analysis on the recorded clock
#      and after conversion to local time
# ---------------------------------------------------------------------------

def stay_cosinor(dev: np.ndarray) -> tuple:
    """First-harmonic (24-hour) cosinor fit to each row of hourly values
    (NaN = clock hour not covered), hours labelled at their centres.
    Returns per-row acrophase (clock hour of the fitted peak, 0-24),
    amplitude and the zero-amplitude F-test p-value."""
    from scipy import stats
    t = 2 * np.pi * HOUR_CENTRES / 24
    c, s = np.cos(t), np.sin(t)
    w = (~np.isnan(dev)).astype(float)
    y = np.where(w > 0, dev, 0.0)
    n = w.sum(axis=1)
    a = np.empty((len(dev), 3, 3))
    a[:, 0, 0] = n
    a[:, 0, 1] = a[:, 1, 0] = (w * c).sum(axis=1)
    a[:, 0, 2] = a[:, 2, 0] = (w * s).sum(axis=1)
    a[:, 1, 1] = (w * c * c).sum(axis=1)
    a[:, 2, 2] = (w * s * s).sum(axis=1)
    a[:, 1, 2] = a[:, 2, 1] = (w * c * s).sum(axis=1)
    b = np.stack([y.sum(axis=1), (y * c).sum(axis=1), (y * s).sum(axis=1)],
                 axis=1)
    beta = np.linalg.solve(a, b[..., None])[..., 0]
    fit = beta[:, [0]] + beta[:, [1]] * c + beta[:, [2]] * s
    ss_res = (w * (y - fit) ** 2).sum(axis=1)
    ybar = y.sum(axis=1) / n
    ss_tot = (w * (y - ybar[:, None]) ** 2).sum(axis=1)
    df2 = n - 3
    with np.errstate(divide="ignore", invalid="ignore"):
        f = ((ss_tot - ss_res) / 2) / (ss_res / df2)
    p = stats.f.sf(f, 2, df2)
    acro = np.mod(np.arctan2(beta[:, 2], beta[:, 1]) * 24 / (2 * np.pi), 24)
    amp = np.hypot(beta[:, 1], beta[:, 2])
    return acro, amp, p


def logit_or(y: np.ndarray, x: np.ndarray, covar, groups: np.ndarray) -> tuple:
    """Odds ratio for binary x (logistic regression, optional covariate),
    95% CI from hospital-clustered standard errors."""
    import statsmodels.api as sm
    cols = [np.ones(len(y)), x.astype(float)]
    if covar is not None:
        cols.append(covar.astype(float))
    try:
        r = sm.Logit(y.astype(float), np.column_stack(cols)).fit(
            disp=0, cov_type="cluster", cov_kwds={"groups": groups})
        lo, hi = r.conf_int()[1]
        return float(np.exp(r.params[1])), float(np.exp(lo)), float(np.exp(hi))
    except Exception as exc:  # separation, too few events
        log(f"E8 logistic model failed: {exc}")
        return float("nan"), float("nan"), float("nan")


def as_bool(x) -> np.ndarray:
    """Boolean mask as a plain numpy array; missing (NA) counts as False."""
    if isinstance(x, (pd.Series, pd.Index)):
        return x.fillna(False).to_numpy(dtype=bool)
    return np.asarray(x, dtype=bool)


def circular_mean_hour(h: np.ndarray) -> float:
    ang = np.angle(np.mean(np.exp(1j * 2 * np.pi * np.asarray(h) / 24)))
    return float(np.mod(ang * 24 / (2 * np.pi), 24))


def window_label(k: int) -> str:
    lo, hi = PEAK_WINDOWS[k]
    return f"{lo:02d}:00-{hi - 1:02d}:59"


def check_e8_peak_time(eroot: Path, pat: pd.DataFrame, e6: dict,
                       e7: dict | None, rep: Report, out: Path) -> dict:
    rep.section("E8  Worked example: heart-rate peak time and in-hospital "
                "mortality, recorded clock against local time")
    rep.add("design modelled on Yang et al., Heart Rhythm 2022 (eICU-CRD; "
            "highest in-hospital mortality with heart-rate peak time "
            "12:00-18:00 on the recorded clock): per-stay cosinor peak time "
            "(acrophase) of hourly heart rate, stays with significant "
            "rhythmicity (zero-amplitude F test p < 0.05), peak time in "
            "6-hour windows, logistic regression for in-hospital mortality")
    rep.add("expected if the eICU-CRD clock is UTC: a recorded peak time of "
            "12:00-17:59 corresponds to about 05:00-13:30 local time (regional "
            "offsets 4.5-7 h), so the stays it labels as afternoon peaks have "
            "morning peaks in local time, and the association attached to "
            "12:00-17:59 changes after conversion; if the clock is local, "
            "conversion should make the association worse, not clearer")
    if e6 is None or "eicu_hr" not in e6:
        raise ValueError("E6 heart-rate stay profiles are not available")
    dev, keep = e6["eicu_hr"]
    p = pat.reset_index(drop=True).copy()
    if len(p) != len(keep):
        raise ValueError("patient table and E6 stay profiles do not align")
    p["row"] = np.arange(len(p))
    hosp = read_table(find_table(eroot, "hospital.csv.gz"),
                      ["hospitalid", "region"])
    region_of = hosp.set_index("hospitalid")["region"]

    # one unit stay per patient: first unit stay of the first hospital stay
    if "uniquepid" in p.columns:
        keys = [k for k in ("uniquepid", "hospitaldischargeyear",
                            "patienthealthsystemstayid", "unitvisitnumber",
                            "patientunitstayid") if k in p.columns]
        if "unitvisitnumber" in p.columns:
            p["unitvisitnumber"] = pd.to_numeric(p["unitvisitnumber"],
                                                 errors="coerce")
        first = p.sort_values(keys).drop_duplicates("uniquepid")
        rep.add(f"unit stays: {len(p):,}; patients (uniquepid): "
                f"{len(first):,}; one stay per patient (first unit stay of "
                f"the first hospital stay, by discharge year)")
    else:
        first = p
        rep.add("uniquepid not available: every unit stay counted")

    pos = np.full(len(p), -1)
    pos[np.flatnonzero(keep)] = np.arange(int(keep.sum()))
    f = first[pos[first["row"].to_numpy()] >= 0].copy()
    rep.add(f"with heart rate covering >= {MIN_STAY_HOURS} clock hours: "
            f"{len(f):,}")
    acro, amp, pval = stay_cosinor(dev[pos[f["row"].to_numpy()]])
    f["acro"], f["amp"], f["p"] = acro, amp, pval
    status = f["hospitaldischargestatus"].astype("string").str.strip()
    f["died"] = np.where(as_bool(status.eq("Expired")), 1.0,
                         np.where(as_bool(status.eq("Alive")), 0.0, np.nan))
    f["region"] = f["hospitalid"].map(region_of).astype("string")

    apache_note = "APACHE IV score not available: adjusted models omitted"
    f["apache"] = np.nan
    try:
        ap = find_table(eroot, "apachePatientResult.csv.gz",
                        "apachepatientresult.csv.gz")
        apr = read_table(ap, ["patientunitstayid", "apacheversion",
                              "apachescore"],
                         dtype={"apacheversion": "string"})
        ver = apr["apacheversion"].str.strip().str.upper()
        sc = pd.to_numeric(apr["apachescore"], errors="coerce")
        apr = apr[as_bool(ver.isin(["IV", "4"])) & as_bool(sc >= 0)]
        score = apr.assign(sc=sc).drop_duplicates("patientunitstayid") \
            .set_index("patientunitstayid")["sc"]
        f["apache"] = f["patientunitstayid"].map(score)
        apache_note = ("adjusted models: APACHE IV score (apachePatientResult, "
                       "apacheversion IV, score >= 0)")
    except (FileNotFoundError, KeyError) as exc:
        log(f"E8: {exc}")
    rep.add(apache_note)

    steps = [("rhythmic (p < 0.05)", as_bool(f["p"] < 0.05)),
             ("known in-hospital outcome", as_bool(f["died"].notna())),
             ("hospital region recorded",
              as_bool(f["region"].isin(list(UTC_OFFSET_MID)))),
             ("APACHE IV score available", as_bool(f["apache"].notna()))]
    rep.add(f"  missing in-hospital outcome: {int(f['died'].isna().sum()):,}; "
            f"no hospital region: {int(f['region'].isna().sum()):,}")
    sel = np.ones(len(f), dtype=bool)
    for label, cond in steps:
        if label.startswith("APACHE") and f["apache"].isna().all():
            continue
        sel &= cond
        rep.add(f"  {label}: {int(sel.sum()):,}")
    co = f[sel].copy()
    if len(co) < 100:
        raise ValueError(f"too few stays for E8 ({len(co)})")

    versions = {"recorded clock": {r: 0.0 for r in UTC_OFFSET_MID},
                "local time, UTC offset (midpoint)": dict(UTC_OFFSET_MID),
                "local time, UTC offset - 1 h":
                    {r: v - 1 for r, v in UTC_OFFSET_MID.items()},
                "local time, UTC offset + 1 h":
                    {r: v + 1 for r, v in UTC_OFFSET_MID.items()}}
    if e7 is not None and "hosp" in e7 and len(e7["hosp"]):
        hd = e7["hosp"]
        hd = hd[hd["signal"] == "hospital_discharge"]
        med = hd.groupby("region")["phase_lag_hours"].median()
        if all(r in med.index for r in UTC_OFFSET_MID):
            versions["local time, hospital-discharge anchor offset"] = \
                {r: float(med[r]) for r in UTC_OFFSET_MID}
    rep.add("conversion to local time: peak time minus a regional offset "
            "(hours): " + "; ".join(
                f"{k}: " + ", ".join(f"{r} {v:+.2f}" for r, v in o.items())
                for k, o in versions.items() if k != "recorded clock"))

    cohorts = {"all rhythmic first stays": co}
    if "apacheadmissiondx" in co.columns:
        dx = co["apacheadmissiondx"].astype("string")
        stroke = as_bool(dx.str.strip().str.lower().eq(STROKE_DX.lower()))
        cohorts["stroke admission diagnosis"] = co[stroke]
        other = dx.str.contains(r"stroke|\bCVA\b", case=False, regex=True,
                                na=False)
        rep.add(f"stroke subgroup: apacheadmissiondx exactly '{STROKE_DX}': "
                f"{int(stroke.sum()):,} stays; other diagnoses containing "
                f"'stroke' or 'CVA' (not included): "
                + ("; ".join(f"{k} ({v:,})" for k, v in
                             dx[as_bool(other) & ~stroke].value_counts()
                             .head(5).items()) or "none"))

    win_rows, model_rows = [], []
    k_risk = PEAK_WINDOWS.index(RISK_WINDOW)
    for cname, cdf in cohorts.items():
        if len(cdf) < 100 or cdf["died"].sum() < 10:
            rep.add(f"{cname}: too few stays or deaths for models "
                    f"({len(cdf):,} stays)")
            continue
        y = cdf["died"].to_numpy()
        groups = pd.factorize(cdf["hospitalid"])[0]
        reg = cdf["region"].to_numpy()
        rep.add(f"{cname}: {len(cdf):,} stays in {len(np.unique(groups))} "
                f"hospitals; in-hospital deaths {int(y.sum()):,} "
                f"({y.mean():.3f})")
        win_rec = None
        for vname, off in versions.items():
            shift = np.array([off[r] for r in reg])
            acro_v = np.mod(cdf["acro"].to_numpy() - shift, 24)
            win = np.minimum((acro_v // 6).astype(int), 3)
            if win_rec is None:
                win_rec = win
            for k in range(len(PEAK_WINDOWS)):
                m = win == k
                win_rows.append({"cohort": cname, "clock": vname,
                                 "peak_window": window_label(k),
                                 "n": int(m.sum()),
                                 "deaths": int(y[m].sum()),
                                 "mortality": float(y[m].mean())
                                 if m.any() else float("nan")})
            x = win == k_risk
            u = logit_or(y, x, None, groups)
            a = logit_or(y, x, cdf["apache"].to_numpy(), groups) \
                if cdf["apache"].notna().all() else (float("nan"),) * 3
            moved = float((win != win_rec).mean())
            model_rows.append({"cohort": cname, "clock": vname,
                               "n": len(cdf), "deaths": int(y.sum()),
                               "share_peak_12_18": float(x.mean()),
                               "mortality_peak_12_18": float(y[x].mean())
                               if x.any() else float("nan"),
                               "mortality_other": float(y[~x].mean())
                               if (~x).any() else float("nan"),
                               "share_moved_window": moved,
                               "or_unadjusted": u[0], "or_unadj_lo": u[1],
                               "or_unadj_hi": u[2], "or_adjusted": a[0],
                               "or_adj_lo": a[1], "or_adj_hi": a[2]})
            shares = ", ".join(
                f"{window_label(k)} {np.mean(win == k):.3f}"
                for k in range(len(PEAK_WINDOWS)))
            rep.add(f"  {vname}: peak-time shares {shares}; stays in a "
                    f"different window than on the recorded clock: {moved:.3f}")
            rep.add(f"    peak {window_label(k_risk)} against all other "
                    f"windows: mortality {model_rows[-1]['mortality_peak_12_18']:.3f} "
                    f"vs {model_rows[-1]['mortality_other']:.3f}; OR {u[0]:.2f} "
                    f"(95% CI {u[1]:.2f}-{u[2]:.2f}); APACHE IV-adjusted OR "
                    f"{a[0]:.2f} (95% CI {a[1]:.2f}-{a[2]:.2f})")
        acro_loc = np.mod(cdf["acro"].to_numpy() - np.array(
            [UTC_OFFSET_MID[r] for r in reg]), 24)
        rep.add(f"  circular mean peak time: recorded clock "
                f"{circular_mean_hour(cdf['acro'].to_numpy()):.1f} h; local "
                f"time (midpoint offset) {circular_mean_hour(acro_loc):.1f} h")
    pd.DataFrame(win_rows).to_csv(out / "e8_peak_time_windows.csv",
                                  index=False)
    pd.DataFrame(model_rows).to_csv(out / "e8_models.csv", index=False)
    return {"windows": pd.DataFrame(win_rows),
            "models": pd.DataFrame(model_rows)}


def check_e2_clock(root: Path, pat: pd.DataFrame, chunksize: int,
                   max_rows: int, rep: Report, out: Path) -> dict:
    rep.section("E2  eICU-CRD time of day")
    yrs = pat["hospitaldischargeyear"].value_counts(dropna=False).sort_index()
    rep.add("hospitaldischargeyear: "
            + "; ".join(f"{k}={v:,}" for k, v in yrs.items()))
    ok = pat["admit_min"].notna()
    rep.add(f"unit stays: {len(pat):,}; parseable unitadmittime24: "
            f"{int(ok.sum()):,}")
    adm_counts = np.bincount((pat.loc[ok, "admit_min"] // 60).astype(int),
                             minlength=24)[:24]
    ua = hour_profile(adm_counts)
    rep.add(f"unit admission clock time: peak hour {ua['peak_hour']:02d}:00, "
            f"busiest 3-hour window {ua['peak_window']} holds "
            f"{ua['peak_window_share']:.3f} (0.125 if flat)")

    admit_map = pat.loc[ok].set_index("patientunitstayid")["admit_min"]
    hosp_codes, hosp_ids = pd.factorize(pat["hospitalid"])
    n_h = len(hosp_ids)
    code_map = pd.Series(hosp_codes, index=pat["patientunitstayid"])
    adm_mat = np.zeros((n_h, BINS))
    a_ok = ok & (hosp_codes >= 0)
    np.add.at(adm_mat, (hosp_codes[a_ok.to_numpy()],
                        (pat.loc[a_ok, "admit_min"] // BIN_MIN).astype(int)
                        .to_numpy()), 1)

    p = find_table(root, "lab.csv.gz")
    log(f"reading {p} (patientunitstayid, labresultoffset)")
    counts = np.zeros(24)
    lab_mat = np.zeros(n_h * BINS)
    n_rows = n_unmatched = 0
    for i, chunk in enumerate(read_table(p, ["patientunitstayid",
                                             "labresultoffset"],
                                         chunksize=chunksize)):
        if max_rows and n_rows >= max_rows:
            break
        if max_rows:
            chunk = chunk.iloc[: max_rows - n_rows]
        n_rows += len(chunk)
        a = chunk["patientunitstayid"].map(admit_map)
        off = pd.to_numeric(chunk["labresultoffset"], errors="coerce")
        clock = (a + off) % 1440
        n_unmatched += int(clock.isna().sum())
        counts += np.bincount((clock.dropna() // 60).astype(int),
                              minlength=24)[:24]
        code = chunk["patientunitstayid"].map(code_map)
        use = clock.notna() & code.notna() & (code >= 0)
        flat = (code[use].astype(int) * BINS
                + (clock[use] // BIN_MIN).astype(int)).to_numpy()
        lab_mat += np.bincount(flat, minlength=n_h * BINS)[: n_h * BINS]
        if (i + 1) % 10 == 0:
            log(f"  ...{n_rows:,} lab rows")
    lab = hour_profile(counts)
    rep.add(f"lab rows read: {n_rows:,}"
            + (" (capped by --max-lab-rows)" if max_rows else "")
            + f"; without a usable admission clock time or offset: "
              f"{n_unmatched:,}")
    rep.add(f"lab result clock time (unitadmittime24 + labresultoffset): "
            f"peak hour {lab['peak_hour']:02d}:00, busiest 3-hour window "
            f"{lab['peak_window']} holds {lab['peak_window_share']:.3f} "
            f"(0.125 if flat); peak/trough ratio "
            f"{lab['peak_trough_ratio']:.2f}")
    rep.add("note: labresultoffset is the result time, not the draw time")
    pd.DataFrame({"hour": range(24), "n_unit_admission": adm_counts,
                  "n_lab_result": counts.astype(int)}) \
        .to_csv(out / "e2_hour.csv", index=False)
    return {"lab_share": counts / counts.sum(),
            "hosp_ids": np.asarray(hosp_ids),
            "lab_mat": lab_mat.reshape(n_h, BINS),
            "adm_mat": adm_mat}


def check_e3_offset_census(pat: pd.DataFrame, rep: Report, out: Path) -> dict:
    rep.section("E3  eICU-CRD concurrency (census built from offsets)")
    off = pd.to_numeric(pat["unitdischargeoffset"], errors="coerce").dropna()
    off = np.sort(off.to_numpy())
    hours = np.arange(CENSUS_HOURS_MAX + 1)
    present = len(off) - np.searchsorted(off, hours * 60, side="right")
    monotone = bool(np.all(np.diff(present) <= 0))
    rep.add(f"unit stays with unitdischargeoffset: {len(off):,}; offsets are "
            f"minutes from each stay's own unit admission (offset 0)")
    rep.add("stays 'present' at hour h after admission: "
            + "; ".join(f"h={h}: {present[h]:,}" for h in (0, 24, 72, 168)))
    rep.add(f"non-increasing in h: {'yes' if monotone else 'NO'} - the "
            f"offset 'census' is the length-of-stay survival curve")
    pd.DataFrame({"hours_since_unit_admission": hours,
                  "stays_present": present}) \
        .to_csv(out / "e3_offset_census.csv", index=False)
    return {"hours": hours, "present": present}


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def style_axis(ax, title: str, key: str | None = None) -> None:
    if key:
        title = f"{title}\n    - - -  {key}"
    ax.set_title(title, loc="left", fontsize=9.5, color=INK)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=8)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def headroom(ax, *values) -> None:
    top = max(float(np.nanmax(np.atleast_1d(v))) for v in values)
    ax.set_ylim(0, top * 1.12)


def ref_line(ax, y, data) -> None:
    ax.axhline(y, color=INK_2, linestyle="--", linewidth=1)
    headroom(ax, data, y)


def ref_steps(ax, x, y, data) -> None:
    ax.plot(x, y, color=INK_2, linestyle="--", linewidth=1,
            drawstyle="steps-mid")
    headroom(ax, data, y)


def make_figure(res: dict, path: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.2), facecolor="white")
    ax = axes.ravel()
    blank = []

    if "m1" in res:
        a = ax[0]
        a.bar(range(7), res["m1"]["weekday_share"], color=SERIES, width=0.7)
        ref_line(a, 1 / 7, res["m1"]["weekday_share"])
        a.set_xticks(range(7), WEEKDAYS)
        a.set_ylabel("Share of admissions", fontsize=8, color=INK_2)
        style_axis(a, f"A  MIMIC-IV scheduled admission by weekday\n"
                      f"    first per patient, n={res['m1']['n']:,}",
                   "if weekday scrambled (1/7)")
    else:
        blank.append(0)

    if "m2" in res:
        a = ax[1]
        a.bar(range(12), res["m2"]["month_share"], color=SERIES, width=0.7)
        ref_steps(a, np.arange(12), res["m2"]["expected"],
                  res["m2"]["month_share"])
        a.set_xticks(range(12), [m[0] for m in MONTHS])
        a.set_ylabel("Share of admissions", fontsize=8, color=INK_2)
        style_axis(a, f"B  MIMIC-IV influenza admission by month\n"
                      f"    first per patient, n={res['m2']['n']:,}",
                   "if season scrambled (days per month)")
    else:
        blank.append(1)

    if "m3" in res:
        a = ax[2]
        a.bar(range(24), res["m3"]["lab_share"], color=SERIES, width=0.75)
        ref_line(a, 1 / 24, res["m3"]["lab_share"])
        a.set_xticks(range(0, 24, 3), [f"{h:02d}" for h in range(0, 24, 3)])
        a.set_xlabel("Hour of day", fontsize=8, color=INK_2)
        a.set_ylabel("Share of result rows", fontsize=8, color=INK_2)
        style_axis(a, "C  MIMIC-IV laboratory specimen time\n"
                      "    labevents.charttime", "if clock scrambled (1/24)")
    else:
        blank.append(2)

    if "m4" in res and np.isfinite(res["m4"]["real_mean"]):
        a = ax[3]
        vals = [res["m4"]["mean_shifted"], res["m4"]["real_mean"]]
        a.bar([0, 1], vals, color=SERIES, width=0.55)
        for i, v in enumerate(vals):
            a.text(i, v, f"{v:.1f}", ha="center", va="bottom", fontsize=8,
                   color=INK)
        headroom(a, vals)
        a.set_xticks([0, 1], ["From shifted\ntimes", "Implied by\nreal span"])
        a.set_ylabel("Mean concurrent ICU stays", fontsize=8, color=INK_2)
        style_axis(a, "D  MIMIC-IV ICU census\n"
                      "    time-weighted mean")
    else:
        blank.append(3)

    if "e2" in res:
        a = ax[4]
        a.bar(range(24), res["e2"]["lab_share"], color=SERIES, width=0.75)
        ref_line(a, 1 / 24, res["e2"]["lab_share"])
        a.set_xticks(range(0, 24, 3), [f"{h:02d}" for h in range(0, 24, 3)])
        a.set_xlabel("Hour of day", fontsize=8, color=INK_2)
        a.set_ylabel("Share of results", fontsize=8, color=INK_2)
        style_axis(a, "E  eICU-CRD laboratory result clock time\n"
                      "    unitadmittime24 + labresultoffset",
                   "if clock scrambled (1/24)")
    else:
        blank.append(4)

    if "e3" in res:
        a = ax[5]
        h, n = res["e3"]["hours"], res["e3"]["present"]
        keep = h <= 168
        a.plot(h[keep], n[keep], color=SERIES, linewidth=2)
        a.set_xticks(range(0, 169, 24))
        a.set_xlabel("Hours since each stay's unit admission", fontsize=8,
                     color=INK_2)
        a.set_ylabel("Stays 'present'", fontsize=8, color=INK_2)
        style_axis(a, "F  eICU-CRD 'census' built from offsets\n"
                      "    = length-of-stay survival curve")
    else:
        blank.append(5)

    for i in blank:
        ax[i].axis("off")
        ax[i].text(0.5, 0.5, "check not run or failed", ha="center",
                   va="center", fontsize=9, color=INK_2,
                   transform=ax[i].transAxes)
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI, facecolor="white")
    plt.close(fig)


def strip_panel(b, lags: pd.DataFrame, order: list, lo_ax: float,
                hi_ax: float) -> int:
    """Per-hospital lags by region: jittered points, median bar with a
    labelled value, points beyond [lo_ax, hi_ax] drawn as triangles at the
    axis edge. Returns the number of points beyond the axis."""
    jit = np.random.default_rng(BOOT_SEED)
    n_out = 0
    for i, reg in enumerate(order):
        x = lags.loc[lags["region"] == reg, "lag_hours"].dropna().to_numpy()
        if not len(x):
            continue
        xs = i + jit.uniform(-0.18, 0.18, len(x))
        inside = (x >= lo_ax) & (x <= hi_ax)
        b.scatter(xs[inside], x[inside], s=14, color=SERIES, alpha=0.6,
                  linewidths=0)
        below, above = x < lo_ax, x > hi_ax
        if below.any():
            b.scatter(xs[below], np.full(below.sum(), lo_ax + 0.15), s=22,
                      marker="v", color=INK_2, linewidths=0)
        if above.any():
            b.scatter(xs[above], np.full(above.sum(), hi_ax - 0.15), s=22,
                      marker="^", color=INK_2, linewidths=0)
        n_out += int(below.sum() + above.sum())
        med = float(np.median(x))
        b.plot([i - 0.3, i + 0.3], [med, med], color=INK, linewidth=2)
        b.text(i + 0.33, med, f"{med:+.1f} h", fontsize=7.5, color=INK,
               va="center", ha="left",
               bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
    b.axhline(0, color=INK_2, linestyle="--", linewidth=1)
    b.set_ylim(lo_ax, hi_ax)
    b.set_xticks(range(len(order)), order)
    b.set_xlim(-0.6, len(order) - 0.05)
    return n_out


def make_figure2(e5: dict, path: Path) -> None:
    order = [r for r in e5["order"] if r != "(missing)"][:4]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), facecolor="white",
                           gridspec_kw={"width_ratios": [1.5, 1]})

    a = ax[0]
    peak = 0.0
    for colour, reg in zip(REGION_COLOURS, order):
        h = e5["hourly"][reg]
        share = h / h.sum()
        peak = max(peak, float(share.max()))
        a.plot(range(24), share, color=colour, linewidth=2, label=reg)
    a.set_xlim(-0.5, 23.5)
    a.set_ylim(0, peak * 1.12)
    a.set_xticks(range(0, 24, 3), [f"{h:02d}" for h in range(0, 24, 3)])
    a.set_xlabel("Recorded clock hour (unitadmittime24 + labresultoffset)",
                 fontsize=8, color=INK_2)
    a.set_ylabel("Share of laboratory results", fontsize=8, color=INK_2)
    a.legend(frameon=False, fontsize=8, loc="upper left")
    style_axis(a, "A  eICU-CRD laboratory result clock time by region")

    b = ax[1]
    lags = e5["lags"]
    lo_ax, hi_ax = -AXIS_CLIP_HOURS[0], AXIS_CLIP_HOURS[1]
    n_out = strip_panel(b, lags, order, lo_ax, hi_ax)
    b.set_ylabel(f"Shift against pooled {e5['ref']} (hours)", fontsize=8,
                 color=INK_2)
    key = "0 if local; later zones > 0 if UTC"
    if n_out:
        key += f"; ▲▼ beyond axis: {n_out}"
    style_axis(b, "B  Per-hospital shift of the laboratory-result profile",
               key)
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI, facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_check(rep: Report, name: str, fn, *args):
    try:
        return fn(*args)
    except Exception as exc:  # report and continue with the other checks
        rep.add(f"FAILED {name}: {type(exc).__name__}: {exc}")
        rep.failed.append(name)
        traceback.print_exc()
        return None


def write_manifest(args, paths: list[Path], out: Path, started: str) -> None:
    lines = [
        f"script: {Path(__file__).name}",
        f"started: {started}",
        f"finished: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"python: {platform.python_version()}",
        f"pandas: {pd.__version__}",
        f"numpy: {np.__version__}",
        f"statsmodels: {statsmodels.__version__}",
        f"matplotlib: {matplotlib.__version__}",
        f"arguments: {vars(args)}",
        "input files:",
    ]
    for p in paths:
        st = p.stat()
        lines.append(f"  {p}  {st.st_size:,} bytes  modified "
                     f"{time.strftime('%Y-%m-%d %H:%M', time.localtime(st.st_mtime))}")
    (out / "run_manifest.txt").write_text("\n".join(lines) + "\n",
                                          encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Paper 26: timestamp checks for MIMIC-IV and eICU-CRD")
    ap.add_argument("--mimic-root", type=Path,
                    default=Path("~/physionet.org/files/mimiciv/3.1"))
    ap.add_argument("--eicu-root", type=Path,
                    default=Path("~/physionet.org/files/eicu-crd/2.0"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("~/paper26_results"))
    ap.add_argument("--chunksize", type=int, default=5_000_000)
    ap.add_argument("--max-lab-rows", type=int, default=0,
                    help="stop each laboratory pass after N rows (0 = all)")
    ap.add_argument("--min-hospital-labs", type=int, default=2000,
                    help="E5: minimum laboratory results per hospital")
    ap.add_argument("--reference-region", default="Northeast",
                    help="E5: reference region for the clock-shift test")
    ap.add_argument("--skip-mimic", action="store_true")
    ap.add_argument("--skip-eicu", action="store_true")
    ap.add_argument("--skip-e6", action="store_true",
                    help="skip E6 (physiological anchor; heavy passes)")
    ap.add_argument("--only-e6", action="store_true",
                    help="run only E6 and figure 3")
    ap.add_argument("--skip-e7", action="store_true",
                    help="skip E7 (administrative anchor: discharge times)")
    ap.add_argument("--only-e7", action="store_true",
                    help="run only E7 and figure 4 (light; no large tables)")
    ap.add_argument("--skip-e8", action="store_true",
                    help="skip E8 (worked example; needs E6)")
    ap.add_argument("--keep-hospital-ids", action="store_true",
                    help="write eICU-CRD hospitalid in hospital-level "
                         "outputs (local use only; default omits it)")
    args = ap.parse_args()
    global KEEP_HOSPITAL_IDS
    KEEP_HOSPITAL_IDS = args.keep_hospital_ids

    mroot = args.mimic_root.expanduser()
    eroot = args.eicu_root.expanduser()
    out = args.out_dir.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    started = time.strftime("%Y-%m-%d %H:%M:%S")
    rep = Report()
    only = args.only_e6 or args.only_e7
    rep.add("Paper 26 timestamp checks (v11)"
            + (" - E6 only" if args.only_e6 else "")
            + (" - E7 only" if args.only_e7 else ""))
    rep.add(f"run {started}")
    rep.add(f"MIMIC-IV root: {mroot}" + (" (skipped)" if args.skip_mimic else ""))
    rep.add(f"eICU-CRD root: {eroot}" + (" (skipped)" if args.skip_eicu else ""))

    res: dict = {}
    used: list[Path] = []
    pat = None

    if not args.skip_mimic and not only:
        if not mroot.exists():
            rep.add(f"FAILED MIMIC-IV: root not found: {mroot}")
            rep.failed.append("MIMIC-IV root")
        else:
            for rel in ("hosp/admissions.csv.gz", "hosp/diagnoses_icd.csv.gz",
                        "hosp/labevents.csv.gz", "hosp/patients.csv.gz",
                        "icu/icustays.csv.gz"):
                try:
                    used.append(find_table(mroot, rel))
                except FileNotFoundError as exc:
                    rep.add(f"NOTE: {exc}")
            adm = run_check(rep, "MIMIC-IV admissions", load_admissions, mroot)
            if adm is not None:
                r = run_check(rep, "M1", check_m1_weekday, adm, rep, out)
                if r is not None:
                    res["m1"] = r
                r = run_check(rep, "M2", check_m2_season, mroot, adm,
                              args.chunksize, rep, out)
                if r is not None:
                    res["m2"] = r
                r = run_check(rep, "M3", check_m3_clock, mroot, adm,
                              args.chunksize, args.max_lab_rows, rep, out)
                if r is not None:
                    res["m3"] = r
                run_check(rep, "M5", check_m5_placeholders, mroot, adm, rep,
                          out)
            r = run_check(rep, "M4", check_m4_concurrency, mroot, rep, out)
            if r is not None:
                res["m4"] = r

    if not args.skip_eicu and not only:
        if not eroot.exists():
            rep.add(f"FAILED eICU-CRD: root not found: {eroot}")
            rep.failed.append("eICU-CRD root")
        else:
            for rel in ("patient.csv.gz", "lab.csv.gz", "hospital.csv.gz"):
                try:
                    used.append(find_table(eroot, rel))
                except FileNotFoundError as exc:
                    rep.add(f"NOTE: {exc}")
            run_check(rep, "E1", check_e1_fields, eroot, rep, out)
            pat = run_check(rep, "eICU patient", load_eicu_patient, eroot)
            if pat is not None:
                r = run_check(rep, "E2", check_e2_clock, eroot, pat,
                              args.chunksize, args.max_lab_rows, rep, out)
                if r is not None:
                    res["e2"] = r
                r = run_check(rep, "E3", check_e3_offset_census, pat, rep, out)
                if r is not None:
                    res["e3"] = r
                run_check(rep, "E4", check_e4_placeholders, pat, rep, out)
                if "e2" in res:
                    r = run_check(rep, "E5", check_e5_timezone, eroot,
                                  res["e2"], res.get("m3"),
                                  args.min_hospital_labs,
                                  args.reference_region, rep, out)
                    if r is not None:
                        res["e5"] = r

    run_e6 = (not args.skip_e6) and not args.only_e7 \
        and not (args.skip_mimic or args.skip_eicu)
    if run_e6:
        if not (mroot.exists() and eroot.exists()):
            rep.add("FAILED E6: MIMIC-IV or eICU-CRD root not found")
            rep.failed.append("E6 roots")
        else:
            for root, rel in ((mroot, "icu/chartevents.csv.gz"),
                              (mroot, "icu/d_items.csv.gz"),
                              (mroot, "icu/icustays.csv.gz"),
                              (eroot, "vitalPeriodic.csv.gz"),
                              (eroot, "patient.csv.gz"),
                              (eroot, "hospital.csv.gz")):
                try:
                    p = find_table(root, rel)
                    if p not in used:
                        used.append(p)
                except FileNotFoundError as exc:
                    rep.add(f"NOTE: {exc}")
            if pat is None:
                pat = run_check(rep, "eICU patient", load_eicu_patient, eroot)
            if pat is not None:
                r = run_check(rep, "E6", check_e6_physiology, mroot, eroot,
                              pat, args.chunksize, args.reference_region,
                              rep, out)
                if r is not None:
                    res["e6"] = r

    run_e7 = (not args.skip_e7) and not args.only_e6 \
        and not (args.skip_mimic or args.skip_eicu)
    if run_e7:
        if not (mroot.exists() and eroot.exists()):
            rep.add("FAILED E7: MIMIC-IV or eICU-CRD root not found")
            rep.failed.append("E7 roots")
        else:
            for root, rel in ((mroot, "hosp/admissions.csv.gz"),
                              (mroot, "icu/icustays.csv.gz"),
                              (eroot, "patient.csv.gz"),
                              (eroot, "hospital.csv.gz")):
                try:
                    p = find_table(root, rel)
                    if p not in used:
                        used.append(p)
                except FileNotFoundError as exc:
                    rep.add(f"NOTE: {exc}")
            if pat is None:
                pat = run_check(rep, "eICU patient", load_eicu_patient, eroot)
            if pat is not None:
                r = run_check(rep, "E7", check_e7_discharge, mroot, eroot,
                              pat, args.reference_region, rep, out)
                if r is not None:
                    res["e7"] = r

    if "e6" in res and not args.skip_e8:
        try:
            p = find_table(eroot, "apachePatientResult.csv.gz",
                           "apachepatientresult.csv.gz")
            if p not in used:
                used.append(p)
        except FileNotFoundError as exc:
            rep.add(f"NOTE: {exc}")
        r = run_check(rep, "E8", check_e8_peak_time, eroot, pat, res["e6"],
                      res.get("e7"), rep, out)
        if r is not None:
            res["e8"] = r

    if any(k in res for k in ("m1", "m2", "m3", "m4", "e2", "e3")):
        run_check(rep, "figure", make_figure, res, out / "paper26_fig1.png")
    if "e5" in res:
        run_check(rep, "figure 2", make_figure2, res["e5"],
                  out / "paper26_fig2.png")
    if "e6" in res:
        run_check(rep, "figure 3", make_figure3, res["e6"],
                  out / "paper26_fig3.png")
    if "e7" in res:
        run_check(rep, "figure 4", make_figure4, res["e7"],
                  out / "paper26_fig4.png")

    rep.section("Status")
    rep.add("all checks completed" if not rep.failed
            else f"FAILED: {', '.join(rep.failed)}")
    rep.add(f"outputs in {out}")
    rep.write(out / "paper26_summary.txt")
    write_manifest(args, used, out, started)
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main())
