# /// script
# requires-python = ">=3.11"
# dependencies = ["pandas", "numpy"]
# ///
"""Synthetic MIMIC-IV / eICU-shaped tables for testing paper26_timestamp_checks_v10.py.

v4 adds the missing values found in the real eICU-CRD tables: hospitals with
no region and stays with no hospitaldischargestatus.

v3 adds (for E8): per-stay heart-rate peak time in local time (25% of stays
with a morning peak at 09:00, the rest at 16:00), in-hospital mortality
higher with a morning peak, uniquepid (some patients with two stays),
unitvisitnumber, apacheadmissiondx and an apachePatientResult table.

v2 adds: labevents.specimen_id (results grouped into specimens), admissions
dischtime and hospital_expire_flag, discharge clock times with a midday peak
in local time (MIMIC-IV dischtime and icustays.outtime; eICU-CRD
hospitaldischargetime24 and unitdischargetime24), eICU-CRD discharge status
and patienthealthsystemstayid.

Usage: make_mock_v4.py <scrambled|preserved> <outdir> [utc|local]

variant 'scrambled' : per-patient shift of a random number of whole days
                      (weekday and season scrambled, clock preserved)
variant 'preserved' : per-patient shift of whole weeks closest to whole years
                      (weekday preserved, season approximately preserved)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

variant, out = sys.argv[1], Path(sys.argv[2])
rng = np.random.default_rng(7)
N_PAT = 15000

mroot = out / "mimiciv" / "3.1"
(mroot / "hosp").mkdir(parents=True, exist_ok=True)
(mroot / "icu").mkdir(parents=True, exist_ok=True)

# --- real-time admissions ------------------------------------------------
rows = []
hadm = 20000000
for sid in range(10000000, 10000000 + N_PAT):
    for _ in range(rng.integers(1, 4)):
        hadm += 1
        u = rng.random()
        atype = ("ELECTIVE" if u < 0.15 else
                 "SURGICAL SAME DAY ADMISSION" if u < 0.25 else
                 "URGENT" if u < 0.40 else "EW EMER.")
        day = pd.Timestamp("2008-01-01") + pd.Timedelta(
            days=int(rng.integers(0, 365 * 12)))
        if atype in ("ELECTIVE", "SURGICAL SAME DAY ADMISSION"):
            while day.dayofweek >= 5 and rng.random() < 0.98:
                day += pd.Timedelta(days=1)
            minute = int(rng.integers(5 * 60, 9 * 60))
            if atype == "SURGICAL SAME DAY ADMISSION" and rng.random() < 0.3:
                minute = 0
        else:
            minute = int(rng.integers(0, 1440))
        rows.append((sid, hadm, atype, day + pd.Timedelta(minutes=minute)))
adm = pd.DataFrame(rows, columns=["subject_id", "hadm_id", "admission_type",
                                  "real_admit"])

# --- shift per subject ---------------------------------------------------
subs = adm["subject_id"].unique()
k_years = rng.integers(100, 190, size=len(subs))
if variant == "scrambled":
    shift_days = rng.integers(36500, 69000, size=len(subs))
else:
    shift_days = (np.round(365.2425 * k_years / 7) * 7).astype(int)
shift = pd.Series(shift_days, index=subs)
adm["shift"] = pd.to_timedelta(adm["subject_id"].map(shift), unit="D")
adm["admittime"] = adm["real_admit"] + adm["shift"]
def midday_minute(k, p=0.5):
    """Local discharge minute: share p between 11:00 and 15:00, rest uniform."""
    return np.where(rng.random(k) < p, rng.integers(11 * 60, 15 * 60, k),
                    rng.integers(0, 1440, k))


los_days = rng.integers(1, 10, len(adm))
disch_local = adm["real_admit"].dt.normalize() \
    + pd.to_timedelta(los_days, unit="D") \
    + pd.to_timedelta(midday_minute(len(adm)), unit="m")
adm["dischtime"] = disch_local + adm["shift"]
adm["hospital_expire_flag"] = (rng.random(len(adm)) < 0.05).astype(int)
adm[["subject_id", "hadm_id", "admittime", "dischtime", "admission_type",
     "hospital_expire_flag"]].assign(
    admittime=lambda d: d["admittime"].dt.strftime("%Y-%m-%d %H:%M:%S"),
    dischtime=lambda d: d["dischtime"].dt.strftime("%Y-%m-%d %H:%M:%S"),
).to_csv(mroot / "hosp" / "admissions.csv.gz", index=False)

# --- diagnoses: influenza seasonal in real time ---------------------------
dx = []
for r in adm.itertuples():
    winter = r.real_admit.month in (12, 1, 2, 3)
    p = 0.06 if winter else 0.004
    ver = 9 if r.real_admit.year < 2015 else 10
    if rng.random() < p:
        dx.append((r.hadm_id, "4871" if ver == 9 else "J101", ver))
    dx.append((r.hadm_id, "4019" if ver == 9 else "I10", ver))
pd.DataFrame(dx, columns=["hadm_id", "icd_code", "icd_version"]) \
    .to_csv(mroot / "hosp" / "diagnoses_icd.csv.gz", index=False)

# --- labevents: morning peak in real clock time ----------------------------
lab, spec = [], []
sid = 0
for r in adm.sample(frac=0.4, random_state=1).itertuples():
    for _ in range(5):
        d = int(rng.integers(0, 5))
        u = rng.random()
        m = 0 if u < 0.03 else int(rng.integers(4 * 60, 7 * 60)) if u < 0.5 \
            else int(rng.integers(0, 1440))
        t = (r.real_admit.normalize() + pd.Timedelta(days=d, minutes=m)
             + r.shift)
        sid += 1
        for _ in range(int(rng.integers(1, 6))):  # results per specimen
            lab.append(t.strftime("%Y-%m-%d %H:%M:%S"))
            spec.append(sid)
lab_df = pd.DataFrame({"labevent_id": range(len(lab)), "specimen_id": spec,
                       "charttime": lab})
lab_df.loc[lab_df.sample(50, random_state=2).index, "charttime"] = np.nan
lab_df.to_csv(mroot / "hosp" / "labevents.csv.gz", index=False)

# --- patients --------------------------------------------------------------
first = adm.groupby("subject_id")["real_admit"].min()
grp = pd.cut(first.dt.year, [2007, 2010, 2013, 2016, 2019],
             labels=["2008 - 2010", "2011 - 2013", "2014 - 2016",
                     "2017 - 2019"])
pd.DataFrame({"subject_id": first.index, "anchor_year_group": grp.values}) \
    .to_csv(mroot / "hosp" / "patients.csv.gz", index=False)

# --- icustays ----------------------------------------------------------------
icu = adm.sample(frac=0.3, random_state=3).copy()
icu["intime"] = icu["admittime"] + pd.to_timedelta(
    rng.integers(1, 24, size=len(icu)), unit="h")
icu_days = rng.integers(1, 7, size=len(icu))
icu["outtime"] = (icu["intime"] - icu["shift"]).dt.normalize() \
    + pd.to_timedelta(icu_days, unit="D") \
    + pd.to_timedelta(midday_minute(len(icu)) * 60
                      + rng.integers(0, 60, len(icu)), unit="s") \
    + icu["shift"]
icu["stay_id"] = range(30000000, 30000000 + len(icu))

# --- chartevents: hourly HR/RR on the local clock (date shift = whole days)
def circ(hour, mean, amp, peak):
    return mean + amp * np.cos(2 * np.pi * (hour - peak) / 24)

ce_parts = []
for r in icu.itertuples():
    t = pd.date_range(r.intime.ceil("h"), r.outtime, freq="h")
    if len(t) == 0:
        continue
    hrs = t.hour.to_numpy()
    eff = rng.normal(0, 10)
    hr = circ(hrs, 85 + eff, 8, 16) + rng.normal(0, 4, len(t))
    rr = circ(hrs, 18, 2, 15) + rng.normal(0, 1.5, len(t))
    ts = t.strftime("%Y-%m-%d %H:%M:%S")
    ce_parts.append(pd.DataFrame({"stay_id": r.stay_id, "itemid": 220045,
                                  "charttime": ts, "valuenum": hr}))
    ce_parts.append(pd.DataFrame({"stay_id": r.stay_id, "itemid": 220210,
                                  "charttime": ts, "valuenum": rr}))
    ce_parts.append(pd.DataFrame({"stay_id": r.stay_id, "itemid": 220046,
                                  "charttime": ts[:3], "valuenum": 120.0}))
pd.concat(ce_parts).to_csv(mroot / "icu" / "chartevents.csv.gz", index=False)
pd.DataFrame({"itemid": [220045, 220046, 220210],
              "label": ["Heart Rate", "Heart rate Alarm - High",
                        "Respiratory Rate"]}) \
    .to_csv(mroot / "icu" / "d_items.csv.gz", index=False)

for c in ("intime", "outtime"):
    icu[c] = icu[c].dt.strftime("%Y-%m-%d %H:%M:%S")
icu[["subject_id", "hadm_id", "stay_id", "intime", "outtime"]] \
    .to_csv(mroot / "icu" / "icustays.csv.gz", index=False)

# --- eICU ----------------------------------------------------------------------
eroot = out / "eicu-crd" / "2.0"
eroot.mkdir(parents=True, exist_ok=True)
tzmode = sys.argv[3] if len(sys.argv) > 3 else "local"
n = 16000
hosp_ids = np.arange(1, 49)
regions = np.array(["Northeast", "Midwest", "South", "West"])
hosp_region = regions[hosp_ids % 4]
utc_off = {"Northeast": 300, "Midwest": 360, "South": 330, "West": 480}
hid = rng.choice(hosp_ids, n)
reg_of_stay = hosp_region[hid - 1]
shift = np.array([utc_off[r] if tzmode == "utc" else 0 for r in reg_of_stay])
admit_local = np.where(rng.random(n) < 0.4,
                       rng.integers(10 * 60, 16 * 60, n),
                       rng.integers(0, 1440, n))
admit_min = (admit_local + shift) % 1440  # recorded clock
unit_disch_local = midday_minute(n, 0.4)
los_min = (rng.integers(0, 7, n) * 1440 + unit_disch_local - admit_local)
los_min = np.where(los_min < 60, los_min + 1440, los_min)
hosp_disch_local = midday_minute(n)
morning_peak = rng.random(n) < 0.25            # E8: local heart-rate peak 09:00
hr_peak = np.where(morning_peak, 9.0, 16.0)
hstatus = np.where(rng.random(n) < np.where(morning_peak, 0.15, 0.04),
                   "Expired", "Alive")
pid = np.arange(n)
second = rng.random(n) < 0.1                   # some patients have two stays
pid[1:][second[1:]] = pid[:-1][second[1:]]
ustatus = np.where((hstatus == "Expired") & (rng.random(n) < 0.6),
                   "Expired", "Alive")
pd.DataFrame({"hospitalid": hosp_ids, "numbedscategory": 1,
              "teachingstatus": 0,
              "region": np.where(hosp_ids % 12 == 0, None, hosp_region)}) \
    .to_csv(eroot / "hospital.csv.gz", index=False)
pat = pd.DataFrame({
    "patientunitstayid": range(1, n + 1),
    "hospitalid": hid,
    "unitadmittime24": [f"{m // 60:02d}:{m % 60:02d}:00" for m in admit_min],
    "unitdischargeoffset": los_min,
    "hospitaldischargeyear": rng.choice([2014, 2015], n),
    "hospitaladmittime24": "08:00:00",
    "unitdischargetime24": [f"{(m + l) % 1440 // 60:02d}:{(m + l) % 60:02d}:00"
                            for m, l in zip(admit_min, los_min)],
    "patienthealthsystemstayid": range(500001, 500001 + n),
    "hospitaldischargestatus": np.where(rng.random(n) < 0.01, None, hstatus),
    "uniquepid": [f"{100000 + i:06d}-x" for i in pid],
    "unitvisitnumber": 1,
    "apacheadmissiondx": np.where(rng.random(n) < 0.2,
                                  "CVA, cerebrovascular accident/stroke",
                                  "Sepsis, pulmonary"),
    "unitdischargestatus": ustatus,
    "hospitaldischargetime24": [
        f"{(m + sh) % 1440 // 60:02d}:{(m + sh) % 60:02d}:00"
        for m, sh in zip(hosp_disch_local, shift)],
})
pat.to_csv(eroot / "patient.csv.gz", index=False)
died = (hstatus == "Expired").astype(float)
score = np.round(50 + 25 * died + rng.normal(0, 15, n)).clip(0, None)
score[rng.random(n) < 0.03] = -1
pd.concat([
    pd.DataFrame({"patientunitstayid": range(1, n + 1), "apacheversion": "IV",
                  "apachescore": score}),
    pd.DataFrame({"patientunitstayid": range(1, n + 1), "apacheversion": "IVa",
                  "apachescore": score + 3}),
]).to_csv(eroot / "apachePatientResult.csv.gz", index=False)
le = []
for i in range(n):
    for _ in range(10):
        d = int(rng.integers(0, max(1, los_min[i] // 1440 + 1)))
        clock = int(rng.integers(5 * 60, 7 * 60)) if rng.random() < 0.5 \
            else int(rng.integers(0, 1440))  # local clock
        off = d * 1440 + clock - admit_local[i]
        if 0 <= off <= los_min[i]:
            le.append((i + 1, off))
pd.DataFrame(le, columns=["patientunitstayid", "labresultoffset"]) \
    .to_csv(eroot / "lab.csv.gz", index=False)
pd.DataFrame(columns=["patientunitstayid", "day1meds", "ventday1"]) \
    .to_csv(eroot / "apachePredVar.csv.gz", index=False)
n_obs = los_min // 30 + 1
stay_rep = np.repeat(np.arange(n), n_obs)
offs = np.concatenate([np.arange(k) * 30 for k in n_obs])
local_h = ((admit_local[stay_rep] + offs) % 1440) / 60
eff = rng.normal(0, 10, n)[stay_rep]
pd.DataFrame({
    "patientunitstayid": stay_rep + 1,
    "observationoffset": offs,
    "heartrate": np.round(circ(local_h, 85 + eff, 8, hr_peak[stay_rep])
                          + rng.normal(0, 4, len(offs)), 1),
    "respiration": np.round(circ(local_h, 18, 2, 15)
                            + rng.normal(0, 1.5, len(offs)), 1),
}).to_csv(eroot / "vitalPeriodic.csv.gz", index=False)
print(f"{variant}: admissions {len(adm):,}, labevents {len(lab):,}, "
      f"icustays {len(icu):,}, eICU stays {n:,}, eICU lab {len(le):,}")
