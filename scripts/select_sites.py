"""Objective confirmation-site selection. Uses only site metadata and daily L2 climatology (no hourly RH response)."""
import sys, glob, json, hashlib
from pathlib import Path
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
meta, discovery_csv, out_csv = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
SEED, CAP = 20261105, 25
R = 6371.0
def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.column_stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])
def chord_to_km(c): return 2*R*np.arcsin(np.clip(c/2, 0, 1))
loc = pd.read_csv(meta/"locations.txt", skiprows=1)
loc.columns = ["site", "lon", "lat", "elev"]
# daily climatology
aod = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(meta/"AOD20_*.csv")))])
sda = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(meta/"SDA20_*.csv")))])
aod["AOD"] = aod.AOD500.fillna(aod.AOD440)
aod["month"] = pd.to_datetime(aod.date, format="%d:%m:%Y").dt.month
a = aod.dropna(subset=["AOD", "AE"]).groupby("site").agg(n_aod=("AOD", "size"), medAOD=("AOD", "median"), medAE=("AE", "median"))
ga = aod.dropna(subset=["AOD"])
a["bb_ratio"] = ga[ga.month.isin([8, 9, 10])].groupby("site").AOD.median() / ga[ga.month.isin([12, 1, 2, 3, 4, 5])].groupby("site").AOD.median()
s = sda.dropna(subset=["FMF"]).groupby("site").agg(n_sda=("FMF", "size"), medFMF=("FMF", "median"))
c = loc.set_index("site").join(a, how="inner").join(s, how="inner").reset_index()
# coast distance
coast = json.load(open(meta/"ne_10m_coastline.geojson", encoding="utf-8"))
pts = []
for ft in coast["features"]:
    g = ft["geometry"]; lines = g["coordinates"] if g["type"] == "MultiLineString" else [g["coordinates"]]
    for ln in lines: pts.extend(ln)
pts = np.array(pts)
c["dist_coast_km"] = chord_to_km(cKDTree(xyz(pts[:, 1], pts[:, 0])).query(xyz(c.lat, c.lon))[0])
# big cities
pp = json.load(open(meta/"ne_10m_populated_places_simple.geojson", encoding="utf-8"))
big = np.array([f["geometry"]["coordinates"] for f in pp["features"] if (f["properties"].get("pop_max") or 0) >= 1_000_000])
c["dist_city1m_km"] = chord_to_km(cKDTree(xyz(big[:, 1], big[:, 0])).query(xyz(c.lat, c.lon))[0])
# discovery exclusion
disc = pd.read_csv(discovery_csv)
dd = chord_to_km(cKDTree(xyz(disc.lat, disc.lon)).query(xyz(c.lat, c.lon))[0])
c["dist_discovery_km"] = dd
box = lambda la0, la1, lo0, lo1: c.lat.between(la0, la1) & c.lon.between(lo0, lo1)
bb = box(-25, 5, -75, -40) | box(-25, 0, 10, 40)
seasia = box(10, 25, 95, 110)
elig = (~c.site.isin(disc.site_id)) & (dd > 100) & (c.elev < 1500) & (c.n_aod >= 300) & (c.n_sda >= 300)
cls = pd.Series("none", index=c.index)
dust = (c.medAE < 0.75) & (c.medAOD >= 0.15)
marine = ~dust & (c.dist_coast_km <= 5) & (c.dist_city1m_km > 50) & (c.medAOD < 0.15) & (c.medFMF < 0.6) & (c.medAE < 1.0)
burn = ~dust & ~marine & bb & (c.bb_ratio >= 2.0) & (c.medFMF >= 0.6)
poll = ~dust & ~marine & ~burn & ~seasia & (c.medFMF >= 0.7) & (c.medAE >= 1.2) & (c.dist_city1m_km <= 50)
for name, m in [("dust", dust), ("marine", marine), ("burning", burn), ("pollution", poll)]:
    cls[m & (cls == "none")] = name
c["class"] = cls; c["eligible"] = elig
pool = c[elig & (c["class"] != "none")].sort_values("n_aod", ascending=False)
# de-duplicate: within a class, drop a site within 25 km of a kept site with more data
keep = []
for k, g in pool.groupby("class"):
    kept = []
    for _, r in g.iterrows():
        if kept and chord_to_km(cKDTree(xyz(np.array([x.lat for x in kept]), np.array([x.lon for x in kept]))).query(xyz([r.lat], [r.lon]))[0])[0] <= 50:
            continue
        kept.append(r)
    kg = pd.DataFrame(kept)
    n_pool = len(kg)
    if len(kg) > CAP: kg = kg.sample(CAP, random_state=SEED)
    kg["n_class_pool"] = n_pool; keep.append(kg)
sel = pd.concat(keep).sort_values(["class", "site"])
cols = ["class", "site", "lat", "lon", "elev", "n_aod", "n_sda", "medAOD", "medAE", "medFMF", "bb_ratio", "dist_coast_km", "dist_city1m_km", "dist_discovery_km", "n_class_pool"]
sel[cols].round(3).to_csv(out_csv, index=False)
print("candidates with both products:", len(c), "| eligible:", int(elig.sum()))
print("eligible by class (before dedup):", c[elig].groupby("class").size().to_dict())
print(sel.groupby("class").agg(n=("site", "size"), pool=("n_class_pool", "first")).to_string())
print(sel[cols].round(2).to_string(index=False))
print("SHA-256:", hashlib.sha256(open(out_csv, "rb").read()).hexdigest())
