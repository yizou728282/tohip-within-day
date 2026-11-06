import sys, io, time, json, urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from sitelist import read_sites
sites = read_sites(sys.argv[1]); out = Path(sys.argv[2]); out.mkdir(exist_ok=True)
Y0, Y1 = 2015, 2024
UA = {"User-Agent": "research-script/0.1"}
def get(url, tries=4, timeout=180):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            err = e; time.sleep(5*(k+1))
    print("FAIL", url, err, flush=True); return None
def aeronet(site, y, prod):
    url = (f"https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?site={site}&year={y}&month=1&day=1"
           f"&year2={y}&month2=12&day2=31&{prod}=1&AVG=10&if_no_html=1")
    txt = get(url)
    if not txt or "AERONET_Site" not in txt: return None
    lines = txt.splitlines(); h = next(i for i,l in enumerate(lines) if l.startswith("AERONET_Site"))
    if len(lines) <= h+1 or not lines[h+1].startswith(site + ","): return None
    df = pd.read_csv(io.StringIO("\n".join(lines[h:])), na_values=[-999, -999.0, "-999.000000"], low_memory=False)
    df.columns = [c.strip() for c in df.columns]
    dcol = [c for c in df.columns if c.startswith("Date")][0]; tcol = [c for c in df.columns if c.startswith("Time")][0]
    df["t"] = pd.to_datetime(df[dcol] + " " + df[tcol], format="%d:%m:%Y %H:%M:%S")
    if prod == "AOD20":
        keep = {"AOD_500nm":"AOD500","AOD_440nm":"AOD440","440-870_Angstrom_Exponent":"AE",
                "Precipitable_Water(cm)":"PW","Solar_Zenith_Angle(Degrees)":"SZA"}
    else:
        keep = {"Total_AOD_500nm[tau_a]":"tau_a","Fine_Mode_AOD_500nm[tau_f]":"tau_f","Coarse_Mode_AOD_500nm[tau_c]":"tau_c",
                "FineModeFraction_500nm[eta]":"FMF","AE-Fine_Mode_500nm[alpha_f]":"alpha_f"}
    d = df[["t"]+list(keep)].rename(columns=keep)
    d["hour"] = d.t.dt.floor("h")
    g = d.drop(columns="t").groupby("hour")
    a = g.mean(); a[f"n_{prod}"] = g.size()
    return a
def openmeteo(lat, lon):
    url = (f"https://archive-api.open-meteo.com/v1/archive?latitude={lat}&longitude={lon}&start_date={Y0}-01-01"
           f"&end_date={Y1}-12-31&hourly=relative_humidity_2m,temperature_2m,dew_point_2m&timezone=GMT")
    js = json.loads(get(url, tries=6, timeout=300)); h = js["hourly"]
    m = pd.DataFrame({"hour": pd.to_datetime(h["time"]), "RH": h["relative_humidity_2m"],
                      "T2m": h["temperature_2m"], "Td2m": h["dew_point_2m"]}).set_index("hour")
    m.attrs["om_elev"] = js.get("elevation"); return m
def do_site(row):
    site = row.site_id; f = out / f"{site}.csv"
    if f.exists(): return site, "cached"
    parts_a, parts_s = [], []
    for y in range(Y0, Y1+1):
        a = aeronet(site, y, "AOD20"); s = aeronet(site, y, "SDA20")
        if a is not None: parts_a.append(a)
        if s is not None: parts_s.append(s)
        time.sleep(1)
    if not parts_a: return site, "no AOD"
    A = pd.concat(parts_a)
    if parts_s: A = A.join(pd.concat(parts_s), how="left")
    m = openmeteo(row.lat, row.lon)
    A = A.join(m, how="inner"); A["site_id"] = site; A["regime"] = row.cls; A["om_elev"] = m.attrs["om_elev"]
    A.to_csv(f); return site, f"{len(A)} hours"
with ThreadPoolExecutor(3) as ex:
    for site, msg in ex.map(do_site, sites.itertuples()):
        print(site, msg, flush=True)
print("DONE")
