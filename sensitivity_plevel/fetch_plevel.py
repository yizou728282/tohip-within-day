"""DEVIATION / SUBSTITUTE for preregistered sensitivity #6 (ERA5 925/850 hPa RH via CDS).

CDS was not accessible (no account or API key on this machine), and the Open-Meteo ERA5 archive does not serve
pressure levels. This script instead downloads hourly GFS (Open-Meteo historical-forecast API, model gfs_global)
2 m, 925 hPa and 850 hPa relative humidity plus the geopotential height of both levels, for 2021-03-01 to
2024-12-31 (the period for which these fields exist in that archive).

A failed site raises and writes nothing (unlike the frozen fetch_rh.py, which reused the previous site's data).

usage: python fetch_plevel.py SITE_LIST OUT_DIR [CLASS ...]
"""
import sys, json, time, urllib.request, urllib.error
from pathlib import Path
import pandas as pd

URL = ("https://historical-forecast-api.open-meteo.com/v1/forecast?latitude={lat}&longitude={lon}"
       "&start_date=2021-03-01&end_date=2024-12-31&models=gfs_global&timezone=GMT"
       "&hourly=relative_humidity_2m,relative_humidity_925hPa,relative_humidity_850hPa,"
       "geopotential_height_925hPa,geopotential_height_850hPa")
UA = {"User-Agent": "research-script/0.1"}


def get(url, tries=8):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            wait = 90 * (k + 1) if e.code == 429 else 20 * (k + 1)
            print(f"  HTTP {e.code}; retry in {wait}s", flush=True); time.sleep(wait)
        except Exception as e:
            print(f"  {e!r}; retry in 30s", flush=True); time.sleep(30)
    raise RuntimeError("download failed after retries")


def main():
    sites, out = pd.read_csv(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    classes = sys.argv[3:] or sorted(sites["class"].unique())
    todo = sites[sites["class"].isin(classes)]
    for s in todo.itertuples():
        f = out / f"{s.site}.csv"
        if f.exists():
            continue
        js = get(URL.format(lat=s.lat, lon=s.lon)); h = js["hourly"]
        d = pd.DataFrame({"hour": pd.to_datetime(h["time"]),
                          "RH2m_gfs": h["relative_humidity_2m"], "RH925": h["relative_humidity_925hPa"],
                          "RH850": h["relative_humidity_850hPa"], "z925": h["geopotential_height_925hPa"],
                          "z850": h["geopotential_height_850hPa"]})
        d = d.dropna(subset=["RH925", "RH850"], how="all")
        d["grid_lat"], d["grid_lon"], d["grid_elev"] = js["latitude"], js["longitude"], js.get("elevation")
        d.to_csv(f, index=False)
        print("done", s.site, len(d), d.hour.min(), d.hour.max(), flush=True)
        time.sleep(12)


if __name__ == "__main__":
    main()
