import sys, json, time, urllib.request
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from sitelist import read_sites
sites = read_sites(sys.argv[1]); out = Path(sys.argv[2])
for r in sites.itertuples():
    f = out / f"{r.site_id}.csv"
    if f.exists(): continue
    url = (f"https://archive-api.open-meteo.com/v1/archive?latitude={r.lat}&longitude={r.lon}&start_date=2014-12-01"
           f"&end_date=2024-12-31&hourly=relative_humidity_2m&timezone=GMT")
    for k in range(6):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "research-script/0.1"}), timeout=300) as resp:
                h = json.loads(resp.read())["hourly"]; break
        except Exception as e:
            print("retry", r.site_id, e, flush=True); time.sleep(10*(k+1))
    pd.DataFrame({"hour": h["time"], "RH": h["relative_humidity_2m"]}).to_csv(f, index=False)
    print(r.site_id, len(h["time"]), flush=True); time.sleep(1)
