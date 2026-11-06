import sys, io, time, urllib.request
from pathlib import Path
import pandas as pd
out = Path(sys.argv[1])
def get(url):
    for k in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "research-script/0.1"}), timeout=600) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            print("retry", url, e, flush=True); time.sleep(20*(k+1))
for prod, cols in [("AOD20", {"AOD_500nm": "AOD500", "AOD_440nm": "AOD440", "440-870_Angstrom_Exponent": "AE"}),
                   ("SDA20", {"FineModeFraction_500nm[eta]": "FMF"})]:
    for y in range(2015, 2025):
        f = out / f"{prod}_{y}.csv"
        if f.exists(): continue
        txt = get(f"https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?year={y}&month=1&day=1&year2={y}&month2=12&day2=31&{prod}=1&AVG=20&if_no_html=1")
        lines = txt.splitlines(); h = next(i for i, l in enumerate(lines) if l.startswith("AERONET_Site"))
        df = pd.read_csv(io.StringIO("\n".join(lines[h:])), na_values=[-999, "-999.", "-999.000000"], low_memory=False)
        df.columns = [c.strip() for c in df.columns]
        dcol = [c for c in df.columns if c.startswith("Date")][0]
        keep = df[["AERONET_Site", dcol] + list(cols)].rename(columns={**cols, "AERONET_Site": "site", dcol: "date"})
        keep.to_csv(f, index=False); print(prod, y, len(keep), keep.site.nunique(), flush=True)
print("DONE")
