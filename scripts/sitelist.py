"""Read a site list (discovery: site_id/regime; confirmation: site/class) into site_id, cls, lat, lon."""
import pandas as pd
CLASS_MAP = {"urban": "pollution", "smoke": "burning"}   # discovery labels -> preregistered class names
def read_sites(path):
    d = pd.read_csv(path).rename(columns={"site": "site_id", "class": "cls", "regime": "cls"})
    d["cls"] = d["cls"].replace(CLASS_MAP)
    return d[["site_id", "cls", "lat", "lon"]]
