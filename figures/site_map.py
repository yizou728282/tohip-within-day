"""Figure 1: map of the 72 confirmation sites and 18 discovery sites, coloured by aerosol class.

usage: python site_map.py CONFIRMATION_SITES DISCOVERY_SITES NE_COASTLINE_GEOJSON OUT_STEM
"""
import sys, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pyproj import Transformer

CONF, DISC, COAST, OUT = sys.argv[1:5]
CLASSES = [("pollution", "Pollution", "#2a78d6"), ("dust", "Dust", "#eb6834"),
           ("marine", "Marine", "#1baf7a"), ("burning", "Biomass burning", "#eda100")]
INK, INK2, GRID, AXIS = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"

conf = pd.read_csv(CONF).rename(columns={"class": "cls"})
disc = pd.read_csv(DISC).rename(columns={"site_id": "site", "regime": "cls", "class": "cls"})
disc["cls"] = disc.cls.replace({"urban": "pollution", "smoke": "burning"})
tr = Transformer.from_crs("EPSG:4326", "ESRI:54030", always_xy=True)      # Robinson

fig, ax = plt.subplots(figsize=(7.09, 4.0), dpi=300)                      # 180 mm wide
plt.rcParams.update({"font.size": 7})

# graticule and outline
for lat in range(-60, 91, 30):
    x, y = tr.transform(np.linspace(-180, 180, 361), np.full(361, lat)); ax.plot(x, y, color=GRID, lw=0.4, zorder=1)
for lon in range(-180, 181, 60):
    x, y = tr.transform(np.full(181, lon), np.linspace(-90, 90, 181))
    ax.plot(x, y, color=AXIS if abs(lon) == 180 else GRID, lw=0.6 if abs(lon) == 180 else 0.4, zorder=1)
x, y = tr.transform(np.linspace(-180, 180, 361), np.full(361, -90.0)); ax.plot(x, y, color=AXIS, lw=0.6, zorder=1)
x, y = tr.transform(np.linspace(-180, 180, 361), np.full(361, 90.0)); ax.plot(x, y, color=AXIS, lw=0.6, zorder=1)

# coastlines
for ft in json.load(open(COAST, encoding="utf-8"))["features"]:
    g = ft["geometry"]; lines = g["coordinates"] if g["type"] == "MultiLineString" else [g["coordinates"]]
    for ln in lines:
        a = np.asarray(ln)
        if len(a) < 2:
            continue
        x, y = tr.transform(a[:, 0], a[:, 1]); ax.plot(x, y, color="#9a998f", lw=0.25, zorder=2)

handles = []
for key, label, col in CLASSES:
    c, d = conf[conf.cls == key], disc[disc.cls == key]
    x, y = tr.transform(c.lon.values, c.lat.values)
    ax.scatter(x, y, s=16, c=col, edgecolors="white", linewidths=0.5, zorder=4)
    x, y = tr.transform(d.lon.values, d.lat.values)
    ax.scatter(x, y, s=26, marker="^", facecolors="white", edgecolors=col, linewidths=1.1, zorder=5)
    handles.append(Line2D([], [], ls="", marker="o", ms=5, mfc=col, mec="white", mew=0.5,
                          label=f"{label}: {len(c)} confirmation, {len(d)} discovery"))
handles += [Line2D([], [], ls="", marker="o", ms=5, mfc=INK2, mec="white", label="Confirmation site (filled circle)"),
            Line2D([], [], ls="", marker="^", ms=5.5, mfc="white", mec=INK2, mew=1.1, label="Discovery site (open triangle)")]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=3, frameon=False,
          fontsize=6.5, handletextpad=0.3, columnspacing=1.2, labelcolor=INK)
ax.set_aspect("equal"); ax.axis("off")
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.12)
fig.savefig(OUT + ".png", dpi=300, facecolor="white")
fig.savefig(OUT + ".pdf", facecolor="white")
print("confirmation:", conf.cls.value_counts().to_dict(), "| discovery:", disc.cls.value_counts().to_dict())
