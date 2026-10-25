"""Final maps and charts of the Syntrichia manuscript, with one common style.

Produces (figures/):
  Fig1_occurrences_genus_top5_species   Fig. 1  genus (A) + the five modelled species (B-F)
  Fig2_records_and_richness             Fig. 2  records (A) and species richness (B) per 50 km hexagon
  Fig4_protected_areas                  Fig. 4  records in protected areas (A-E)
  FigS1_records_per_species             Fig. S1 records per species
  summary.json, hex_grid_50km_summary.csv, TableS3_species_in_protected_areas.csv

Each figure is saved as PNG (600 dpi), PDF (vector) and TIFF (RGB, 600 dpi, LZW).

Common style
- Projection: Albers equal-area conic centred on Chile (no appreciable shape distortion along
  the country, equal-area cells).
- Base map: Natural Earth 1:10m. Pale blue-grey sea, warm grey neighbouring countries,
  white Chile with a grey outline. Graticule every 5 degrees, 500 km scale bar.
- Typeface: Liberation Sans (Arial metrics), 7-8 pt.
- Colour: blue = records / main magnitude; orange = richness / outside protected areas.
  One-hue ramps (monotonic, minimum contrast 2:1).
- Hexagonal grid: regular hexagons, 50 km between opposite sides (~2,165 km2),
  built in the equal-area projection.

Inputs
- data/figshare/Syntrichia_occurrences_with_cleaning_flags.csv (06_figshare_dataset.R;
  the records of the maps are those with usedChileMaps == TRUE)
- protected_area_layers/*.gpkg (11_download_protected_area_layers.R)
- Natural Earth 1:10m countries (downloaded on first run into data/natural_earth/)

Usage:  python python/make_figures.py
"""
import math, json, os
import numpy as np, pandas as pd, geopandas as gpd
from shapely.geometry import Polygon, box, LineString, Point
from shapely import make_valid
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D

# ------------------------------------------------------------------ paths
DATA = "data/figshare/Syntrichia_occurrences_with_cleaning_flags.csv"
NE_URL = "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip"
NE_DIR = "data/natural_earth"
NE = NE_DIR + "/ne_10m_admin_0_countries.shp"
ANP = "protected_area_layers/"
OUT = "figures"
os.makedirs(OUT, exist_ok=True); os.makedirs(NE_DIR, exist_ok=True)

# ------------------------------------------------------------------ style
AEA = "+proj=aea +lat_1=-22 +lat_2=-50 +lat_0=-36 +lon_0=-71 +datum=WGS84 +units=m +no_defs"
FONT = "Liberation Sans"
S = 1.12          # common scaling factor for texts (not legends)
plt.rcParams.update({"font.family": FONT, "font.size": 7.5, "axes.titlesize": 8,
                     "legend.fontsize": 9, "legend.title_fontsize": 9.5,
                     "pdf.fonttype": 42, "svg.fonttype": "none",
                     "mathtext.fontset": "custom", "mathtext.rm": FONT, "mathtext.bf": FONT + ":bold", "mathtext.it": FONT + ":italic"})
C = dict(sea="#bcd4e6", land="#efe9df", chile="#ffffff", coast="#6a6f73", border="#c9cdd1",
         grid="#aec6d9", ink="#0b0b0b", ink2="#4a4f55", muted="#7d8288", empty="#d5dade",
         blue="#2a78d6", aqua="#1baf7a", orange="#eb6834", context="#b5bbc1")
BLUES = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]      # ordinal, validated
ORANGES = ["#ef9b70", "#e57a48", "#d95926", "#ad4317", "#7a2d0c"]    # ordinal, validated
HEX_W = 50_000          # m between opposite sides
LON0, LON1, LAT0, LAT1 = -77.0, -64.0, -56.3, -16.9    # extent of the single maps

# ------------------------------------------------------------------ data
ds = pd.read_csv(DATA, low_memory=False)
occ = ds[ds.usedChileMaps.astype(bool)].copy()
assert len(occ) == 1857
pts = gpd.GeoDataFrame(occ, geometry=gpd.points_from_xy(occ.decimalLongitude, occ.decimalLatitude), crs=4326).to_crs(AEA)
# S. prostrata = S. serrulata (Gallego et al. 2009, Muller 2009, GBIF). Already applied in
# analysisSpecies by 06_figshare_dataset.R; kept here as a safeguard
SYN = {"Syntrichia prostrata": "Syntrichia serrulata"}
pts["analysisSpecies"] = pts.analysisSpecies.replace(SYN)
sp = pts[pts.analysisSpecies.notna()].copy()

if not os.path.exists(NE):
    import urllib.request, zipfile, io
    print("Downloading Natural Earth 1:10m countries ...")
    zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(NE_URL).read())).extractall(NE_DIR)
world = gpd.read_file(NE)
frame_ll = box(-82, -60, -55, -14)
chile = world[world.ADMIN == "Chile"].clip(box(-76.5, -60, -60, -14)).to_crs(AEA)   # without oceanic islands
neigh = world[world.ADMIN.isin(["Argentina", "Peru", "Bolivia", "Brazil", "Paraguay", "Uruguay"])].clip(frame_ll).to_crs(AEA)

def ext(lon0=LON0, lon1=LON1, lat0=LAT0, lat1=LAT1):
    g = gpd.GeoSeries([box(lon0, lat0, lon1, lat1)], crs=4326).to_crs(AEA)
    # projected box that contains the geographic box
    return g.total_bounds

def graticule(ax, b, lab=True, lab_lon=True):
    for lat in range(-55, -15, 5):
        ln = gpd.GeoSeries([LineString([(lon, lat) for lon in np.linspace(-85, -50, 80)])], crs=4326).to_crs(AEA)
        ln.plot(ax=ax, color=C["grid"], linewidth=0.35, zorder=1.5)
        if lab:
            x, y = gpd.GeoSeries([Point(-76.2, lat)], crs=4326).to_crs(AEA).iloc[0].coords[0]
            if b[1] < y < b[3]:
                ax.text(b[0] + (b[2]-b[0])*0.015, y, f"{-lat}°S", fontsize=6*S, color=C["muted"], va="bottom", ha="left", zorder=6)
    for lon in (-75, -70, -65):
        ln = gpd.GeoSeries([LineString([(lon, lat) for lat in np.linspace(-60, -14, 80)])], crs=4326).to_crs(AEA)
        ln.plot(ax=ax, color=C["grid"], linewidth=0.35, zorder=1.5)
        if lab_lon == "out":                        # label below the frame (narrow panels)
            xs, ys = np.array(ln.iloc[0].coords).T
            xb = np.interp(b[1], ys, xs)
            if b[0] < xb < b[2]:
                ax.annotate(f"{-lon}°W", xy=(xb, b[1]), xytext=(0, -2), textcoords="offset points",
                            fontsize=6*S, color=C["muted"], ha="center", va="top", annotation_clip=False)
        elif lab and lab_lon:
            x, y = gpd.GeoSeries([Point(lon, LAT0 + 0.25)], crs=4326).to_crs(AEA).iloc[0].coords[0]
            if b[0] < x < b[2]:
                ax.text(x, b[1] + (b[3]-b[1])*0.006, f"{-lon}°W", fontsize=6*S, color=C["muted"], ha="center", va="bottom", zorder=6)

def scalebar(ax, b, km=500, frac_x=0.55, frac_y=0.17, mid=True):
    x0 = b[0] + (b[2]-b[0])*frac_x; y0 = b[1] + (b[3]-b[1])*frac_y
    L = km*1000
    ax.add_patch(Rectangle((x0, y0), L/2, (b[3]-b[1])*0.004, fc=C["ink2"], ec=C["ink2"], lw=0.4, zorder=7))
    ax.add_patch(Rectangle((x0+L/2, y0), L/2, (b[3]-b[1])*0.004, fc="white", ec=C["ink2"], lw=0.4, zorder=7))
    for xx, t in [(x0, "0")] + ([(x0+L/2, f"{km//2}")] if mid else []) + [(x0+L, f"{km} km")]:
        ax.text(xx, y0 + (b[3]-b[1])*0.008, t, fontsize=5.5*S, color=C["ink2"], ha="center", va="bottom", zorder=7)

def north(ax, b, fx=0.9, fy=0.955):
    x = b[0] + (b[2]-b[0])*fx; y = b[1] + (b[3]-b[1])*fy
    ax.annotate("N", xy=(x, y), xytext=(x, y - (b[3]-b[1])*0.03), ha="center", va="top", fontsize=7*S, color=C["ink2"],
                arrowprops=dict(arrowstyle="-|>", color=C["ink2"], lw=0.6), zorder=7)

def basemap(ax, b, lab=True, bar=True, lab_lon=True):
    ax.set_facecolor(C["sea"])
    neigh.plot(ax=ax, color=C["land"], edgecolor="#d8cfc0", linewidth=0.5, zorder=1)
    graticule(ax, b, lab, lab_lon)
    chile.plot(ax=ax, color=C["chile"], edgecolor=C["coast"], linewidth=0.35, zorder=2)
    ax.set_xlim(b[0], b[2]); ax.set_ylim(b[1], b[3]); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_color(C["border"]); s.set_linewidth(0.5)
    if bar: scalebar(ax, b); north(ax, b)

def single_map():
    b = ext()
    w = 3.5; h = w * (b[3]-b[1]) / (b[2]-b[0])
    fig = plt.figure(figsize=(w, h)); ax = fig.add_axes([0, 0, 1, 1])
    return fig, ax, b

def legend(ax, handles, title, loc=(0.985, 0.27)):
    lg = ax.legend(handles=handles, title=title, loc="center right", bbox_to_anchor=loc, frameon=True,
                   framealpha=0.95, edgecolor=C["border"], fancybox=False, borderpad=0.6, labelspacing=0.45,
                   handlelength=1.4, handleheight=1.1, alignment="left")
    lg.get_frame().set_linewidth(0.5); lg.set_zorder(8)
    return lg

def save(fig, name):
    for a in fig.axes:
        if a.get_xlabel() in ("Easting [metre]",): a.set_xlabel("")
        if a.get_ylabel() in ("Northing [metre]",): a.set_ylabel("")
    fig.savefig(f"{OUT}/{name}.png", dpi=600); fig.savefig(f"{OUT}/{name}.pdf")
    plt.close(fig)
    from PIL import Image                                  # TIFF RGB 600 dpi LZW for the journal
    im = Image.open(f"{OUT}/{name}.png"); bg = Image.new("RGB", im.size, "white")
    bg.paste(im, mask=im.split()[3] if im.mode == "RGBA" else None)
    bg.save(f"{OUT}/{name}.tiff", compression="tiff_lzw", dpi=(600, 600))

# ------------------------------------------------------------------ hexagonal grid (equal area)
def hex_grid(bounds, w):
    R = w / math.sqrt(3); dy = 1.5 * R
    xmin, ymin, xmax, ymax = bounds
    polys = []; row = 0; y = ymin - R
    while y < ymax + R:
        x = xmin - w + (w/2 if row % 2 else 0)
        while x < xmax + w:
            polys.append(Polygon([(x + R*math.sin(math.radians(60*k)), y + R*math.cos(math.radians(60*k))) for k in range(6)]))
            x += w
        y += dy; row += 1
    return gpd.GeoDataFrame(geometry=polys, crs=AEA)

grid = hex_grid(chile.total_bounds, HEX_W)
hit = gpd.sjoin(pts[["geometry"]], grid, predicate="intersects").index_right.unique()
grid = grid[grid.intersects(chile.union_all()) | grid.index.isin(hit)].reset_index(drop=True)
def per_cell(p):
    j = gpd.sjoin(p[["geometry", "analysisSpecies"]], grid, predicate="intersects")
    return j[~j.index.duplicated(keep="first")]
ja, js = per_cell(pts), per_cell(sp)
grid["n"] = grid.index.map(ja.groupby("index_right").size()).fillna(0).astype(int)
grid["rich"] = grid.index.map(js.groupby("index_right").analysisSpecies.nunique()).fillna(0).astype(int)
assert grid.n.sum() == len(pts)
summary = {"hexagons": len(grid), "with_records": int((grid.n > 0).sum()), "max_records": int(grid.n.max()),
           "max_richness": int(grid.rich.max()), "hex_area_km2": round(grid.geometry.area.iloc[0]/1e6, 1)}
print(summary)

def classed(ax, col, bins, labels, cols, title, below=False):
    grid[grid[col] == 0].plot(ax=ax, facecolor="none", edgecolor=C["empty"], linewidth=0.3, zorder=3)
    for (lo, hi), c in zip(bins, cols):
        s = grid[(grid[col] >= lo) & (grid[col] <= hi)]
        if len(s): s.plot(ax=ax, facecolor=c, edgecolor="white", linewidth=0.3, zorder=4)
    h = [Patch(facecolor=c, edgecolor="none", label=l) for l, c in zip(labels, cols)]
    h.append(Patch(facecolor="none", edgecolor=C["empty"], label="No records"))
    if not below:
        legend(ax, h, title); return
    lg = ax.legend(handles=h, title=title, loc="upper center", bbox_to_anchor=(0.5, -0.012), ncol=2, frameon=False,
                   handlelength=1.6, handleheight=1.25, columnspacing=1.3, alignment="center",
                   fontsize=11, title_fontsize=11.5)

# ---- Fig. 2: records (A) and species richness (B) per hexagon
def panel_title(ax, letter, text, fs_l=8*S, fs_t=7*S):
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(1, 5), textcoords="offset points",
                fontsize=fs_l, fontweight="bold", ha="left", va="baseline", color=C["ink"])
    ax.annotate(text, xy=(0, 1), xycoords="axes fraction", xytext=(1 + 1.05*fs_l, 5), textcoords="offset points",
                fontsize=fs_t, ha="left", va="baseline", color=C["ink"])
b2 = ext()
ratio = (b2[3]-b2[1]) / (b2[2]-b2[0])
tit2, leg2 = 0.34, 1.32
pw2 = min(3.55, (9.3 - tit2 - leg2) / ratio)
H2 = pw2*ratio + tit2 + leg2
fig, axs = plt.subplots(1, 2, figsize=(2*pw2 + 0.06, H2))
fig.subplots_adjust(left=0.004, right=0.996, top=1 - tit2/H2, bottom=leg2/H2, wspace=0.03)
nm, rm = summary["max_records"], summary["max_richness"]
basemap(axs[0], b2, lab=True, bar=True)
classed(axs[0], "n", [(1, 1), (2, 5), (6, 15), (16, 40), (41, nm)], ["1", "2–5", "6–15", "16–40", f"41–{nm}"], BLUES, "Records per hexagon", below=True)
panel_title(axs[0], "A", "Records", 12.5, 11.5)
basemap(axs[1], b2, lab=False, bar=False)
classed(axs[1], "rich", [(1, 2), (3, 4), (5, 7), (8, 10), (11, rm)], ["1–2", "3–4", "5–7", "8–10", f"11–{rm}"], ORANGES, "Species per hexagon", below=True)
panel_title(axs[1], "B", "Species richness", 12.5, 11.5)
save(fig, "Fig2_records_and_richness")
print("Fig2 size in:", round(2*pw2 + 0.06, 2), round(H2, 2))

# ---- Fig. 1: whole genus (A) and the five species with most records (B-F);
#      these are the five modelled species, in the same order as Fig. 3 (SDM)
top = sp.analysisSpecies.value_counts().head(5)
b1 = ext(-76.6, -64.6, LAT0, LAT1)
pw = 7.2 / 6; ph = pw * (b1[3]-b1[1]) / (b1[2]-b1[0]); tit1, leg1 = 0.30, 0.55
H1 = ph + tit1 + leg1
fig, axs = plt.subplots(1, 6, figsize=(7.2, H1))
fig.subplots_adjust(left=0.004, right=0.996, top=1 - tit1/H1, bottom=leg1/H1, wspace=0.035)
panels = [("A", "Syntrichia", None, len(pts))] + [(L, f"S. {s_.split()[1]}", s_, n_) for L, (s_, n_) in zip("BCDEF", top.items())]
for i, (ax, (L, name, s_, n_)) in enumerate(zip(axs, panels)):
    basemap(ax, b1, lab=(i == 0), bar=False, lab_lon="out")
    if s_ is None:
        pts.plot(ax=ax, color=C["blue"], markersize=7.2, edgecolor="white", linewidth=0.25, zorder=4)
    else:
        pts.plot(ax=ax, color=C["context"], markersize=1.6, linewidth=0, zorder=3)
        sp[sp.analysisSpecies == s_].plot(ax=ax, color=C["blue"], markersize=7.2, edgecolor="white", linewidth=0.25, zorder=4)
    ax.annotate(L, xy=(0, 1), xycoords="axes fraction", xytext=(1, 5), textcoords="offset points",
                fontsize=10.5, fontweight="bold", ha="left", va="baseline", color=C["ink"])
    ax.annotate(name, xy=(0, 1), xycoords="axes fraction", xytext=(1 + 1.05*10.5, 5), textcoords="offset points",
                fontsize=8.9, style="italic", ha="left", va="baseline", color=C["ink"])
    ax.text(0.975, 0.52, f"n = {n_:,}", transform=ax.transAxes, ha="right", va="top", fontsize=8.4, color=C["ink2"], zorder=8)
    ax.set_xlabel(""); ax.set_ylabel("")
y35 = gpd.GeoSeries([Point(-67, -35)], crs=4326).to_crs(AEA).iloc[0].y
y47 = gpd.GeoSeries([Point(-67, -47.5)], crs=4326).to_crs(AEA).iloc[0].y
scalebar(axs[0], b1, km=300, frac_x=0.6, frac_y=(y47 - b1[1]) / (b1[3] - b1[1]), mid=False)
fig.legend(handles=[Line2D([], [], marker="o", ls="", markersize=7, markerfacecolor=C["blue"], markeredgecolor="white", markeredgewidth=0.3, label="Records of the taxon in the panel"),
                    Line2D([], [], marker="o", ls="", markersize=4.5, markerfacecolor=C["context"], markeredgecolor="none", label="Other $\\it{Syntrichia}$ records")],
           loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, frameon=False, fontsize=10, columnspacing=2.0, handlelength=1.2)
save(fig, "Fig1_occurrences_genus_top5_species")
print("Fig1 size in:", 7.2, round(H1, 2))

# ---- Fig. S1: records per species
t = sp.analysisSpecies.value_counts().sort_values()
modelled = {"Syntrichia epilosa", "Syntrichia anderssonii", "Syntrichia pseudorobusta", "Syntrichia breviseta", "Syntrichia squarripila"}
fig, ax = plt.subplots(figsize=(4.2, 0.145*len(t) + 0.7))
cols = [C["blue"] if s in modelled else "#bcc1c6" for s in t.index]
ax.barh(range(len(t)), t.values, color=cols, height=0.72, zorder=3)
for i, v in enumerate(t.values): ax.text(v + 2, i, str(v), va="center", fontsize=6*S, color=C["ink2"])
ax.set_yticks(range(len(t))); ax.set_yticklabels([f"S. {s.split()[1]}" for s in t.index], fontsize=6.5*S, style="italic")
ax.set_xlabel("Records in Chile", fontsize=7*S, color=C["ink2"]); ax.tick_params(axis="x", labelsize=6.5*S, colors=C["ink2"], length=2)
ax.tick_params(axis="y", length=0)
ax.grid(axis="x", color="#e1e0d9", linewidth=0.5, zorder=0); ax.set_axisbelow(True)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color("#c3c2b7"); ax.set_ylim(-0.7, len(t)-0.3); ax.set_xlim(0, t.max()*1.1)
ax.legend(handles=[Patch(color=C["blue"], label="Modelled species"), Patch(color="#bcc1c6", label="Other species")],
          loc="lower right", frameon=False, fontsize=9)
fig.tight_layout(); save(fig, "FigS1_records_per_species")

# ---- Fig. 4: protected areas, 5 panels (A-D by type, E combined)
def layer(*names):
    g = pd.concat([gpd.read_file(ANP + n + ".gpkg")[["geometry"]] for n in names])
    g["geometry"] = g.geometry.apply(make_valid)
    g = g.explode(index_parts=False)                      # make_valid can leave stray points or lines
    g = g[g.geom_type.isin(["Polygon", "MultiPolygon"])]
    return g.to_crs(AEA)
def inside(g):
    return pts.index.isin(gpd.sjoin(pts[["geometry"]], g, predicate="intersects").index)
types = [("A", "State protected areas", layer("state_pa_mbn")),
         ("B", "Nature sanctuaries", layer("nature_sanctuaries")),
         ("C", "Private protected areas", layer("private_pa_a", "private_pa_b")),
         ("D", "Priority sites", layer("priority_sites_mma2020"))]
ins = {t: inside(g) for _, t, g in types}
in_off = ins["State protected areas"] | ins["Nature sanctuaries"] | ins["Private protected areas"]
in_pri = ins["Priority sites"]
official = pd.concat([g for _, t, g in types if t != "Priority sites"])
priority = types[3][2]
# Areas as solid patches (no outline, no transparency). A-D: one colour per type,
# different from the point colours (blue, orange, aqua) and from the base map.
TYPE_COL = {"State protected areas": "#7b62a8",    # soft purple
            "Nature sanctuaries": "#c9a227",       # mustard yellow
            "Private protected areas": "#c45a4e",  # soft brick red
            "Priority sites": "#8aac2e"}           # darker lime green
E_OFF, E_PRI = "#4e9a2e", "#c4ddb2"       # panel E: official areas (green) and priority sites (light green)
CHILE_LAND = chile.union_all()
def draw_pa(ax, g, col, z=2.5):
    # drawing only: clipped to the land outline of Chile (the layers include marine
    # parts and do not match the Natural Earth coastline). The figures use the full layers.
    gd = gpd.clip(g.assign(geometry=g.simplify(200)), CHILE_LAND)
    gd = gd[gd.geom_type.isin(["Polygon", "MultiPolygon"])]
    gd.plot(ax=ax, facecolor=col, edgecolor="none", linewidth=0, zorder=z)
    chile.boundary.plot(ax=ax, color=C["coast"], linewidth=0.35, zorder=2.8)   # coastline on top of the areas
def draw_pts(ax, mask, col, ms=7.2, z=4):
    p = pts[mask]
    if len(p): p.plot(ax=ax, color=col, markersize=ms, edgecolor="white", linewidth=0.25, zorder=z)
def pct(a, b):
    return "$\\bf{" + f"{100*a/b:.1f}" + "\\%}$"
def stats_text(ax, s_, n_, gap=4.5):
    """Species on top and records below, separated by a gap (in pt)."""
    kw = dict(xycoords="axes fraction", textcoords="offset points", ha="right", va="top",
              fontsize=5.8*S, color=C["ink2"], linespacing=1.1, zorder=8)
    ax.annotate(f"{s_}/{nsp} species\n{pct(s_, nsp)}", xy=(0.96, 0.515), xytext=(0, 0), **kw)
    two_lines = 2 * 5.8*S * 1.1 * 1.2          # approximate height of two lines (pt)
    ax.annotate(f"{n_} records\n{pct(n_, len(pts))}", xy=(0.96, 0.515), xytext=(0, -(two_lines + gap)), **kw)
def stats(m):
    return int(pts[m].analysisSpecies.dropna().nunique()), int(m.sum())

pw = 7.2 / 5; ph = pw * (b1[3]-b1[1]) / (b1[2]-b1[0]); leg_h = 0.92; tit_h = 0.22
fig, axs = plt.subplots(1, 5, figsize=(7.2, ph + leg_h + tit_h))
fig.subplots_adjust(left=0.004, right=0.996, top=1 - tit_h/(ph+leg_h+tit_h), bottom=leg_h/(ph+leg_h+tit_h), wspace=0.035)
nsp = sp.analysisSpecies.nunique()
for i, (ax, (L, t, g)) in enumerate(zip(axs[:4], types)):
    basemap(ax, b1, lab=(i == 0), bar=False, lab_lon="out")
    draw_pa(ax, g, TYPE_COL[t])
    draw_pts(ax, ins[t], C["blue"], z=5)        # A-D: only records inside
    s_, n_ = stats(ins[t])
    panel_title(ax, L, t)
    stats_text(ax, s_, n_)
ax = axs[4]; basemap(ax, b1, lab=False, bar=False, lab_lon="out")
draw_pa(ax, official, E_OFF, z=2.5)    # E: official protected areas only (priority sites are not protected areas, Law 21.600 art. 3 No. 31)
status = np.where(in_off, "off", np.where(in_pri, "pri", "out"))
draw_pts(ax, ~in_off, C["orange"], z=4); draw_pts(ax, in_off, C["blue"], z=5)
s_, n_ = stats(in_off)
panel_title(ax, "E", "Protected areas (A–C)")
stats_text(ax, s_, n_)
y35 = gpd.GeoSeries([Point(-67, -35)], crs=4326).to_crs(AEA).iloc[0].y
scalebar(axs[0], b1, km=300, frac_x=0.6, frac_y=(y35 - b1[1]) / (b1[3] - b1[1]), mid=False)
mk = lambda c, lab: Line2D([], [], marker="o", ls="", markersize=7, markerfacecolor=c, markeredgecolor="white", markeredgewidth=0.3, label=lab)
h = [Patch(facecolor=TYPE_COL[t], edgecolor="none", label=f"{L}  {t}") for L, t, _ in types]
h += [Patch(facecolor=E_OFF, edgecolor="none", label="E  Protected areas (A–C combined)"),
      mk(C["blue"], "Record inside the area shown"),
      mk(C["orange"], "Record outside protected areas (E)")]
lg = fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, leg_h/(ph+leg_h+tit_h) - 0.025), ncol=3, frameon=False,
                handlelength=1.5, handleheight=1.1, columnspacing=1.4, fontsize=8.2)
save(fig, "Fig4_protected_areas")
f8 = {t: dict(zip(["species", "records"], stats(ins[t]))) for t in ins}
f8.update(official=dict(zip(["species", "records"], stats(in_off))),
          any=dict(zip(["species", "records"], stats(in_off | in_pri))), species_total=int(nsp), records_total=len(pts))
print(f8)

# ---- Table S3: records of each species in each type of area (continental Chile)
sp_idx = pts.index[pts.analysisSpecies.notna()]
t3 = pd.DataFrame({"Species": pts.loc[sp_idx, "analysisSpecies"].values,
                   "State protected areas": ins["State protected areas"][pts.index.get_indexer(sp_idx)],
                   "Nature sanctuaries": ins["Nature sanctuaries"][pts.index.get_indexer(sp_idx)],
                   "Private protected areas": ins["Private protected areas"][pts.index.get_indexer(sp_idx)],
                   "Priority sites": ins["Priority sites"][pts.index.get_indexer(sp_idx)],
                   "Official (any)": in_off[pts.index.get_indexer(sp_idx)]})
tS3 = t3.groupby("Species").agg(Records=("Species", "size"), **{c: (c, "sum") for c in t3.columns[1:]}).reset_index()
tS3 = tS3.sort_values(["Records", "Species"], ascending=[False, True])
tS3["Species"] = "S. " + tS3.Species.str.replace("^Syntrichia ", "", regex=True)
tS3.to_csv(f"{OUT}/TableS3_species_in_protected_areas.csv", index=False)
print(tS3.head(8).to_string(index=False))

summary.update(fig4=f8)
json.dump(summary, open(f"{OUT}/summary.json", "w"), indent=1)
grid.drop(columns="geometry").to_csv(f"{OUT}/hex_grid_50km_summary.csv", index=False)
print("OK")
