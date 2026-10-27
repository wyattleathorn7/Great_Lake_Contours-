"""Build the GitHub-hosted Great Lakes depth-gradient layer.

Google Earth Web's file-import validator rejects NetworkLink, ScreenOverlay and
Region. This build emits the same bathymetry as a flat Document of plain
GroundOverlays pointing at absolute GitHub Pages URLs, so it loads by URL and
imports as a file without any unsupported elements.

One GroundOverlay per lake, each rendered at 4096 px on its long side over the
lake's own NOAA LLD grid extent. That is finer than a single basin-wide sheet
(0.076-0.150 km/px) and needs no Regions at all.

Land and every island are cut out from the GSHHG full-resolution shoreline.
Each water pixel is drawn by exactly one overlay (see OWNER rule below) so
overlapping grid extents never stack alpha into a dark seam.

Outputs (repo root):
    Great_Lakes_Depth_Gradient.kml   entry file, absolute Pages URLs
    Great_Lakes_Depth_Gradient.kmz   same doc with images embedded, relative
    depth_gradient/*.png             the overlays
"""

import os
import shutil
import sys
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import rasterio
import shapefile
from PIL import Image, ImageDraw, ImageFont
from rasterio.crs import CRS
from rasterio.features import rasterize
from rasterio.transform import from_bounds
from rasterio.warp import Resampling, reproject
from shapely.geometry import Point, box, shape
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BASE = os.path.dirname(ROOT)

GSHHG = os.path.join(BASE, "NOAA-Nautical-Overlay", "bathy", "gshhg",
                     "GSHHS_shp", "f")
BATHY = os.path.join(BASE, "NOAA-Nautical-Overlay", "bathy")
IMGDIR = os.path.join(ROOT, "depth_gradient")

PAGES = "https://wyattleathorn7.github.io/Great_Lake_Contours-/"

LAKES = ["superior", "michigan", "huron", "erie", "ontario"]
LABEL = {
    "superior": "Lake Superior",
    "michigan": "Lake Michigan",
    "huron": "Lake Huron",
    "erie": "Lake Erie",
    "ontario": "Lake Ontario",
}

# Lower index wins where grid extents overlap on water. Only three overlaps
# exist, and each is assigned to the authoritative grid:
#   superior/michigan  Straits of Mackinac  -> superior
#   superior/huron     St. Marys River      -> superior
#   huron/ontario      Niagara / Ontario W  -> ontario (the Huron grid carries
#                       bogus 84-103 m values in that corner; Lake Ontario's
#                       west basin is only ~10-30 m deep)
OWNER = ["superior", "michigan", "ontario", "huron", "erie"]

LONG_SIDE = 4096
OVERLAY_ALPHA = 220
DEPTH_MIN, DEPTH_MAX = 0.0, 410.0

# One scale for every lake: dark blue = 0 m (shore), dark red = 410 m (the
# deepest point in the basin). A shallow lake simply never reaches the warm end,
# so Erie (62 m) stops at green while Superior runs through orange to red.
BATHY_STOPS = [
    (0.00, (10, 22, 72)),      # darkest navy - shoreline
    (0.03, (14, 44, 122)),
    (0.06, (19, 72, 168)),
    (0.09, (24, 102, 196)),
    (0.12, (29, 133, 208)),    # blue ends here: 0-49 m
    (0.17, (34, 163, 206)),    # cyan
    (0.24, (40, 187, 180)),    # teal
    (0.32, (58, 198, 140)),    # sea green
    (0.41, (95, 202, 88)),     # green
    (0.50, (155, 203, 52)),    # yellow-green
    (0.59, (208, 198, 38)),    # yellow
    (0.67, (231, 172, 27)),    # amber
    (0.74, (236, 128, 20)),    # orange
    (0.81, (226, 78, 16)),     # orange-red
    (0.89, (186, 30, 12)),     # red
    (0.96, (143, 4, 9)),
    (1.00, (104, 0, 7)),       # dark red - deepest
]

SEEDS = [
    ("Superior", -88.0, 47.5),
    ("Michigan", -87.0, 44.2),
    ("Huron", -82.5, 44.8),
    ("Erie", -81.5, 42.2),
    ("Ontario", -77.8, 43.6),
    ("StClair", -82.45, 42.43),
]

KML_NS = "http://www.opengis.net/kml/2.2"


def log(*a):
    print(*a, flush=True)


def grid_path(name):
    return os.path.join(BATHY, f"grid_{name}", f"{name}_lld", f"{name}_lld.tif")


# ----------------------------------------------------------------- geometry
def _records(shp):
    return list(shapefile.Reader(shp))


def _geom(sr):
    g = shape(sr.shape.__geo_interface__)
    return g if g.is_valid else g.buffer(0)


def load_water():
    """GSHHG Great Lakes water with every island and pond cut out."""
    lakes = []
    for sr in _records(os.path.join(GSHHG, "GSHHS_f_L2.shp")):
        g = _geom(sr)
        if g.is_empty:
            continue
        if any(g.contains(Point(x, y)) for _n, x, y in SEEDS):
            lakes.append(g)
    if not lakes:
        raise RuntimeError("no GSHHG L2 polygon contains a lake seed")
    water = unary_union(lakes)
    log(f"  lake polygons: {len(lakes)}  parts: {len(list(water.geoms))}")

    islands = []
    for sr in _records(os.path.join(GSHHG, "GSHHS_f_L3.shp")):
        g = _geom(sr)
        if g.is_empty:
            continue
        c = g.centroid
        if c.x < -95 or c.x > -73:
            continue
        if water.contains(c):
            islands.append(g)
    log(f"  islands cut out: {len(islands)}")
    if islands:
        water = water.difference(unary_union(islands))

    ponds = []
    try:
        for sr in _records(os.path.join(GSHHG, "GSHHS_f_L4.shp")):
            g = _geom(sr)
            if g.is_empty:
                continue
            c = g.centroid
            if -95 <= c.x <= -73 and water.contains(c):
                ponds.append(g)
    except Exception as e:
        log(f"  L4 skipped: {e}")
    if ponds:
        log(f"  ponds restored as water: {len(ponds)}")
        water = unary_union([water] + ponds)

    if not water.is_valid:
        water = water.buffer(0)
    return water


# ------------------------------------------------------------------ raster
def lake_extent(name):
    ds = rasterio.open(grid_path(name))
    b = ds.bounds
    ext = (b.left, b.bottom, b.right, b.top)
    ds.close()
    return ext


def pixel_dims(w, s, e, n):
    lat_mid = 0.5 * (s + n)
    km_x = (e - w) * 111.32 * np.cos(np.radians(lat_mid))
    km_y = (n - s) * 110.57
    if km_x >= km_y:
        pw, ph = LONG_SIDE, max(2, int(round(LONG_SIDE * km_y / km_x)))
    else:
        pw, ph = max(2, int(round(LONG_SIDE * km_x / km_y))), LONG_SIDE
    return pw, ph


def coverage(water, bbox, pw, ph, ss=2):
    """Water coverage fraction per pixel, 2x supersampled for clean edges."""
    w, s, e, n = bbox
    tr = from_bounds(w, s, e, n, pw * ss, ph * ss)
    m = rasterize([(water, 1)], out_shape=(ph * ss, pw * ss), transform=tr,
                  fill=0, default_value=1, all_touched=False, dtype=np.uint8)
    q = (m[0::ss, 0::ss].astype(np.uint16) + m[0::ss, 1::ss] +
         m[1::ss, 0::ss] + m[1::ss, 1::ss])
    return (q / float(ss * ss)).astype(np.float32)


def depth_and_owner(bbox, pw, ph):
    """Depth in metres plus the lake index that owns each pixel.

    The LLD grids encode lakebed as negative values and land as positive ones.
    Positive cells are flattened to 0 m rather than dropped as nodata: that is
    the usual bathymetric convention (0 m at the shoreline), it stops bilinear
    resampling from eroding a NaN halo inland along every shore, and it lets
    the GSHHG shoreline - not the survey mask - define the water edge.
    """
    w, s, e, n = bbox
    dst_tr = from_bounds(w, s, e, n, pw, ph)
    depth = np.zeros((ph, pw), dtype=np.float32)
    owner = np.full((ph, pw), -1, dtype=np.int16)

    for i, name in enumerate(OWNER):
        gw, gs, ge, gn = lake_extent(name)
        if ge <= w or gw >= e or gn <= s or gs >= n:
            continue
        ds = rasterio.open(grid_path(name))
        v = ds.read(1).astype(np.float32)
        nod = -9999.0 if ds.nodata is None else float(ds.nodata)
        src = np.where((v < 0) & (v != nod), v, 0.0)
        dst = np.zeros((ph, pw), dtype=np.float32)
        reproject(src, dst,
                  src_transform=ds.transform, src_crs=ds.crs,
                  dst_transform=dst_tr, dst_crs=CRS.from_epsg(4326),
                  resampling=Resampling.bilinear)
        ds.close()
        take = (dst < 0) & (owner < 0)
        depth[take] = -dst[take]
        owner[take] = i

    # Shoreline and unsurveyed water carry no depth, so give every remaining
    # pixel to the first grid in OWNER whose extent contains it. That rule is a
    # pure function of position, so all five renders agree on who draws what.
    lons = w + (np.arange(pw) + 0.5) * (e - w) / pw
    lats = n - (np.arange(ph) + 0.5) * (n - s) / ph
    for i, name in enumerate(OWNER):
        gw, gs, ge, gn = lake_extent(name)
        if ge <= w or gw >= e or gn <= s or gs >= n:
            continue
        ci = (lons >= gw) & (lons <= ge)
        ri = (lats >= gs) & (lats <= gn)
        if not ci.any() or not ri.any():
            continue
        owner[ri[:, None] & ci[None, :] & (owner < 0)] = i
    return depth, owner


def lut(n=256):
    pos = np.array([s[0] for s in BATHY_STOPS])
    col = np.array([s[1] for s in BATHY_STOPS], dtype=np.float64)
    xs = np.linspace(0.0, 1.0, n)
    out = np.empty((n, 3), dtype=np.uint8)
    for ch in range(3):
        out[:, ch] = np.clip(np.interp(xs, pos, col[:, ch]), 0, 255).astype(np.uint8)
    return out


def bleed(rgba, iterations=8):
    """Push water colour outward under transparent pixels (no dark fringes)."""
    rgb = rgba[:, :, :3].astype(np.float32)
    filled = rgba[:, :, 3] > 0
    if filled.all():
        return rgba
    out = rgb.copy()
    for _ in range(iterations):
        if filled.all():
            break
        f = filled.astype(np.float32)
        tot = np.zeros_like(out)
        cnt = np.zeros(filled.shape, dtype=np.float32)
        tot[1:, :, :] += out[:-1, :, :] * f[:-1, :, None]
        cnt[1:, :] += f[:-1, :]
        tot[:-1, :, :] += out[1:, :, :] * f[1:, :, None]
        cnt[:-1, :] += f[1:, :]
        tot[:, 1:, :] += out[:, :-1, :] * f[:, :-1, None]
        cnt[:, 1:] += f[:, :-1]
        tot[:, :-1, :] += out[:, 1:, :] * f[:, 1:, None]
        cnt[:, :-1] += f[:, 1:]
        grow = (~filled) & (cnt > 0)
        out[grow] = tot[grow] / cnt[grow, None]
        filled = filled | grow
    res = rgba.copy()
    res[:, :, :3] = np.round(out).astype(np.uint8)
    return res


def render_lake(name, idx, water):
    bbox = lake_extent(name)
    pw, ph = pixel_dims(*bbox)
    depth, owner = depth_and_owner(bbox, pw, ph)
    mine = owner == OWNER.index(name)
    cov = coverage(water, bbox, pw, ph)

    l = lut()
    rgba = np.zeros((ph, pw, 4), dtype=np.uint8)
    valid = mine & (cov > 0)
    t = np.zeros((ph, pw), dtype=np.float32)
    t[valid] = np.clip((depth[valid] - DEPTH_MIN) / (DEPTH_MAX - DEPTH_MIN), 0, 1)
    i = (t * 255).astype(np.int32)
    rgba[valid, :3] = l[i[valid]]
    rgba[valid, 3] = np.round(OVERLAY_ALPHA * cov[valid]).astype(np.uint8)
    rgba = bleed(rgba)

    png = os.path.join(IMGDIR, f"depth_{name}.png")
    Image.fromarray(rgba, "RGBA").save(png, optimize=True)
    dmax = float(depth[valid].max()) if valid.any() else float("nan")
    log(f"  {name:9s} {pw}x{ph}  water px {int(valid.sum()):,}  "
        f"max depth {dmax:.1f} m  {os.path.getsize(png)/1e3:.0f} kB")
    return bbox


# --------------------------------------------------------------------- KML
def _f(v):
    return f"{v:.6f}"


def kml_doc(prefix, with_images):
    desc = description_html(with_images)
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           f'<kml xmlns="{KML_NS}">',
           '<Document>',
           '<name>Great Lakes Depth Gradient</name>',
           f'<description><![CDATA[{desc}]]></description>']
    if with_images:
        out.append('<Style><ListStyle><listItemType>check</listItemType>'
                   '</ListStyle></Style>')
    for name in LAKES:
        w, s, e, n = EXTENTS[name]
        out += [
            '<GroundOverlay>',
            f'<name>Depth - {LABEL[name]}</name>',
            '<description>NOAA LLD 3-arcsecond bathymetry as a continuous depth '
            'gradient. Land, islands and ponds cut out (transparent).</description>',
            '<visibility>1</visibility>',
            '<drawOrder>40</drawOrder>',
            '<Icon>',
            f'<href>{prefix}depth_gradient/depth_{name}.png</href>',
            '</Icon>',
            '<LatLonBox>',
            f'<north>{_f(n)}</north><south>{_f(s)}</south>'
            f'<east>{_f(e)}</east><west>{_f(w)}</west>',
            '</LatLonBox>',
            '</GroundOverlay>',
        ]
    out += ['</Document>', '</kml>']
    return "\n".join(out)


def description_html(with_images):
    l = lut()
    sw = []
    for v in range(0, 401, 50):
        c = l[int(round((v - DEPTH_MIN) / (DEPTH_MAX - DEPTH_MIN) * 255))]
        sw.append(f'<span style="display:inline-block;width:14px;height:14px;'
                  f'background:rgb({c[0]},{c[1]},{c[2]});border:1px solid #888;'
                  f'vertical-align:middle"></span> {v} m')
    img = ('<p><img src="' + PAGES + 'depth_gradient/legend.png" '
           'width="420"></p>') if with_images else ''
    return (
        '<h2>Great Lakes Depth Gradient</h2>'
        '<p>Bathymetry rendered as a continuous colour gradient. All land and '
        'every island are cut out, so the overlay shows water only.</p>'
        '<p><b>Depth (metres) - dark blue is shallow, dark red is deep:</b><br>'
        + " &nbsp; ".join(sw) + '</p>'
        '<p>One shared scale for all five lakes, so colour means the same depth '
        'everywhere: Erie (62 m maximum) never leaves the blue-green range, '
        'while Superior reaches dark red at 405 m.</p>'
        + img +
        '<p><b>Sources:</b> NOAA Great Lakes LLD bathymetric grids '
        '(3 arc-second) and GSHHG v2.3.7 full-resolution shoreline '
        '(lakes L2, islands L3, ponds L4).</p>'
        '<p>Five overlays, one per lake, each at up to 4096 px on the long '
        'side. Built from plain image overlays only, so it loads in both '
        'Google Earth Web and Google Earth Pro.</p>')


# ------------------------------------------------------------------ legend
def _font(size, bold=False):
    for c in ("/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold
              else "/System/Library/Fonts/Supplemental/Arial.ttf",
              "/System/Library/Fonts/Helvetica.ttc",
              "/System/Library/Fonts/Supplemental/Verdana.ttf"):
        if os.path.exists(c):
            try:
                return ImageFont.truetype(c, size)
            except Exception:
                pass
    return ImageFont.load_default()


def make_legend(path):
    W, H = 420, 186
    img = Image.new("RGBA", (W, H), (255, 255, 255, 235))
    d = ImageDraw.Draw(img)
    d.text((14, 12), "Great Lakes Depth Gradient", font=_font(20, True),
           fill=(20, 30, 60, 255))
    d.text((14, 38), "Lake depth below surface (metres)", font=_font(13),
           fill=(60, 70, 90, 255))
    x0, y0, x1, y1 = 16, 74, W - 16, 100
    l = lut()
    for i in range(x1 - x0):
        c = l[int(round(i / (x1 - x0 - 1) * 255))]
        d.line([(x0 + i, y0), (x0 + i, y1)],
               fill=(int(c[0]), int(c[1]), int(c[2]), 255))
    fs = _font(12)
    for v in range(0, 401, 50):
        x = x0 + v / DEPTH_MAX * (x1 - x0)
        d.line([(x, y1), (x, y1 + 6)], fill=(40, 50, 70, 255))
        tw = d.textlength(str(v), font=fs)
        d.text((x - tw / 2, y1 + 8), str(v), font=fs, fill=(40, 50, 70, 255))
    d.text((16, 132), "Dark blue = shallow, dark red = deep (0-410 m, all lakes)",
           font=fs, fill=(60, 70, 90, 255))
    d.text((16, 148), "NOAA LLD 3-arcsec bathymetry + GSHHG full-res shoreline",
           font=fs, fill=(60, 70, 90, 255))
    d.text((16, 164), "Land and every island cut out (transparent).",
           font=fs, fill=(60, 70, 90, 255))
    img.save(path)


# -------------------------------------------------------------------- main
EXTENTS = {}


def _iter(root, tag):
    return list(root.iter("{http://www.opengis.net/kml/2.2}" + tag))


def _hrefs(root):
    return [h.text for h in _iter(root, "href")]


def _check_earth_web(root, where):
    bad = [t for t in ("NetworkLink", "ScreenOverlay", "Region")
           if _iter(root, t)]
    if bad:
        raise RuntimeError(f"{where}: Earth Web rejects {bad}")
    n = len(_iter(root, "GroundOverlay"))
    if n != len(LAKES):
        raise RuntimeError(f"{where}: expected {len(LAKES)} overlays, got {n}")


def verify_kml(path):
    """No forbidden elements, and only absolute Pages URLs."""
    root = ET.parse(path).getroot()
    _check_earth_web(root, os.path.basename(path))
    for h in _hrefs(root):
        if not h.startswith(PAGES):
            raise RuntimeError(f"{os.path.basename(path)}: href not absolute: {h}")
    log(f"  verified {os.path.basename(path)}: {len(LAKES)} overlays, "
        f"0 forbidden elements, all hrefs absolute")


def verify_kmz(path):
    """Every relative href must resolve inside the archive, as a real PNG."""
    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
        root = ET.fromstring(z.read("doc.kml"))
        _check_earth_web(root, os.path.basename(path))
        hrefs = _hrefs(root)
        if not hrefs:
            raise RuntimeError("kmz: no image hrefs")
        for h in hrefs:
            if h.startswith("http"):
                raise RuntimeError(f"kmz: href should be relative, got {h}")
            if h not in names:
                raise RuntimeError(f"kmz: href {h} is not in the archive")
            if z.read(h)[:8] != b"\x89PNG\r\n\x1a\n":
                raise RuntimeError(f"kmz: {h} is not a PNG")
    log(f"  verified {os.path.basename(path)}: all {len(hrefs)} images resolve "
        f"inside the archive")


def main():
    os.makedirs(IMGDIR, exist_ok=True)
    log("== Great Lakes depth gradient (GitHub-hosted, Earth Web clean) ==")

    log("[1] water geometry (GSHHG L2 - L3 + L4)")
    water = load_water()

    log("[2] render per-lake overlays at 4096 px long side")
    for idx, name in enumerate(LAKES):
        EXTENTS[name] = render_lake(name, idx, water)

    log("[3] legend + KML")
    make_legend(os.path.join(IMGDIR, "legend.png"))
    kml = os.path.join(ROOT, "Great_Lakes_Depth_Gradient.kml")
    with open(kml, "w") as f:
        f.write(kml_doc(PAGES, with_images=True))
    verify_kml(kml)
    log(f"  {kml}  ({os.path.getsize(kml)/1e3:.1f} kB)")

    log("[4] offline KMZ (embedded images, relative paths)")
    kmz = os.path.join(ROOT, "Great_Lakes_Depth_Gradient.kmz")
    tmp = os.path.join(IMGDIR, "_kmz")
    shutil.rmtree(tmp, ignore_errors=True)
    # The KMZ mirrors the repo layout, so doc.kml keeps the same
    # "depth_gradient/<file>.png" hrefs as the hosted KML and resolves them
    # against the archive instead of the web.
    os.makedirs(os.path.join(tmp, "depth_gradient"))
    with open(os.path.join(tmp, "doc.kml"), "w") as f:
        f.write(kml_doc("", with_images=False))
    entries = ["doc.kml"]
    for name in LAKES + ["legend"]:
        src = f"depth_{name}.png" if name != "legend" else "legend.png"
        arc = f"depth_gradient/{src}"
        shutil.copy(os.path.join(IMGDIR, src), os.path.join(tmp, arc))
        entries.append(arc)
    with zipfile.ZipFile(kmz, "w", zipfile.ZIP_DEFLATED) as z:
        for rel in entries:
            z.write(os.path.join(tmp, rel), rel)
    shutil.rmtree(tmp)
    verify_kmz(kmz)
    log(f"  {kmz}  ({os.path.getsize(kmz)/1e6:.2f} MB)")

    log("done.")


if __name__ == "__main__":
    sys.exit(main())
