# Great Lake Contours

Static layers for Google Earth, served from GitHub Pages:

```
https://wyattleathorn7.github.io/Great_Lake_Contours-/
```

## Great Lakes Depth Gradient

A continuous bathymetry colour gradient over all five Great Lakes with **all
land and every island cut out** (transparent), so the lakebed shows through
instead of a flat blue blob.

| File | Use |
| --- | --- |
| `Great_Lakes_Depth_Gradient.kml` | Entry file. Open this URL in Google Earth Web or Google Earth Pro. |
| `Great_Lakes_Depth_Gradient.kmz` | Same layer with the images embedded, for offline use. |
| `depth_gradient/depth_<lake>.png` | The five per-lake overlays. |
| `depth_gradient/legend.png` | Depth scale. |
| `scripts/build_depth_gradient.py` | Rebuilds all of the above. |

Direct link to the entry file:

```
https://wyattleathorn7.github.io/Great_Lake_Contours-/Great_Lakes_Depth_Gradient.kml
```

### Why this is built the way it is

Google Earth Web's file importer rejects `NetworkLink`, `ScreenOverlay` and
`Region`, which is why a tiled, multi-resolution KMZ will not import. This layer
therefore uses nothing but five plain `GroundOverlay` elements:

- One overlay per lake, so no single image has to span the whole basin.
- Each image is 4096 px on its long side, which is 0.076-0.150 km/px - finer
  than a single basin-wide sheet, and still inside Earth Web's texture budget.
- Every `href` is an absolute `https://wyattleathorn7.github.io/...` URL.
- No `Region`, so there is no level-of-detail machinery to fail.

### Method

- **Bathymetry:** NOAA Great Lakes LLD 3 arc-second grids
  (`grid_<lake>/<lake>_lld.tif`, EPSG:4269). Lakebed is negative, land is
  positive; positive cells are treated as 0 m so the shoreline reads as 0 m
  rather than being punched out by nodata.
- **Shoreline:** GSHHG v2.3.7 full-resolution. L2 lake polygons, with all L3
  islands removed and L4 ponds restored as water.
- **Colour ramp:** pale cyan at 0 m to navy at 410 m, drawn at 82% opacity.
- **Ownership:** where two LLD grids overlap on water, the pixel is given to one
  grid only, so overlapping overlays never stack alpha into a dark seam. Only
  three overlaps exist and each is assigned to the authoritative survey:

  | Overlap | Assigned to |
  | --- | --- |
  | Straits of Mackinac (superior/michigan) | superior |
  | St. Marys River (superior/huron) | superior |
  | Niagara / Lake Ontario west end (huron/ontario) | ontario |

  The last one matters because the Huron grid carries bogus 84-103 m values in
  its south-east corner, where Lake Ontario is only about 10-30 m deep.

### Numbers

| Lake | Image | Max depth found | PNG |
| --- | --- | --- | --- |
| Superior | 4096 x 2583 | 404.6 m | 1.8 MB |
| Michigan | 2328 x 4096 | 275.8 m | 1.1 MB |
| Huron | 4033 x 4096 | 224.1 m | 1.9 MB |
| Erie | 4096 x 1825 | 62.5 m | 0.2 MB |
| Ontario | 4096 x 1609 | 243.8 m | 0.7 MB |

Every maximum matches the published NOAA figure for the lake. Checks against the
full-resolution GSHHG mask:

- water not covered: 0.16% (the Thousand Islands, north-east of Lake Ontario,
  which is outside the NOAA grid extent)
- covered but not water: 0.94% (antialiasing along the shoreline)
- pixels drawn by two overlays: ~10
- islands verified transparent: Isle Royale, Beaver Island, Manitoulin,
  Mackinac Island, Pelee Island, St. Joseph Island

Rebuild with:

```sh
python3 scripts/build_depth_gradient.py
```

Requires the `NOAA-Nautical-Overlay` data tree as a sibling directory of this
repository, plus numpy, rasterio, shapely, pyshp and Pillow.

## Other layers

- `Lake_Huron_Complete_Contours.kml` - network links to the 28 Lake Huron
  contour KMZ parts (white contour lines, depths in feet).
- `Great_Lakes_Contours_XYZ/` - XYZ raster tile pyramid of the contour lines.
