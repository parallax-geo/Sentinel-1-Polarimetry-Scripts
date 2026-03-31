#!/usr/bin/env python3
"""
STEP 2: Preprocess Sentinel-1 SLC → C2 Covariance Matrix (v9 ALL IW)
======================================================================

STRATEGY: Process each sub-swath (IW1, IW2, IW3) SEPARATELY through SNAP
(which works reliably), then MOSAIC the 3 geocoded outputs into a single
seamless GeoTIFF using GDAL.

This avoids the SNAP "Cannot construct DataBuffer" / OutOfMemory error
when trying to merge 3 sub-swaths inside a single SNAP graph.

PER SUB-SWATH — SNAP GPT Graph (confirmed working):
┌──────────────────────────────────────────────────────────────┐
│  Read SLC .zip                                               │
│    → Apply-Orbit-File                                        │
│    → TOPSAR-Split (IWx, VV+VH)                               │
│    → Calibration (outputImageInComplex=true)                  │
│    → TOPSAR-Deburst                                           │
│    → Polarimetric-Matrices (C2)                               │
│    → Multilook (5az × 1rg)                                    │
│    → Polarimetric-Speckle-Filter (Refined Lee 5×5)            │
│    → Write BEAM-DIMAP                                         │
└──────────────────────────────────────────────────────────────┘
  × 3 (IW1, IW2, IW3)

GEOCODING — GDAL (per sub-swath):
┌──────────────────────────────────────────────────────────────┐
│  Read tie-point grids (latitude.img, longitude.img)          │
│    → Build GCPs → gdalwarp TPS → EPSG:4326                  │
└──────────────────────────────────────────────────────────────┘
  × 3

MOSAIC — GDAL (combine 3 sub-swaths):
┌──────────────────────────────────────────────────────────────┐
│  gdal.Warp([IW1.tif, IW2.tif, IW3.tif]) → merged.tif       │
│    → gdal.Translate clip to AOI                              │
│    → Final: C11.tif, C12_real.tif, C12_imag.tif, C22.tif    │
└──────────────────────────────────────────────────────────────┘

Usage:
    conda activate intnotebook
    python step2_preprocess_slc.py
"""

import os
import sys
import glob
import shutil
import subprocess
import re
import numpy as np
import geopandas as gpd
from osgeo import gdal, osr
gdal.UseExceptions()

# =========================================================================
# CONFIGURATION
# =========================================================================
WORK_DIR = os.path.expanduser("~/Desktop/Muddasir")
DOWNLOADS_DIR = os.path.join(WORK_DIR, "downloads")
AOI_SHP = os.path.join(WORK_DIR, "aoi", "Battambang.shp")
OUTPUT_DIR = os.path.join(WORK_DIR, "outputs", "preprocessed")
os.makedirs(OUTPUT_DIR, exist_ok=True)

GPT_PATH = shutil.which("gpt")
if GPT_PATH is None:
    for candidate in [
        os.path.expanduser("~/esa-snap/bin/gpt"),
        os.path.expanduser("~/snap/bin/gpt"),
        "/usr/local/snap/bin/gpt",
        "/opt/snap/bin/gpt",
    ]:
        if os.path.isfile(candidate):
            GPT_PATH = candidate
            break
if GPT_PATH is None:
    print("❌ ERROR: SNAP GPT not found!"); sys.exit(1)
print(f"Using GPT: {GPT_PATH}")

ALL_SUBSWATHS = ["IW1", "IW2", "IW3"]
AZIMUTH_LOOKS = 5
RANGE_LOOKS = 1
PIXEL_SPACING_M = 30.0
PIXEL_SPACING_DEG = PIXEL_SPACING_M / 111000.0
GPT_MEMORY = "8G"


def get_aoi_bounds(shp_path):
    aoi = gpd.read_file(shp_path)
    if aoi.crs and aoi.crs.to_epsg() != 4326:
        aoi = aoi.to_crs(epsg=4326)
    bounds = aoi.total_bounds
    buf = 0.05
    return (bounds[0]-buf, bounds[1]-buf, bounds[2]+buf, bounds[3]+buf)


# =========================================================================
# SNAP GRAPH — single sub-swath (confirmed working)
# =========================================================================
def create_graph_single_iw(slc_zip, output_path, subswath):
    """Exact same graph that worked for IW1 — just parameterized."""
    return f"""<graph id="Graph1_SLC_to_C2_{subswath}">
  <version>1.0</version>
  <node id="Read">
    <operator>Read</operator>
    <sources/>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <file>{slc_zip}</file>
    </parameters>
  </node>
  <node id="Apply-Orbit-File">
    <operator>Apply-Orbit-File</operator>
    <sources><sourceProduct refid="Read"/></sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <orbitType>Sentinel Precise (Auto Download)</orbitType>
      <polyDegree>3</polyDegree>
      <continueOnFail>true</continueOnFail>
    </parameters>
  </node>
  <node id="TOPSAR-Split">
    <operator>TOPSAR-Split</operator>
    <sources><sourceProduct refid="Apply-Orbit-File"/></sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <subswath>{subswath}</subswath>
      <selectedPolarisations>VV,VH</selectedPolarisations>
    </parameters>
  </node>
  <node id="Calibration">
    <operator>Calibration</operator>
    <sources><sourceProduct refid="TOPSAR-Split"/></sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <selectedPolarisations>VV,VH</selectedPolarisations>
      <outputSigmaBand>false</outputSigmaBand>
      <outputBetaBand>false</outputBetaBand>
      <outputGammaBand>false</outputGammaBand>
      <outputImageInComplex>true</outputImageInComplex>
    </parameters>
  </node>
  <node id="TOPSAR-Deburst">
    <operator>TOPSAR-Deburst</operator>
    <sources><sourceProduct refid="Calibration"/></sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <selectedPolarisations>VV,VH</selectedPolarisations>
    </parameters>
  </node>
  <node id="Polarimetric-Matrices">
    <operator>Polarimetric-Matrices</operator>
    <sources><sourceProduct refid="TOPSAR-Deburst"/></sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <matrix>C2</matrix>
    </parameters>
  </node>
  <node id="Multilook">
    <operator>Multilook</operator>
    <sources><sourceProduct refid="Polarimetric-Matrices"/></sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <nAzLooks>{AZIMUTH_LOOKS}</nAzLooks>
      <nRgLooks>{RANGE_LOOKS}</nRgLooks>
      <outputIntensity>false</outputIntensity>
    </parameters>
  </node>
  <node id="Polarimetric-Speckle-Filter">
    <operator>Polarimetric-Speckle-Filter</operator>
    <sources><sourceProduct refid="Multilook"/></sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <filter>Refined Lee Filter</filter>
      <filterSize>5</filterSize>
      <numLooksStr>1</numLooksStr>
    </parameters>
  </node>
  <node id="Write">
    <operator>Write</operator>
    <sources><sourceProduct refid="Polarimetric-Speckle-Filter"/></sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <file>{output_path}</file>
      <formatName>BEAM-DIMAP</formatName>
    </parameters>
  </node>
</graph>"""


def run_gpt(graph_file, label):
    cmd = [GPT_PATH, graph_file, f"-J-Xmx{GPT_MEMORY}", "-J-Xms2G", "-c", "4G"]
    print(f"      CMD: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        if result.returncode != 0:
            print(f"\n  ❌ {label} FAILED (code {result.returncode})")
            stderr = result.stderr[-3000:] if result.stderr else ""
            print(f"  STDERR:\n{stderr}")
            return False
        print(f"  ✅ {label} completed!")
        return True
    except subprocess.TimeoutExpired:
        print(f"  ❌ {label} timed out"); return False
    except Exception as e:
        print(f"  ❌ {label} error: {e}"); return False


# =========================================================================
# GEOCODING — tie-point grids → GCPs → gdalwarp (confirmed working)
# =========================================================================
def read_tiepoint_grids(data_dir, dim_file):
    tp_dir = os.path.join(data_dir, "tie_point_grids")
    lat_file = lon_file = None
    if os.path.isdir(tp_dir):
        for f in os.listdir(tp_dir):
            fl = f.lower()
            if fl.startswith("latitude") and fl.endswith(".img"):
                lat_file = os.path.join(tp_dir, f)
            elif fl.startswith("longitude") and fl.endswith(".img"):
                lon_file = os.path.join(tp_dir, f)
    if not lat_file or not lon_file:
        print(f"    ❌ lat/lon not found in {tp_dir}"); return None
    ds_lat = gdal.Open(lat_file); ds_lon = gdal.Open(lon_file)
    if ds_lat is None or ds_lon is None: return None
    lat_arr = ds_lat.GetRasterBand(1).ReadAsArray().astype(np.float64)
    lon_arr = ds_lon.GetRasterBand(1).ReadAsArray().astype(np.float64)
    ds_lat = None; ds_lon = None
    tp_rows, tp_cols = lat_arr.shape
    print(f"    Tie-points: {tp_cols}×{tp_rows}, Lat: {lat_arr.min():.3f}→{lat_arr.max():.3f}, Lon: {lon_arr.min():.3f}→{lon_arr.max():.3f}")
    with open(dim_file, 'r') as f: content = f.read()
    offset_x = float(re.findall(r'<OFFSET_X>(.*?)</OFFSET_X>', content)[0])
    offset_y = float(re.findall(r'<OFFSET_Y>(.*?)</OFFSET_Y>', content)[0])
    step_x = float(re.findall(r'<STEP_X>(.*?)</STEP_X>', content)[0])
    step_y = float(re.findall(r'<STEP_Y>(.*?)</STEP_Y>', content)[0])
    return lat_arr, lon_arr, offset_x, offset_y, step_x, step_y


def build_gcps(lat_arr, lon_arr, ox, oy, sx, sy):
    tp_rows, tp_cols = lat_arr.shape
    gcps = []
    for r in range(tp_rows):
        for c in range(tp_cols):
            px = ox + c * sx; ln = oy + r * sy
            lo = float(lon_arr[r, c]); la = float(lat_arr[r, c])
            if np.isfinite(la) and np.isfinite(lo) and la != 0 and lo != 0:
                gcps.append(gdal.GCP(lo, la, 0.0, px, ln))
    print(f"    Built {len(gcps)} GCPs")
    return gcps


def geocode_band(img_file, out_path, gcps, gcp_wkt):
    """Geocode a single .img band to EPSG:4326 (no AOI clip yet)."""
    temp_gcps = out_path.replace(".tif", "_tmp_gcps.tif")
    ds_in = gdal.Open(img_file)
    data = ds_in.GetRasterBand(1).ReadAsArray().astype(np.float32)
    cols, rows = ds_in.RasterXSize, ds_in.RasterYSize
    ds_in = None
    driver = gdal.GetDriverByName("GTiff")
    ds_tmp = driver.Create(temp_gcps, cols, rows, 1, gdal.GDT_Float32)
    ds_tmp.SetGCPs(gcps, gcp_wkt)
    ds_tmp.GetRasterBand(1).WriteArray(data)
    ds_tmp.GetRasterBand(1).SetNoDataValue(0)
    ds_tmp.FlushCache(); ds_tmp = None; del data
    try:
        gdal.Warp(out_path, temp_gcps, options=gdal.WarpOptions(
            format="GTiff", dstSRS="EPSG:4326",
            xRes=PIXEL_SPACING_DEG, yRes=PIXEL_SPACING_DEG,
            resampleAlg="bilinear",
            creationOptions=["COMPRESS=LZW", "BIGTIFF=IF_SAFER"],
            dstNodata=0, tps=True))
        print(f"      ✅ Warped")
    except Exception as e:
        print(f"      ❌ Warp failed: {e}")
        if os.path.isfile(temp_gcps): os.remove(temp_gcps)
        return False
    if os.path.isfile(temp_gcps): os.remove(temp_gcps)
    return os.path.isfile(out_path)


def geocode_iw(data_dir, dim_file, out_dir):
    """Geocode all 4 C2 bands of one sub-swath."""
    os.makedirs(out_dir, exist_ok=True)
    tp = read_tiepoint_grids(data_dir, dim_file)
    if tp is None: return False
    lat, lon, ox, oy, sx, sy = tp
    gcps = build_gcps(lat, lon, ox, oy, sx, sy)
    if len(gcps) < 4: return False
    srs = osr.SpatialReference(); srs.ImportFromEPSG(4326)
    gcp_wkt = srs.ExportToWkt()
    c2_map = {}
    for f in sorted(os.listdir(data_dir)):
        if not f.endswith(".img"): continue
        fl = f.lower()
        if "c11" in fl: c2_map["C11"] = os.path.join(data_dir, f)
        elif "c12_real" in fl: c2_map["C12_real"] = os.path.join(data_dir, f)
        elif "c12_imag" in fl: c2_map["C12_imag"] = os.path.join(data_dir, f)
        elif "c22" in fl: c2_map["C22"] = os.path.join(data_dir, f)
    if len(c2_map) < 4:
        print(f"    ❌ Only {len(c2_map)} C2 bands found"); return False
    ok = 0
    for name in ["C11", "C12_real", "C12_imag", "C22"]:
        print(f"    Geocoding {name}...")
        if geocode_band(c2_map[name], os.path.join(out_dir, f"{name}.tif"), gcps, gcp_wkt):
            ok += 1
    return ok == 4


# =========================================================================
# MOSAIC — merge 3 geocoded sub-swaths into one seamless GeoTIFF
# =========================================================================
def mosaic_subswaths(iw_dirs, final_c2_dir, aoi_bounds):
    """
    For each C2 band, merge the 3 sub-swath GeoTIFFs into one mosaic,
    then clip to AOI. GDAL handles this trivially.
    """
    os.makedirs(final_c2_dir, exist_ok=True)
    minx, miny, maxx, maxy = aoi_bounds
    ok = 0

    for band in ["C11", "C12_real", "C12_imag", "C22"]:
        # Collect available sub-swath files for this band
        inputs = []
        for iw_dir in iw_dirs:
            f = os.path.join(iw_dir, f"{band}.tif")
            if os.path.isfile(f):
                inputs.append(f)

        if not inputs:
            print(f"  ❌ No files for {band}"); continue

        print(f"\n  Mosaicing {band}: {len(inputs)} sub-swath(s)...")

        temp_mosaic = os.path.join(final_c2_dir, f"{band}_mosaic_tmp.tif")
        final_path = os.path.join(final_c2_dir, f"{band}.tif")

        # Step 1: Merge all sub-swath TIFs into one
        try:
            gdal.Warp(temp_mosaic, inputs, options=gdal.WarpOptions(
                format="GTiff", dstSRS="EPSG:4326",
                xRes=PIXEL_SPACING_DEG, yRes=PIXEL_SPACING_DEG,
                resampleAlg="bilinear",
                creationOptions=["COMPRESS=LZW", "BIGTIFF=IF_SAFER"],
                dstNodata=0))
            print(f"    ✅ Merged {len(inputs)} sub-swaths")
        except Exception as e:
            print(f"    ❌ Mosaic failed: {e}"); continue

        # Step 2: Clip to AOI
        try:
            gdal.Translate(final_path, temp_mosaic,
                          options=gdal.TranslateOptions(
                              format="GTiff",
                              projWin=[minx, maxy, maxx, miny],
                              projWinSRS="EPSG:4326",
                              creationOptions=["COMPRESS=LZW"],
                              noData=0))
            print(f"    ✅ Clipped to AOI")
        except Exception as e:
            print(f"    ⚠️ Clip failed ({e}), keeping full mosaic")
            shutil.move(temp_mosaic, final_path)

        if os.path.isfile(temp_mosaic): os.remove(temp_mosaic)

        if os.path.isfile(final_path):
            ds = gdal.Open(final_path)
            if ds:
                sz = os.path.getsize(final_path) / (1024*1024)
                print(f"    📐 {ds.RasterXSize} × {ds.RasterYSize} px, {sz:.1f} MB")
                ds = None
            ok += 1

    return final_c2_dir if ok == 4 else None


# =========================================================================
# PROCESS ONE SCENE
# =========================================================================
def process_scene(slc_zip, idx, total, aoi_bounds):
    scene_name = os.path.basename(slc_zip).replace(".zip", "")
    short = scene_name[:55] + "..." if len(scene_name) > 55 else scene_name

    print(f"\n{'='*70}")
    print(f"  SCENE {idx}/{total}: {short}")
    print(f"{'='*70}")

    scene_dir = os.path.join(OUTPUT_DIR, scene_name)
    final_c2 = os.path.join(scene_dir, "C2")

    # Skip if already done
    if os.path.isfile(os.path.join(final_c2, "C11.tif")):
        ds = gdal.Open(os.path.join(final_c2, "C11.tif"))
        if ds and ds.GetGeoTransform() != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
            print(f"  ⏭️ Already done. Skipping.")
            ds = None; return final_c2
        if ds: ds = None
        shutil.rmtree(final_c2, ignore_errors=True)

    os.makedirs(scene_dir, exist_ok=True)

    # ── Process each sub-swath separately ──
    iw_geocoded_dirs = []

    for iw in ALL_SUBSWATHS:
        iw_inter = os.path.join(scene_dir, f"intermediate_{iw}")
        iw_dim = iw_inter + ".dim"
        iw_data = iw_inter + ".data"
        iw_geo = os.path.join(scene_dir, f"geocoded_{iw}")

        # Check if this sub-swath is already geocoded
        if os.path.isdir(iw_geo) and all(
            os.path.isfile(os.path.join(iw_geo, f"{b}.tif"))
            for b in ["C11", "C12_real", "C12_imag", "C22"]
        ):
            print(f"\n  ⏭️ {iw} already geocoded.")
            iw_geocoded_dirs.append(iw_geo)
            continue

        # Stage 1: SNAP GPT for this sub-swath
        if os.path.isfile(iw_dim) and os.path.isdir(iw_data):
            print(f"\n  ⏭️ {iw} SNAP done. Skipping to geocode.")
        else:
            print(f"\n  ── {iw}: SNAP GPT ──")
            print(f"  Read→Orbit→Split({iw})→Cal(complex)→Deburst→C2→ML→Speckle")
            graph_xml = create_graph_single_iw(slc_zip, iw_inter, iw)
            graph_file = os.path.join(scene_dir, f"graph_{iw}.xml")
            with open(graph_file, "w") as f: f.write(graph_xml)
            print(f"  Running SNAP GPT for {iw} (10-30 min)...")
            if not run_gpt(graph_file, f"{iw}"):
                print(f"  ⚠️ {iw} failed — will mosaic without it")
                continue
            if not os.path.isfile(iw_dim):
                print(f"  ❌ {iw} .dim not created"); continue

        # Stage 2: Geocode this sub-swath
        print(f"\n  ── {iw}: GDAL Geocode ──")
        if geocode_iw(iw_data, iw_dim, iw_geo):
            iw_geocoded_dirs.append(iw_geo)
            print(f"  ✅ {iw} geocoded → {iw_geo}")
        else:
            print(f"  ❌ {iw} geocoding failed")

    if not iw_geocoded_dirs:
        print(f"\n  ❌ No sub-swaths successfully processed!")
        return None

    # ── Stage 3: Mosaic all sub-swaths ──
    print(f"\n  ══ MOSAIC: {len(iw_geocoded_dirs)} sub-swaths → single GeoTIFF ══")
    result = mosaic_subswaths(iw_geocoded_dirs, final_c2, aoi_bounds)

    if result:
        print(f"\n  ✅ Scene {idx} complete! → {result}")
        # Cleanup intermediate sub-swath files
        for iw in ALL_SUBSWATHS:
            for sub in [f"intermediate_{iw}.dim", f"intermediate_{iw}.data",
                       f"geocoded_{iw}", f"graph_{iw}.xml"]:
                p = os.path.join(scene_dir, sub)
                if os.path.isdir(p): shutil.rmtree(p, ignore_errors=True)
                elif os.path.isfile(p): os.remove(p)
        print(f"  🧹 Cleaned intermediate files")

    return result


# =========================================================================
# MAIN
# =========================================================================
def main():
    print("=" * 70)
    print("  STEP 2: S1 SLC → C2 (v9 — 3×IW separate + GDAL mosaic)")
    print("  Each sub-swath processed separately, then mosaiced")
    print("=" * 70)

    slc_zips = sorted(glob.glob(os.path.join(DOWNLOADS_DIR, "S1*.zip")))
    if not slc_zips:
        slc_zips = sorted(glob.glob(os.path.join(DOWNLOADS_DIR, "*.zip")))
    if not slc_zips:
        print(f"\n❌ No zip files in: {DOWNLOADS_DIR}"); sys.exit(1)

    print(f"\n{len(slc_zips)} scene(s):")
    for i, z in enumerate(slc_zips, 1):
        print(f"  {i}. {os.path.basename(z)} ({os.path.getsize(z)/(1024**2):.0f} MB)")

    if not os.path.isfile(AOI_SHP):
        print(f"\n❌ AOI not found: {AOI_SHP}"); sys.exit(1)

    aoi = get_aoi_bounds(AOI_SHP)
    print(f"\nAOI: {aoi[0]:.4f},{aoi[1]:.4f} → {aoi[2]:.4f},{aoi[3]:.4f}")

    results = {}
    for i, z in enumerate(slc_zips, 1):
        r = process_scene(z, i, len(slc_zips), aoi)
        results[os.path.basename(z)] = r
        if r is None:
            print(f"\n  ⚠️ Scene {i} failed. Fix, re-run (done scenes auto-skip).")
            break

    print("\n" + "=" * 70 + "\n  SUMMARY\n" + "=" * 70)
    ok = sum(1 for v in results.values() if v)
    for s, p in results.items():
        print(f"  {'✅' if p else '❌'} {s[:60]}")
        if p: print(f"     → {p}")
    print(f"\n  {ok}/{len(results)} done.")
    if ok > 0: print("\n  NEXT: python step3_decompose.py")


if __name__ == "__main__":
    main()




