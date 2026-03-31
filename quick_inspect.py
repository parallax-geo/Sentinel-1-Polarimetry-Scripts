#!/usr/bin/env python3
"""Inspect the DIMAP tie-point grids and folder structure."""
import os, glob

base = os.path.expanduser("~/Desktop/Muddasir/outputs/preprocessed")
scene = glob.glob(os.path.join(base, "S1A*"))[0]
data_dir = os.path.join(scene, "intermediate_c2.data")
dim_file = os.path.join(scene, "intermediate_c2.dim")

# 1. Show ALL folder contents including subdirs
print("=== .data folder contents ===")
for root, dirs, files in os.walk(data_dir):
    level = root.replace(data_dir, '').count(os.sep)
    indent = '  ' * level
    print(f"{indent}📁 {os.path.basename(root)}/")
    for f in sorted(files):
        sz = os.path.getsize(os.path.join(root, f))
        if sz > 1024*1024:
            print(f"{indent}  📄 {f} ({sz/(1024*1024):.1f} MB)")
        else:
            print(f"{indent}  📄 {f} ({sz/1024:.0f} KB)")

# 2. Extract tie-point grid info from .dim XML
print(f"\n=== Tie-Point Grid Info from .dim ===")
with open(dim_file, 'r') as f:
    content = f.read()

# Find tie point grid sections
import re
# Look for Tie_Point_Grid_Info blocks
tp_blocks = re.findall(r'<Tie_Point_Grid_Info>.*?</Tie_Point_Grid_Info>', content, re.DOTALL)
print(f"Found {len(tp_blocks)} tie-point grid block(s)")
for block in tp_blocks:
    print(block[:500])
    print("---")

# Also look for any lat/lon references
for tag in ['NCOLS', 'NROWS', 'OFFSET_X', 'OFFSET_Y', 'STEP_X', 'STEP_Y', 
            'TIE_POINT_GRID_NAME', 'latitude', 'longitude']:
    matches = re.findall(f'<{tag}>(.*?)</{tag}>', content, re.IGNORECASE)
    if matches:
        print(f"  {tag}: {matches[:5]}")

# 3. Try opening .img files with GDAL
print(f"\n=== GDAL open test ===")
from osgeo import gdal
gdal.UseExceptions()

for f in sorted(os.listdir(data_dir)):
    if f.endswith('.img'):
        fp = os.path.join(data_dir, f)
        try:
            ds = gdal.Open(fp)
            print(f"  {f}: {ds.RasterXSize}x{ds.RasterYSize}, bands={ds.RasterCount}")
            print(f"    GT: {ds.GetGeoTransform()}")
            print(f"    GCPs: {len(ds.GetGCPs())}")
            ds = None
        except Exception as e:
            print(f"  {f}: ERROR → {e}")