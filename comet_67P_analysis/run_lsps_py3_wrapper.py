#!/usr/bin/env python3
"""
Python3 wrapper for LSPR routines so we can plate-solve using the existing code
without requiring a Python 2 interpreter. This imports LSPR.py functions and calls
them with the measured centroids.
"""
import os
import sys
import shutil
import importlib
import glob

# compute repository root relative to this script and add to sys.path so
# we can import LSPR_fix (or LSPR) without any absolute hard-coded paths
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, '..'))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

LSPR = importlib.import_module('LSPR_fix')

data_dir = os.path.join(script_dir, 'data')
os.makedirs(data_dir, exist_ok=True)

# Prefer FITS and centroid files inside the comet_67P_analysis folder; fall back to repo root
repo_root = os.path.abspath(os.path.join(script_dir, '..'))
fits_dir_local = script_dir
fits_dir_root = repo_root
centroid_file = os.path.join(data_dir, 'comet_67P_centroids.txt')
# look for centroid file in several candidate locations (data/, analysis/, comet67P/, repo root)
candidates = [
    centroid_file,
    os.path.join(script_dir, 'comet_67P_centroids.txt'),
    os.path.join(script_dir, 'comet67P', 'comet_67P_centroids.txt'),
    os.path.join(repo_root, 'comet_67P_centroids.txt'),
    os.path.join(repo_root, 'comet67P', 'comet_67P_centroids.txt'),
]
found = None
for c in candidates:
    if c and os.path.exists(c):
        found = c
        break
if found and found != centroid_file:
    try:
        shutil.copy(found, centroid_file)
    except Exception:
        pass

# friendly early check: if centroids are not present, print guidance and exit
if not os.path.exists(centroid_file):
    print('No centroid file found at', centroid_file)
    print('Run the centroid extractor first, for example:')
    print('  python centroid_batch.py --fits-dir .\\comet67P --out comet_67P_centroids.txt')
    print('Then move the produced file into data/:')
    print('  Move-Item .\\comet_67P_centroids.txt .\\data\\comet_67P_centroids.txt')
    sys.exit(1)

out_agg = os.path.join(data_dir, 'comet_67P_astrometry.txt')

# Find LSPR star lists anywhere under the repository (prefer files inside analysis folder)
def find_starlist(name):
    # search inside analysis folder first, then repo tree
    candidates = [
        os.path.join(script_dir, name),
        os.path.join(script_dir, 'Method-of-Gauss-Orbital-Determination', name),
        os.path.join(repo_root, name),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    # try recursive glob search under repo_root
    hits = glob.glob(os.path.join(repo_root, '**', name), recursive=True)
    if hits:
        return hits[0]
    return None

star1 = find_starlist('LSPRstars1.txt')
star2 = find_starlist('LSPRstars2.txt')

star_map = {
    '20220217191741-313-RA.wcs.proc.fits': star1,
    '20220220190400-656-RA.wcs.proc.fits': star2,
    '20220221190037-229-RA.wcs.proc.fits': star1,
    '20220225192324-727-RA.wcs.proc.fits': star2,
}

results = []
for L in open(centroid_file).read().splitlines():
    if not L.strip() or L.startswith('#'):
        continue
    parts = [p.strip() for p in L.split(',')]
    fname = parts[0]
    x = float(parts[1]); y = float(parts[2])
    starfile = star_map.get(fname)
    if not starfile or not os.path.exists(starfile):
        print('Skipping', fname, 'no star list')
        continue
    outname = os.path.join(data_dir, fname + '.lspr.out.txt')
    print('Plate-solving', fname)
    try:
        # LSPR.run expects strings; call run() which internally writes the output file
        LSPR.run(starfile, x, y, outname)
        txt = open(outname).read()
        ra_line = [l for l in txt.splitlines() if l.strip().startswith('RA=')][0]
        dec_line = [l for l in txt.splitlines() if l.strip().startswith('Dec=')][0]
        ra = ra_line.split('=')[1].strip(); dec = dec_line.split('=')[1].strip()
        results.append((fname, ra, dec, outname))
    except Exception as e:
        print('  ERROR running LSPR on', fname, e)

with open(out_agg, 'w') as fo:
    fo.write('# filename, RA(h:m:s), DEC(±d:m:s), lspr_outfile\n')
    for r in results:
        fo.write('{},{},{},{}\n'.format(*r))

print('Wrote', out_agg)
print('Data files are stored in', data_dir)
