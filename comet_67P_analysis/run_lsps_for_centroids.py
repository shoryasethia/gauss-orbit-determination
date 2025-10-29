#!/usr/bin/env python3
"""
Wrapper to run the repository LSPR_fix (Python3-adapted) plate-solver on measured centroids.

This version is safe to run from inside `comet_67P_analysis/`. It looks for
`comet_67P_centroids.txt` inside `comet_67P_analysis/data/` and writes
`comet_67P_astrometry.txt` back into `comet_67P_analysis/data/`.
"""
import os
import sys
import importlib

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, '..'))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

# prefer reading centroids from data/ inside the analysis folder
fits_dir = script_dir
data_dir = os.path.join(script_dir, 'data')
os.makedirs(data_dir, exist_ok=True)
centroid_file = os.path.join(data_dir, 'comet_67P_centroids.txt')
out_agg = os.path.join(data_dir, 'comet_67P_astrometry.txt')

LSPR = importlib.import_module('LSPR_fix')

# mapping: filename -> star-list file to use (star lists are in the repo root)
star_map = {
    '20220217191741-313-RA.wcs.proc.fits': os.path.join(repo_root, 'LSPRstars1.txt'),
    '20220220190400-656-RA.wcs.proc.fits': os.path.join(repo_root, 'LSPRstars2.txt'),
    '20220221190037-229-RA.wcs.proc.fits': os.path.join(repo_root, 'LSPRstars1.txt'),
    '20220225192324-727-RA.wcs.proc.fits': os.path.join(repo_root, 'LSPRstars2.txt'),
}

lines = []
if os.path.exists(centroid_file):
    lines = open(centroid_file).read().splitlines()
else:
    print('No centroid file found at', centroid_file)
    sys.exit(1)

results = []
for L in lines:
    if L.startswith('#') or not L.strip():
        continue
    parts = [p.strip() for p in L.split(',')]
    fname = parts[0]
    try:
        x = float(parts[1]); y = float(parts[2])
    except Exception:
        print('Skipping malformed centroid line:', L)
        continue
    starfile = star_map.get(fname)
    if not starfile or not os.path.exists(starfile):
        print('Skipping', fname, 'no star list at', starfile)
        continue
    outname = os.path.join(data_dir, fname + '.lspr.out.txt')
    print('Plate-solving', fname)
    try:
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
