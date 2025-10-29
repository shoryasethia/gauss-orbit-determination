#!/usr/bin/env python3
"""Run the full comet 67P pipeline (centroids, plate-solve, Gauss, jackknife/mag, report).
This is the single, obvious entrypoint to reproduce the full analysis end-to-end.
"""
import os
import shutil
import subprocess
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(script_dir, '..'))
data_dir = os.path.join(script_dir, 'data')
results_dir = os.path.join(script_dir, 'results')
os.makedirs(data_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

print('Workspace root cleaned; moved generated files into data/ and results/')
# Move any previously-generated .txt and .png files found in the root into data/ or results/
for f in os.listdir(root):
    if f.lower().endswith('.txt') or f.lower().endswith('.png'):
        src = os.path.join(root, f)
        if os.path.isfile(src):
            if f.lower().endswith('.png'):
                dst = os.path.join(results_dir, f)
            else:
                dst = os.path.join(data_dir, f)
            try:
                shutil.move(src, dst)
            except Exception:
                pass

# Step 1: ensure centroids exist
centroid_parent = os.path.join(root, 'comet_67P_centroids.txt')
centroid_parent_local = os.path.join(script_dir, 'comet_67P_centroids.txt')
centroid_data = os.path.join(data_dir, 'comet_67P_centroids.txt')
if os.path.exists(centroid_parent) and not os.path.exists(centroid_data):
    shutil.copy(centroid_parent, centroid_data)
elif os.path.exists(centroid_parent_local) and not os.path.exists(centroid_data):
    shutil.copy(centroid_parent_local, centroid_data)

# Step 2: run centroid extractor if needed
if not os.path.exists(centroid_data):
    print('Running centroid extractor...')
    # prefer a 'comet67P' subfolder inside comet_67P_analysis, then the analysis folder itself, then repo root
    fits_src = None
    # look for a folder named comet67P (case-insensitive)
    for d in os.listdir(script_dir):
        p = os.path.join(script_dir, d)
        if os.path.isdir(p) and d.lower() == 'comet67p':
            fits_src = p
            break
    if fits_src is None:
        # use analysis folder if it contains FITS, else fallback to repo root
        has_local_fits = any(f.lower().endswith('.fits') or f.lower().endswith('.fit') for f in os.listdir(script_dir))
        if has_local_fits:
            fits_src = script_dir
        else:
            fits_src = root
    subprocess.check_call([sys.executable, os.path.join(script_dir, 'centroid_batch.py'), '--fits-dir', fits_src, '--out', 'comet_67P_centroids.txt'])
    # move produced file
    # produced file is written into the fits_src location
    produced = os.path.join(fits_src, 'comet_67P_centroids.txt')
    if os.path.exists(produced):
        shutil.move(produced, centroid_data)

# Step 3: run plate-solve wrapper
print('Running plate-solver wrapper...')
subprocess.check_call([sys.executable, os.path.join(script_dir, 'run_lsps_py3_wrapper.py')])

# Step 4: run analysis
print('Running analyze_comet_67p.py to refresh results...')
subprocess.check_call([sys.executable, os.path.join(script_dir, 'analyze_comet_67p.py')])

# Step 5: run jackknife + magnitude
print('Running jackknife and magnitude summary...')
subprocess.check_call([sys.executable, os.path.join(script_dir, 'run_jackknife_and_magnitude.py')])

print('Wrote', os.path.join(results_dir, 'RESULTS.md'))
print('Done.')
