#!/usr/bin/env python3
"""
Move generated files from workspace root into pipeline data/results, re-run analysis,
and write RESULTS.md in results/ summarizing outputs and embedding images.
"""
import os, shutil, subprocess, sys

script_dir = os.path.dirname(os.path.abspath(__file__))
root = script_dir
data_dir = os.path.join(script_dir, 'data')
results_dir = os.path.join(script_dir, 'results')
os.makedirs(data_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

# Original FITS (keep these in root)
original_fits = set([
    '20220217191741-313-RA.wcs.proc.fits',
    '20220220190400-656-RA.wcs.proc.fits',
    '20220221190037-229-RA.wcs.proc.fits',
    '20220225192324-727-RA.wcs.proc.fits'
])

# Move TXT and PNG files from root into data/ or results/
for fname in os.listdir(root):
    fpath = os.path.join(root, fname)
    if not os.path.isfile(fpath):
        continue
    if fname in original_fits:
        continue
    lower = fname.lower()
    try:
        if lower.endswith('.txt'):
            # move txt to data if it looks like astrometry/obs, else to results
            if 'astrometry' in lower or 'centroid' in lower or 'observ' in lower or (fname.startswith('2022') and '.lspr.out' in lower):
                shutil.move(fpath, os.path.join(data_dir, fname))
            else:
                shutil.move(fpath, os.path.join(results_dir, fname))
        elif lower.endswith('.png'):
            shutil.move(fpath, os.path.join(results_dir, fname))
    except Exception as e:
        print('move failed', fname, e)

print('Workspace root cleaned; moved generated files into data/ and results/')

# Re-run analysis to refresh results
print('Running analyze_comet_67p.py to refresh results...')
ret = subprocess.run([sys.executable, os.path.join(script_dir, 'analyze_comet_67p.py')], cwd=script_dir)
if ret.returncode != 0:
    print('analyze_comet_67p.py failed with code', ret.returncode)
else:
    print('Analysis finished successfully')

# Build RESULTS.md from representative summary and include images
summary_file = os.path.join(results_dir, '67P_results_summary.txt')
results_md = os.path.join(results_dir, 'RESULTS.md')
with open(results_md, 'w', encoding='utf-8') as md:
    md.write('# Comet 67P — Analysis Results\n\n')
    if os.path.exists(summary_file):
        md.write('## Orbital Elements (representative)\n\n')
        md.write('''```
''')
        md.write(open(summary_file).read())
        md.write('\n```\n\n')
    else:
        md.write('No summary found. Check results directory.\n\n')
    md.write('## Plots\n\n')
    for img in sorted(os.listdir(results_dir)):
        if img.lower().endswith('.png'):
            md.write(f'![{img}](./{img})\n\n')

print('Wrote', results_md)
print('Done.')
