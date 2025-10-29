#!/usr/bin/env python3
"""Compute jackknife uncertainties from Gauss combo outputs and simple magnitudes from centroid fluxes.

Produces:
 - results/67P_uncertainties.txt
 - results/67P_magnitude.txt
 - appends summaries to results/RESULTS.md
"""
import os
import re
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, 'data')
results_dir = os.path.join(script_dir, 'results')
os.makedirs(results_dir, exist_ok=True)

combos_file = os.path.join(results_dir, '67P_results_combos.txt')
centroids_file = os.path.join(data_dir, 'comet_67P_centroids.txt')
results_out = os.path.join(results_dir, '67P_uncertainties.txt')
mag_out = os.path.join(results_dir, '67P_magnitude.txt')
results_md = os.path.join(results_dir, 'RESULTS.md')

def parse_combos(path):
    """Parse the combos file and return list of element tuples (a,e,Omega,i,omega,M)."""
    elems = []
    if not os.path.exists(path):
        return elems
    with open(path, 'r', encoding='utf-8') as fh:
        text = fh.read()
    # find lines like: a=..., e=..., Omega=..., i=..., w=..., M=...
    pattern = re.compile(r'a=([0-9Ee+\-.]+) AU, e=([0-9Ee+\-.]+), Omega=([0-9Ee+\-.]+) deg, i=([0-9Ee+\-.]+) deg, w=([0-9Ee+\-.]+) deg, M=([0-9Ee+\-.]+) deg')
    for m in pattern.finditer(text):
        a = float(m.group(1))
        e = float(m.group(2))
        Omega = float(m.group(3))
        i = float(m.group(4))
        w = float(m.group(5))
        M = float(m.group(6))
        elems.append((a,e,Omega,i,w,M))
    return elems

def jackknife_stats(values):
    """Compute jackknife mean and standard error for a 1D list.
    Returns (mean, std_error)
    """
    vals = np.array(values)
    n = len(vals)
    if n == 0:
        return (np.nan, np.nan)
    mean = np.mean(vals)
    # jackknife pseudo-values
    jack_means = np.array([np.mean(np.delete(vals, i)) for i in range(n)])
    jk_mean = np.mean(jack_means)
    # standard error estimate
    se = np.sqrt((n-1)/n * np.sum((jack_means - jk_mean)**2))
    return mean, se

def compute_jackknife(elems_list):
    # elems_list: list of tuples (a,e,Omega,i,w,M)
    if not elems_list:
        return None
    arr = np.array(elems_list)
    stats = {}
    labels = ['a_AU','e','Omega_deg','i_deg','w_deg','M_deg']
    for idx,label in enumerate(labels):
        vals = arr[:, idx]
        stats[label] = jackknife_stats(vals)
    return stats

def parse_centroids(path):
    # expected CSV per-line: filename, x_pix, y_pix, flux_sum, sigma_x, sigma_y
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            parts = [p.strip() for p in line.split(',') if p.strip()]
            if len(parts) < 4:
                continue
            fname = parts[0]
            try:
                x = float(parts[1]); y = float(parts[2]); flux = float(parts[3])
            except Exception:
                continue
            rows.append((fname, x, y, flux))
    return rows

def compute_instrumental_mags(centroids):
    # convert flux -> instrumental mag using -2.5 log10(flux)
    res = []
    for fname,x,y,flux in centroids:
        if flux <= 0:
            mag = np.nan
        else:
            mag = -2.5 * np.log10(flux)
        res.append((fname, x, y, flux, mag))
    return res

def write_uncertainties(stats, outpath):
    with open(outpath, 'w', encoding='utf-8') as fh:
        if stats is None:
            fh.write('# No Gauss combo elements found; cannot compute jackknife uncertainties.\n')
            return
        fh.write('# Jackknife uncertainties computed from Gauss 3-observation combinations\n')
        for k,v in stats.items():
            mean,se = v
            fh.write(f'{k}: mean={mean:.6f}, jackknife_std_error={se:.6f}\n')

def write_magnitudes(mags, outpath):
    with open(outpath, 'w', encoding='utf-8') as fh:
        fh.write('# filename, x_pix, y_pix, flux_sum, instrumental_mag\n')
        for fname,x,y,flux,mag in mags:
            fh.write(f'{fname},{x:.3f},{y:.3f},{flux:.6f},{mag if not np.isnan(mag) else "nan"}\n')

def append_results_md(stats, mags, mdpath):
    lines = []
    lines.append('\n## Jackknife uncertainties and instrumental magnitudes\n')
    if stats is None:
        lines.append('\nNo Gauss combo results were available to compute jackknife uncertainties.\n')
    else:
        lines.append('\nJackknife summary (mean ± jackknife standard error):\n\n')
        for k,v in stats.items():
            mean,se = v
            lines.append(f'- {k}: {mean:.6f} ± {se:.6f}\n')
    if mags:
        lines.append('\nInstrumental magnitudes derived from centroid flux sums (no zeropoint applied):\n\n')
        for fname,x,y,flux,mag in mags:
            lines.append(f'- {os.path.basename(fname)}: flux={flux:.6f} -> inst_mag={mag if not np.isnan(mag) else "nan"}\n')
    with open(mdpath, 'a', encoding='utf-8') as fh:
        fh.writelines(lines)

def main():
    elems = parse_combos(combos_file)
    stats = compute_jackknife(elems) if elems else None
    write_uncertainties(stats, results_out)
    centroids = parse_centroids(centroids_file)
    mags = compute_instrumental_mags(centroids)
    write_magnitudes(mags, mag_out)
    append_results_md(stats, mags, results_md)
    print('✓ Jackknife and magnitude outputs written:')
    print('  -', results_out)
    print('  -', mag_out)
    print('  - Appended summary to', results_md)

if __name__ == '__main__':
    main()
