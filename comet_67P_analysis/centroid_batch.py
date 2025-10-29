#!/usr/bin/env python3
"""
Batch centroid extraction for comet images.

Usage:
  python centroid_batch.py --fits-dir ".." --out comet_67P_centroids.txt

This script tries to locate the comet in each FITS by doing a simple
flux-weighted centroid inside a small box. It prefers to use WCS/CRPIX
to center the search box when available, otherwise uses image center.

Outputs a text file with lines:
  filename, x_pix, y_pix, flux_sum, sigma_x, sigma_y

Notes: This is a simple automatic tool intended to produce starting
centroids for plate solving. For highest accuracy run the interactive
`centroidprogram.py` from the original repo and plate-solve with LSPR.
"""

import argparse
import os
import math
from astropy.io import fits
import numpy as np


def flux_weighted_centroid(data, x0, y0, halfbox=25):
    ny, nx = data.shape
    x0 = int(round(x0))
    y0 = int(round(y0))
    x1 = max(0, x0 - halfbox)
    x2 = min(nx, x0 + halfbox)
    y1 = max(0, y0 - halfbox)
    y2 = min(ny, y0 + halfbox)
    cut = data[y1:y2, x1:x2].astype(float)
    # subtract local median background
    med = np.median(cut)
    cut_sub = cut - med
    cut_sub[cut_sub < 0] = 0.0
    total = cut_sub.sum()
    if total <= 0:
        return None
    ys, xs = np.indices(cut_sub.shape)
    cx = (cut_sub * xs).sum() / total + x1
    cy = (cut_sub * ys).sum() / total + y1
    # estimate uncertainties as sqrt(variance / flux)
    varx = (cut_sub * (xs - (cx - x1))**2).sum() / total
    vary = (cut_sub * (ys - (cy - y1))**2).sum() / total
    sigma_x = math.sqrt(varx) if varx > 0 else 0.0
    sigma_y = math.sqrt(vary) if vary > 0 else 0.0
    return cx, cy, total, sigma_x, sigma_y


def find_seed_from_header(hdr, nx, ny):
    # try CRPIX values
    crpix1 = hdr.get('CRPIX1')
    crpix2 = hdr.get('CRPIX2')
    if crpix1 is not None and crpix2 is not None:
        return float(crpix1 - 1), float(crpix2 - 1)  # FITS to 0-index
    # fallback: image center
    return nx / 2.0, ny / 2.0


def process_file(path, halfbox=25):
    try:
        with fits.open(path, memmap=False) as hdul:
            hdr = hdul[0].header
            data = hdul[0].data
            if data is None:
                # try next extension
                for hdu in hdul:
                    if hdu.data is not None:
                        data = hdu.data
                        hdr = hdu.header
                        break
            if data is None:
                return (path, None, 'no image data')
            # ensure 2D
            if data.ndim > 2:
                # take first extension or average across planes
                data = data[0]
            ny, nx = data.shape
            x0, y0 = find_seed_from_header(hdr, nx, ny)
            result = flux_weighted_centroid(data, x0, y0, halfbox=halfbox)
            if result is None:
                # try bigger box
                result = flux_weighted_centroid(data, x0, y0, halfbox=halfbox*2)
                if result is None:
                    return (path, None, 'centroid failed')
            cx, cy, total, sx, sy = result
            return (path, (cx, cy, total, sx, sy), None)
    except Exception as e:
        return (path, None, str(e))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fits-dir', default='..', help='Directory containing FITS files')
    p.add_argument('--out', default='comet_67P_centroids.txt')
    p.add_argument('--halfbox', type=int, default=25)
    args = p.parse_args()

    fits_dir = args.fits_dir
    files = [f for f in os.listdir(fits_dir) if f.lower().endswith('.fits') or f.lower().endswith('.fit')]
    files = sorted(files)
    if not files:
        print('No FITS files found in', fits_dir)
        return
    outlines = []
    for fn in files:
        path = os.path.join(fits_dir, fn)
        print('Processing', fn)
        path, result, err = process_file(path, halfbox=args.halfbox)
        if err:
            print('  ERROR:', err)
            outlines.append('{},{},{},{},{}'.format(fn, 'nan', 'nan', '0.0', err))
        else:
            cx, cy, total, sx, sy = result
            outlines.append('{},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}'.format(fn, cx, cy, total, sx, sy))

    with open(os.path.join(fits_dir, args.out), 'w') as fh:
        fh.write('# filename, x_pix, y_pix, flux_sum, sigma_x, sigma_y\n')
        for L in outlines:
            fh.write(L + '\n')

    print('Wrote', os.path.join(fits_dir, args.out))


if __name__ == '__main__':
    main()
