#!/usr/bin/env python3
"""Accurate O-C residuals using poliastro two-body propagation.

Requires: poliastro, astropy, numpy

This script:
- parses the representative result in results/67P_results_summary.txt to get r2,r2dot and the combo used
- builds a poliastro Orbit from r2/r2dot at the epoch corresponding to the middle observation of the combo
- propagates to each observation time and computes topocentric RA/DEC by subtracting Earth's heliocentric position
- prints O-C residuals (arcseconds) and an RMS summary
"""
import os
import re
import numpy as np
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric

from poliastro.twobody import Orbit
from poliastro.bodies import Sun

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, 'data')
results_dir = os.path.join(script_dir, 'results')

summary_file = os.path.join(results_dir, '67P_results_summary.txt')
astrometry_file = os.path.join(data_dir, 'comet_67P_astrometry.txt')

if not os.path.exists(summary_file) or not os.path.exists(astrometry_file):
    print('Missing required files; run pipeline first')
    raise SystemExit(1)

# parse r2, r2dot and combo
with open(summary_file, 'r', encoding='utf-8') as fh:
    txt = fh.read()

combo_m = re.search(r'Representative result from combo:\s*\(([^\)]+)\)', txt)
pos_m = re.search(r'Position vector \(AU\): \[([^\]]+)\]', txt)
vel_m = re.search(r'Velocity vector \(AU/day\): \[([^\]]+)\]', txt)
if not combo_m or not pos_m or not vel_m:
    print('Could not parse summary file; expected combo, position and velocity')
    raise SystemExit(1)

combo = tuple(int(x.strip()) for x in combo_m.group(1).split(','))
r2 = np.array([float(x) for x in pos_m.group(1).split(',')])
r2dot = np.array([float(x) for x in vel_m.group(1).split(',')])

print('Using representative combo', combo)
print('r2 (AU):', r2)
print('r2dot (AU/day):', r2dot)

# read astrometry to get times and observed RA/DEC
obs = []
with open(astrometry_file, 'r', encoding='utf-8', errors='replace') as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = [p.strip() for p in line.split(',')]
        fname = parts[0]
        ra = parts[1]
        dec = parts[2]
        obs.append((fname, ra, dec))

def parse_time(fname):
    base = os.path.basename(fname)
    dt = base.split('-')[0]
    tstr = f'{dt[0:4]}-{dt[4:6]}-{dt[6:8]} {dt[8:10]}:{dt[10:12]}:{dt[12:14]}.000'
    return Time(tstr, format='iso', scale='utc')

# Determine epoch: r2 from Gauss is the vector at the middle observation (t2)
mid_index = combo[1]
epoch_time = parse_time(obs[mid_index][0])

# Build poliastro Orbit from vectors: convert AU->km and AU/day -> km/s
AU_km = 149597870.7
r_km = (r2 * AU_km) * u.km
v_kms = (r2dot * AU_km / 86400.0) * (u.km / u.s)

orbit = Orbit.from_vectors(Sun, r_km, v_kms, epoch=epoch_time)
print('Built Orbit at epoch', orbit.epoch.iso)

# propagate and compute O-C
residuals = []
with solar_system_ephemeris.set('builtin'):
    for i, (fname, ra_obs, dec_obs) in enumerate(obs):
        t_obs = parse_time(fname)
        # propagate orbit to t_obs
        tof = (t_obs - orbit.epoch).to(u.s)
        r_prop, v_prop = orbit.propagate(tof).rv()
        # r_prop is heliocentric vector in km
        r_comp_au = (r_prop.to(u.km).value / AU_km)
        # get Earth barycentric and Sun barycentric positions to form Earth heliocentric
        earth_bary = get_body_barycentric('earth', t_obs)
        sun_bary = get_body_barycentric('sun', t_obs)
        earth_helio = np.array([earth_bary.x.value - sun_bary.x.value,
                               earth_bary.y.value - sun_bary.y.value,
                               earth_bary.z.value - sun_bary.z.value])
        # topocentric vector (AU)
        topoc = r_comp_au - earth_helio
        # convert to RA/DEC using astropy SkyCoord
        coord = SkyCoord(x=topoc[0]*u.AU, y=topoc[1]*u.AU, z=topoc[2]*u.AU, frame='icrs', representation_type='cartesian')
        ra_comp = coord.icrs.ra.to_string(unit=u.hour, sep=':', pad=True)
        dec_comp = coord.icrs.dec.to_string(unit=u.deg, sep=':', pad=True, alwayssign=True)
        c_obs = SkyCoord(ra_obs + ' ' + dec_obs, unit=(u.hourangle, u.deg), frame='icrs')
        c_comp = SkyCoord(ra_comp + ' ' + dec_comp, unit=(u.hourangle, u.deg), frame='icrs')
        sep_arcsec = c_obs.separation(c_comp).arcsecond
        # RA/Dec residual components (approx)
        dra = (c_obs.ra.radian - c_comp.ra.radian) * np.cos(c_comp.dec.radian) * 206264.806
        ddec = (c_obs.dec.radian - c_comp.dec.radian) * 206264.806
        residuals.append((fname, ra_obs, dec_obs, ra_comp, dec_comp, dra, ddec, sep_arcsec))

print('\nFile, obs_RA, obs_DEC, comp_RA, comp_DEC, dRA*cos(Dec) ("), dDec ("), sep (\")')
for r in residuals:
    print(f"{os.path.basename(r[0])}, {r[1]}, {r[2]}, {r[3]}, {r[4]}, {r[5]:.3f}, {r[6]:.3f}, {r[7]:.3f}")

seps = np.array([r[7] for r in residuals])
print(f"\nRMS separation = {np.sqrt(np.mean(seps**2)):.3f} arcsec")

print('\nNote: this propagation uses a two-body model (heliocentric) and converts to apparent RA/DEC by subtracting Earth heliocentric position. For highest accuracy include light-time correction and observer location (topocenter).')
