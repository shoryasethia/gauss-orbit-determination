#!/usr/bin/env python3
"""Compute O-C residuals (arcseconds) for the four observations using the representative
Gauss solution saved in results/67P_results_summary.txt.

This is a quick diagnostic: it propagates the representative r2/r2dot (assumed heliocentric
position and velocity in AU and AU/day) to the observation epochs, converts to RA/DEC
as seen from Earth, and compares to the observed RA/DEC in data/comet_67P_astrometry.txt.

This script uses astropy for coordinate transforms and assumes small light-time errors
(basic correction applied using instantaneous positions).
"""
import os
import numpy as np
from math import radians, degrees
from astropy.time import Time
from astropy.coordinates import SkyCoord, CartesianRepresentation, CartesianDifferential
from astropy import units as u
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, 'data')
results_dir = os.path.join(script_dir, 'results')

# Read representative solution (r2 and r2dot) from results_summary
summary_file = os.path.join(results_dir, '67P_results_summary.txt')
astrometry_file = os.path.join(data_dir, 'comet_67P_astrometry.txt')

if not os.path.exists(summary_file) or not os.path.exists(astrometry_file):
    print('Missing required files:', summary_file, astrometry_file)
    raise SystemExit(1)

# parse r2 and r2dot (simple text parsing)
with open(summary_file, 'r', encoding='utf-8') as fh:
    txt = fh.read()

import re
m_r = re.search(r'Position vector \(AU\): \[([^\]]+)\]', txt)
m_v = re.search(r'Velocity vector \(AU/day\): \[([^\]]+)\]', txt)
if not m_r or not m_v:
    print('Could not parse r2/r2dot from', summary_file)
    raise SystemExit(1)

r2 = np.array([float(x) for x in m_r.group(1).split(',')])
r2dot = np.array([float(x) for x in m_v.group(1).split(',')])
print('Representative r2 (AU):', r2)
print('Representative r2dot (AU/day):', r2dot)

# Read astrometry: filename, RA(h:m:s), DEC(±d:m:s), ...
obs = []
with open(astrometry_file, 'r', encoding='utf-8', errors='replace') as fh:
    for line in fh:
        line=line.strip()
        if not line or line.startswith('#'):
            continue
        parts=[p.strip() for p in line.split(',')]
        fname=parts[0]
        ra=parts[1]
        dec=parts[2]
        obs.append((fname,ra,dec))

# Need observation times: parse from filenames
def parse_time_from_fname(fname):
    # filename begins with YYYYMMDDhhmmss-...
    base=os.path.basename(fname)
    dt=base.split('-')[0]
    tstr=f'{dt[0:4]}-{dt[4:6]}-{dt[6:8]} {dt[8:10]}:{dt[10:12]}:{dt[12:14]}.000'
    return Time(tstr, format='iso', scale='utc')

# We'll compare heliocentric positions propagated linearly (small dt) and compute topocentric
# vector by subtracting Earth's barycentric position at the same time (approx). This is crude
# but sufficient for a quick O-C diagnostic.

with solar_system_ephemeris.set('builtin'):
    for fname, ra_obs, dec_obs in obs:
        t = parse_time_from_fname(fname)
        # Earth's barycentric position (AU)
        earth = get_body_barycentric_posvel('earth', t)[0]
        earth_xyz = np.array([earth.x.value, earth.y.value, earth.z.value])
        # Propagate comet r2 from representative epoch (we don't have epoch stored; assume r2 corresponds to obs[0] time)
        # For a quick check, propagate linearly from r2 using r2dot*(dt days)
        # Use t0 = first observation time
    
# assume representative r2 corresponds approximately to the first observation time in obs
t0 = parse_time_from_fname(obs[0][0])
print('Assuming representative solution epoch ~', t0.iso)

print('\nO - C residuals (approx, arcseconds):')
print('file, observed_RA, observed_DEC, computed_RA, computed_DEC, dRA*cos(Dec) ("), dDec (")')

with solar_system_ephemeris.set('builtin'):
    for fname, ra_obs, dec_obs in obs:
        t_obs = parse_time_from_fname(fname)
        dt_days = (t_obs - t0).to(u.day).value
        # linear propagation (r + v*dt)
        r_prop = r2 + r2dot * dt_days
        # vector from Earth to comet (approx): r_prop - earth_pos
        earth = get_body_barycentric_posvel('earth', t_obs)[0]
        earth_xyz = np.array([earth.x.value, earth.y.value, earth.z.value])
        topoc = r_prop - earth_xyz
        # convert to astropy SkyCoord in ICRS
        cart = CartesianRepresentation(topoc[0]*u.AU, topoc[1]*u.AU, topoc[2]*u.AU)
        coord = SkyCoord(cart, frame='icrs')
        # Ensure we have spherical coordinates for ra/dec
        ra_comp = coord.ra.to_string(unit=u.hour, sep=':', pad=True)
        dec_comp = coord.dec.to_string(unit=u.deg, sep=':', pad=True, alwayssign=True)
        # compute residuals in arcsec
        c_obs = SkyCoord(ra_obs + ' ' + dec_obs, unit=(u.hourangle, u.deg), frame='icrs')
        c_comp = SkyCoord(ra_comp + ' ' + dec_comp, unit=(u.hourangle, u.deg), frame='icrs')
        sep = c_obs.separation(c_comp).arcsecond
        # separate RA/Dec residuals (small-angle approx)
        dra = (c_obs.ra.radian - c_comp.ra.radian) * np.cos(c_comp.dec.radian) * 206264.806
        ddec = (c_obs.dec.radian - c_comp.dec.radian) * 206264.806
        print(f"{os.path.basename(fname)}, {ra_obs}, {dec_obs}, {ra_comp}, {dec_comp}, {dra:.3f}, {ddec:.3f}")

print('\nNote: this is a quick diagnostic (linear propagation). For accurate O-C you should')
print(' propagate the orbit properly from mean anomaly and include light-time correction and aberration.\n')
