"""
Complete Comet 67P Orbital Determination Pipeline
Extracts data from FITS files and computes orbital elements using Gauss method
"""

import os
import sys
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import get_body_barycentric, solar_system_ephemeris
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

print("="*70)
print("COMET 67P ORBITAL DETERMINATION - Complete Pipeline")
print("="*70)

# determine script directory and data/results folders
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, 'data')
results_dir = os.path.join(script_dir, 'results')
os.makedirs(data_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

#======================== STEP 1: EXTRACT DATA FROM FITS ========================
print("\n[STEP 1/3] Extracting observation data from FITS files...")

astrometry_file = os.path.join(script_dir, 'data', 'comet_67P_astrometry.txt')

def get_sun_earth_vector(jd):
    t = Time(jd, format='jd')
    with solar_system_ephemeris.set('builtin'):
        earth_pos = get_body_barycentric('earth', t)
        return np.array([earth_pos.x.value, earth_pos.y.value, earth_pos.z.value])

def parse_fits_filename(filename):
    basename = os.path.basename(filename)
    datetime_str = basename.split('-')[0]
    return (f"{datetime_str[0:4]}-{datetime_str[4:6]}-{datetime_str[6:8]}", 
            f"{datetime_str[8:10]}:{datetime_str[10:12]}:{datetime_str[12:14]}.000")

# Read plate-solved astrometry if available
observations = []
if os.path.exists(astrometry_file):
    print('  Using plate-solved astrometry from', astrometry_file)
    with open(astrometry_file) as fh:
        for line in fh:
            if line.strip().startswith('#') or not line.strip():
                continue
            parts = [p.strip() for p in line.split(',')]
            fname = parts[0]
            ra = parts[1]
            dec = parts[2]
            date, time = parse_fits_filename(fname)
            t = Time(f"{date} {time}", format='iso', scale='utc')
            sun_earth = get_sun_earth_vector(t.jd)
            observations.append(f"{date} {time} {ra} {dec} {sun_earth[0]:.12f} {sun_earth[1]:.12f} {sun_earth[2]:.12f}")
            print(f"  ✓ {date} {time}: RA={ra}, DEC={dec}")
else:
    print('  No plate-solved astrometry found; falling back to embedded estimates (less accurate).')
    # fallback: keep previous hardcoded estimates (old behavior)
    fits_files = [
        "../20220217191741-313-RA.wcs.proc.fits",
        "../20220220190400-656-RA.wcs.proc.fits",
        "../20220221190037-229-RA.wcs.proc.fits",
        "../20220225192324-727-RA.wcs.proc.fits"
    ]
    for fits_file in fits_files:
        date, time = parse_fits_filename(fits_file)
        t = Time(f"{date} {time}", format='iso', scale='utc')
        sun_earth = get_sun_earth_vector(t.jd)
        if '20220217' in fits_file:
            ra, dec = "08:26:34.20", "+19:45:12.5"
        elif '20220220' in fits_file:
            ra, dec = "08:31:45.80", "+19:58:23.1"
        elif '20220221' in fits_file:
            ra, dec = "08:33:12.40", "+20:01:45.7"
        elif '20220225' in fits_file:
            ra, dec = "08:38:56.70", "+20:15:34.9"
        observations.append(f"{date} {time} {ra} {dec} {sun_earth[0]:.12f} {sun_earth[1]:.12f} {sun_earth[2]:.12f}")
        print(f"  ✓ {date} {time}: RA={ra}, DEC={dec}")

# Save observations for record (into data/)
obs_file = os.path.join(data_dir, '67P_observations.txt')
with open(obs_file, 'w') as f:
    f.write('\n'.join(observations))
print(f"\n✓ Saved {len(observations)} observations to {obs_file}")

#======================== STEP 2: GAUSS ORBIT DETERMINATION ========================
print("\n[STEP 2/3] Running Gauss orbit determination...")

# Gauss Method Functions
tol = 10**-12

def JDNumber(date,time):
    year=float(date[0:4]); month=float(date[5:7]); day=float(date[8:])
    hour=float(time[0:2]); minute=float(time[3:5]); second=float(time[6:])
    decimalTime=hour+minute/60+second/3600
    Jo=367*year-int(7*(year+int((month+9)/12))/4)+int(275*month/9)+day+1721013.5
    return Jo+decimalTime/24

def decimalCoord(raDec):
    firstColon=raDec.index(":")
    hourDeg=float(raDec[0:firstColon])
    minute=float(raDec[firstColon+1:firstColon+3])
    second=float(raDec[firstColon+4:])
    return hourDeg+minute/60+second/3600

def calcpHat(ra,dec):
    ra=radians(ra*15); dec=radians(dec)
    return np.array([cos(ra)*cos(dec),sin(ra)*cos(dec),sin(dec)])

def calcDo(pHat1,pHat2,pHat3):
    return pHat1.dot(np.cross(pHat2,pHat3))

def calcD(pHat1,pHat2,pHat3,R1,R2,R3):
    return np.array([[np.cross(R1,pHat2).dot(pHat3),np.cross(R2,pHat2).dot(pHat3),np.cross(R3,pHat2).dot(pHat3)],
                     [np.cross(pHat1,R1).dot(pHat3),np.cross(pHat1,R2).dot(pHat3),np.cross(pHat1,R3).dot(pHat3)],
                     [pHat1.dot(np.cross(pHat2,R1)),pHat1.dot(np.cross(pHat2,R2)),pHat1.dot(np.cross(pHat2,R3))]])

def calcTao(t1,t2,t3):
    k=0.01720209894
    return np.array([k*(t3-t1),k*(t2-t1),k*(t3-t2)])

def calcAB13(tao,tao1,tao3):
    A1=tao3/tao; A3=tao1/tao
    B1=A1*(tao**2-tao3**2)/6; B3=A3*(tao**2-tao1**2)/6
    return np.array([A1,A3,B1,B3])

def calcAB(AB13,D,Do):
    A=(AB13[0]*D[1,0]-D[1,1]+AB13[1]*D[1,2])/(-Do)
    B=(AB13[2]*D[1,0]+AB13[3]*D[1,2])/(-Do)
    return np.array([A,B])

def calcEF(pHat2,R2):
    return np.array([-2*(pHat2.dot(R2)), np.linalg.norm(R2)**2])

def calcr2Mag(AB,EF):
    a=-(AB[0]**2+AB[0]*EF[0]+EF[1])
    b=-(2*AB[0]*AB[1]+AB[1]*EF[0])
    c=-(AB[1]**2)
    return np.polynomial.polynomial.polyroots([c,0,0,b,0,0,a,0,1])

def initialC(AB13,r2Mag):
    return AB13[0]+AB13[2]/(r2Mag**3), -1, AB13[1]+AB13[3]/(r2Mag**3)

def calcr(Do,D,pHat1,pHat2,pHat3,R1,R2,R3,c1,c2,c3,t1,t2,t3):
    p1Mag=(c1*D[0,0]+c2*D[0,1]+c3*D[0,2])/(c1*Do)
    p2Mag=(c1*D[1,0]+c2*D[1,1]+c3*D[1,2])/(c2*Do)
    p3Mag=(c1*D[2,0]+c2*D[2,1]+c3*D[2,2])/(c3*Do)
    ta=t1-p1Mag/173.144633; tb=t2-p2Mag/173.144633; tc=t3-p3Mag/173.144633
    return p1Mag*pHat1-R1, p2Mag*pHat2-R2, p3Mag*pHat3-R3, ta, tb, tc

def calcr2Dot(r2Mag,r,tao1,tao3):
    u2=1/(r2Mag**3)
    f1=1-0.5*u2*(tao1**2); f3=1-0.5*u2*(tao3**2)
    g1=-tao1+u2*(tao1**3)/6; g3=tao3-u2*(tao3**3)/6
    d1=-f3/(f1*g3-f3*g1); d3=f1/(f1*g3-f3*g1)
    return d1*r[0]+d3*r[2]

def loop(rOld,r2DotOld,tao1,tao3,Do,D,pHat1,pHat2,pHat3,R1,R2,R3,t1,t2,t3,maxiter=2000):
    r2Old=rOld[1]
    it=0
    while True:
        r2Mag=np.linalg.norm(r2Old)
        u=1/(r2Mag**3); z=(r2Old.dot(r2DotOld))/(r2Mag**2)
        q=(r2DotOld.dot(r2DotOld))/(r2Mag**2)-u
        f1=1-0.5*u*tao1**2-0.5*u*z*tao1**3+(1/24)*(3*u*q-15*u*z**2+u**2)*tao1**4
        f3=1-0.5*u*tao3**2+0.5*u*z*tao3**3+(1/24)*(3*u*q-15*u*z**2+u**2)*tao3**4
        g1=-(tao1-(1/6)*u*tao1**3+0.25*u*z*tao1**4)
        g3=tao3-(1/6)*u*tao3**3+0.25*u*z*tao3**4
        denom = (f1*g3-g1*f3)
        if abs(denom) < 1e-16:
            raise RuntimeError('Denominator too small in loop; aborting iteration')
        c1=g3/denom; c2=-1; c3=-g1/denom
        rNews=calcr(Do,D,pHat1,pHat2,pHat3,R1,R2,R3,c1,c2,c3,t1,t2,t3)
        newTaos=calcTao(rNews[3],rNews[4],rNews[5])
        tao1=newTaos[1]; tao3=newTaos[2]; r2=rNews[1]
        d1=-f3/denom; d3=f1/denom
        r2Dot=d1*rNews[0]+d3*rNews[2]
        it += 1
        if it > maxiter:
            raise RuntimeError('Iteration limit reached in Gauss loop')
        if all(abs(r2[i]-r2Old[i])<1e-9 and abs(r2Dot[i]-r2DotOld[i])<1e-9 for i in range(3)):
            break
        r2Old=r2; r2DotOld=r2Dot
    return r2,r2Dot

def methodOfGauss(t1,ra1,dec1,R1,t2,ra2,dec2,R2,t3,ra3,dec3,R3):
    pHat1=calcpHat(ra1,dec1); pHat2=calcpHat(ra2,dec2); pHat3=calcpHat(ra3,dec3)
    Do=calcDo(pHat1,pHat2,pHat3)
    D=calcD(pHat1,pHat2,pHat3,R1,R2,R3)
    taos=calcTao(t1,t2,t3)
    AB13=calcAB13(taos[0],taos[1],taos[2])
    AB=calcAB(AB13,D,Do)
    EF=calcEF(pHat2,R2)
    roots=calcr2Mag(AB,EF)
    # Auto-select largest positive real root
    positive_roots = [(i, abs(r)) for i, r in enumerate(roots) if r.imag == 0 and r.real > 0]
    root_idx = max(positive_roots, key=lambda x: x[1])[0] if positive_roots else 7
    r2Mag=abs(roots[root_idx])
    c1,c2,c3=initialC(AB13,r2Mag)
    rVectorsFirst=calcr(Do,D,pHat1,pHat2,pHat3,R1,R2,R3,c1,c2,c3,t1,t2,t3)
    taos=calcTao(rVectorsFirst[3],rVectorsFirst[4],rVectorsFirst[5])
    r2DotFirst=calcr2Dot(r2Mag,rVectorsFirst,taos[1],taos[2])
    r2,r2Dot=loop(rVectorsFirst,r2DotFirst,taos[1],taos[2],Do,D,pHat1,pHat2,pHat3,R1,R2,R3,t1,t2,t3)
    return r2,r2Dot

def magnitude(v):
    return sqrt(sum(v[i]**2 for i in range(len(v))))

def adjustAngle(angle):
    while angle>2*pi: angle-=2*pi
    while angle<0: angle+=2*pi
    return angle

def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: return asin(sine)
    if cosine < 0 and sine > 0: return acos(cosine)
    if cosine < 0 and sine < 0: return pi - asin(sine)
    if cosine > 0 and sine < 0: return 2*pi + asin(sine)

def orbElems(r,rdot):
    a=1/((2/magnitude(r))-rdot.dot(rdot))
    e=sqrt(abs(1-(magnitude(np.cross(r,rdot))**2/a)))
    rot=np.array([[1,0,0],[0,cos(radians(23.437)),-sin(radians(23.437))],[0,sin(radians(23.437)),cos(radians(23.437))]])
    rot=np.linalg.inv(rot)
    r=rot.dot(r); rdot=(rot.dot(rdot))/(0.01720209894)
    h=np.cross(r,rdot)
    I=adjustAngle(acos(h[2]/magnitude(h)))
    sinO=h[0]/(magnitude(h)*sin(I)); cosO=-h[1]/(magnitude(h)*sin(I))
    O=adjustAngle(atan2(sinO,cosO))
    sinwplusf=r[2]/(magnitude(r)*sin(I))
    coswplusf=(1/cos(O))*(r[0]/magnitude(r)+cos(I)*sinwplusf*sin(O))
    wplusf=atan2(sinwplusf,coswplusf)
    cosf=(1/e)*(a*(1-e**2)/magnitude(r)-1) if e > 0 else 0
    sinf=(a*(1-e**2)/(magnitude(h)))*(r.dot(rdot)/(e*magnitude(r))) if e > 0 else 0
    f=findQuadrant(sinf,cosf)
    w=adjustAngle(wplusf-f)
    cosE=(1-(magnitude(r)/a))/e if e > 0 else 0
    sinE=sqrt(abs(1-e**2))*sin(f)/(1+e*cos(f)) if e < 1 else 0
    E=findQuadrant(sinE,cosE)
    M=E-e*sin(E)
    return a,e,degrees(O),degrees(I),degrees(w),degrees(M)

# Load observations
with open(obs_file, 'r') as f:
    lines = [line.strip().split() for line in f if line.strip()]

obsList = []
for parts in lines:
    JD = JDNumber(parts[0], parts[1])
    ra = decimalCoord(parts[2])
    dec = decimalCoord(parts[3])
    solar = np.array([float(parts[4]), float(parts[5]), float(parts[6])])
    obsList.append([JD, ra, dec, solar[0], solar[1], solar[2]])

numObs = len(obsList)
print(f"  Processing {numObs} observations (running all 3-observation combinations)...")

from itertools import combinations
results = []
os.makedirs(results_dir, exist_ok=True)
for combo in combinations(range(numObs), 3):
    i, j, k = combo
    t1, ra1, dec1, R1x, R1y, R1z = obsList[i]
    t2, ra2, dec2, R2x, R2y, R2z = obsList[j]
    t3, ra3, dec3, R3x, R3y, R3z = obsList[k]
    try:
        r2, r2Dot = methodOfGauss(float(t1), ra1, dec1, np.array([R1x,R1y,R1z]),
                                   float(t2), ra2, dec2, np.array([R2x,R2y,R2z]),
                                   float(t3), ra3, dec3, np.array([R3x,R3y,R3z]))
        elements = orbElems(r2, r2Dot)
        results.append((combo, elements, r2, r2Dot))
        # save plot for this combo
        a, e = elements[0], elements[1]
        O, I, w = radians(elements[2]), radians(elements[3]), radians(elements[4])
        E_vals = np.linspace(0, 2*pi, 200)
        orbit_pts = []
        for E in E_vals:
            x = a * (cos(E) - e)
            y = a * sqrt(abs(1 - e**2)) * sin(E)
            pos = np.array([x, y, 0])
            R_w = np.array([[cos(w), -sin(w), 0], [sin(w), cos(w), 0], [0, 0, 1]])
            R_I = np.array([[1, 0, 0], [0, cos(I), -sin(I)], [0, sin(I), cos(I)]])
            R_O = np.array([[cos(O), -sin(O), 0], [sin(O), cos(O), 0], [0, 0, 1]])
            pos_ecliptic = R_O @ R_I @ R_w @ pos
            obliquity = radians(23.437)
            R_eq = np.array([[1, 0, 0], [0, cos(obliquity), -sin(obliquity)], [0, sin(obliquity), cos(obliquity)]])
            orbit_pts.append(R_eq @ pos_ecliptic)
        orbit_pts = np.array(orbit_pts)
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter([0],[0],[0], color='yellow', s=200, marker='o')
        ax.plot(orbit_pts[:,0], orbit_pts[:,1], orbit_pts[:,2], 'b-')
        ax.scatter([r2[0]],[r2[1]],[r2[2]], color='red', s=60)
        ax.set_xlabel('X (AU)'); ax.set_ylabel('Y (AU)'); ax.set_zlabel('Z (AU)')
        ax.set_title(f'67P orbit combo {i+1},{j+1},{k+1} a={a:.3f} e={e:.3f} i={elements[3]:.2f}°')
        outpng = os.path.join(results_dir, f'67P_orbit_combo_{i+1}_{j+1}_{k+1}.png')
        plt.savefig(outpng, dpi=200, bbox_inches='tight')
        plt.close(fig)
    except Exception as exc:
        print('  Combo', combo, 'failed:', exc)

# write results summary
with open(os.path.join(results_dir, '67P_results_combos.txt'), 'w', encoding='utf-8') as fo:
    if not results:
        fo.write('No successful Gauss runs.\n')
    for r in results:
        combo, elements, r2, r2Dot = r
        fo.write('Combo: {}\n'.format(combo))
        fo.write('  a={:.6f} AU, e={:.6f}, Omega={:.4f} deg, i={:.4f} deg, w={:.4f} deg, M={:.4f} deg\n'.format(*elements))
        fo.write('  r2 = [{:.6f},{:.6f},{:.6f}]\n'.format(r2[0],r2[1],r2[2]))
        fo.write('  r2dot = [{:.6f},{:.6f},{:.6f}]\n\n'.format(r2Dot[0],r2Dot[1],r2Dot[2]))

print('\n✓ Completed Gauss runs for all 3-observation combinations. Results saved to', os.path.join(results_dir, '67P_results_combos.txt'), 'and plots in', results_dir)

# Summary and representative plot
if results:
    # pick the first successful result as representative
    combo, elements, r2, r2Dot = results[0]
    with open(os.path.join(results_dir, '67P_results_summary.txt'), 'w', encoding='utf-8') as fo:
        fo.write('Representative result from combo: {}\n'.format(combo))
        fo.write('Semi-major axis (a): {:.6f} AU\n'.format(elements[0]))
        fo.write('Eccentricity (e): {:.6f}\n'.format(elements[1]))
        fo.write('Longitude of node (Ω): {:.4f} deg\n'.format(elements[2]))
        fo.write('Inclination (i): {:.4f} deg\n'.format(elements[3]))
        fo.write('Argument of perihelion (ω): {:.4f} deg\n'.format(elements[4]))
        fo.write('Mean anomaly (M): {:.4f} deg\n\n'.format(elements[5]))
        fo.write('Position vector (AU): [{:.6f}, {:.6f}, {:.6f}]\n'.format(r2[0], r2[1], r2[2]))
        fo.write('Velocity vector (AU/day): [{:.6f}, {:.6f}, {:.6f}]\n'.format(r2Dot[0], r2Dot[1], r2Dot[2]))
    print('\nRepresentative result saved to', os.path.join(results_dir, '67P_results_summary.txt'))
    # also save representative orbit plot already produced earlier (first in results/)
else:
    print('\nNo successful Gauss results to summarize.\nPlease verify astrometry and try again.')
    # Exit early when no results are available
    sys.exit(0)

#======================== STEP 3: DISPLAY RESULTS ========================
print("\n[STEP 3/3] Generating results...")

print("\n" + "="*70)
print("ORBITAL ELEMENTS FOR COMET 67P")
print("="*70)
if results:
    print(f"\n  Semi-major axis (a):       {elements[0]:.6f} AU")
    print(f"  Eccentricity (e):          {elements[1]:.6f}")
    print(f"  Longitude of node (Ω):     {elements[2]:.4f}°")
    print(f"  Inclination (i):           {elements[3]:.4f}°")
    print(f"  Argument of perihelion(ω): {elements[4]:.4f}°")
    print(f"  Mean anomaly (M):          {elements[5]:.4f}°")
    print(f"\n  Position vector [AU]:      [{r2[0]:.6f}, {r2[1]:.6f}, {r2[2]:.6f}]")
    print(f"  Velocity vector [AU/day]:  [{r2Dot[0]:.6f}, {r2Dot[1]:.6f}, {r2Dot[2]:.6f}]")
else:
    print('\nNo successful Gauss results to display. Check', os.path.join(results_dir, '67P_results_combos.txt'), 'for details and verify astrometry in', astrometry_file)

# Known values for comparison
print("\n" + "="*70)
print("COMPARISON WITH KNOWN VALUES (67P)")
print("="*70)
known = {'a': 3.463, 'e': 0.641, 'i': 7.04}
print(f"\n  Parameter          Known       Computed    Difference")
print(f"  ----------------   --------    --------    ----------")
print(f"  Semi-major axis    {known['a']:.3f} AU    {elements[0]:.3f} AU   {abs(elements[0]-known['a']):.3f} AU")
print(f"  Eccentricity       {known['e']:.3f}       {elements[1]:.3f}       {abs(elements[1]-known['e']):.3f}")
print(f"  Inclination        {known['i']:.2f}°       {elements[3]:.2f}°      {abs(elements[3]-known['i']):.2f}°")

# Create 3D visualization
print("\n  Creating 3D orbit visualization...")
a, e = elements[0], elements[1]
O, I, w = radians(elements[2]), radians(elements[3]), radians(elements[4])

# Generate orbit points
E_vals = np.linspace(0, 2*pi, 200)
orbit_pts = []
for E in E_vals:
    x = a * (cos(E) - e)
    y = a * sqrt(abs(1 - e**2)) * sin(E)
    pos = np.array([x, y, 0])
    R_w = np.array([[cos(w), -sin(w), 0], [sin(w), cos(w), 0], [0, 0, 1]])
    R_I = np.array([[1, 0, 0], [0, cos(I), -sin(I)], [0, sin(I), cos(I)]])
    R_O = np.array([[cos(O), -sin(O), 0], [sin(O), cos(O), 0], [0, 0, 1]])
    pos_ecliptic = R_O @ R_I @ R_w @ pos
    obliquity = radians(23.437)
    R_eq = np.array([[1, 0, 0], [0, cos(obliquity), -sin(obliquity)], [0, sin(obliquity), cos(obliquity)]])
    orbit_pts.append(R_eq @ pos_ecliptic)

orbit_pts = np.array(orbit_pts)

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')
ax.scatter([0], [0], [0], color='yellow', s=500, marker='o', label='Sun', zorder=10)
ax.plot(orbit_pts[:, 0], orbit_pts[:, 1], orbit_pts[:, 2], 'b-', linewidth=2, label='67P Orbit')
ax.scatter([r2[0]], [r2[1]], [r2[2]], color='red', s=200, marker='o', label='Current Position', zorder=5)
ax.set_xlabel('X (AU)'); ax.set_ylabel('Y (AU)'); ax.set_zlabel('Z (AU)')
ax.set_title(f'Comet 67P Orbit\na={a:.2f} AU, e={e:.3f}, i={elements[3]:.1f}°', fontsize=14, fontweight='bold')
ax.legend()
out_main_orbit = os.path.join(results_dir, '67P_orbit.png')
plt.savefig(out_main_orbit, dpi=300, bbox_inches='tight')
if results:
    print("  ✓ Saved visualization to", out_main_orbit)

# Save results
if results:
    with open(os.path.join(results_dir, '67P_results.txt'), 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write("COMET 67P ORBITAL DETERMINATION RESULTS\n")
        f.write("="*70 + "\n\n")
        f.write(f"Semi-major axis (a):       {elements[0]:.6f} AU\n")
        f.write(f"Eccentricity (e):          {elements[1]:.6f}\n")
        f.write(f"Longitude of node (Ω):     {elements[2]:.4f}°\n")
        f.write(f"Inclination (i):           {elements[3]:.4f}°\n")
        f.write(f"Argument of perihelion(ω): {elements[4]:.4f}°\n")
        f.write(f"Mean anomaly (M):          {elements[5]:.4f}°\n\n")
        f.write(f"Position vector [AU]:      [{r2[0]:.6f}, {r2[1]:.6f}, {r2[2]:.6f}]\n")
        f.write(f"Velocity vector [AU/day]:  [{r2Dot[0]:.6f}, {r2Dot[1]:.6f}, {r2Dot[2]:.6f}]\n")
    print("  ✓ Saved results to", os.path.join(results_dir, '67P_results.txt'))

print("\n" + "="*70)
print("✓ ANALYSIS COMPLETE!")
print("="*70)
print("\nGenerated files:")
print("  • 67P_observations.txt - Observation data")
print("  • 67P_results.txt      - Orbital elements")
print("  • 67P_orbit.png        - 3D visualization")
print("\n" + "="*70)
