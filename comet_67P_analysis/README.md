# PH556

This folder contains a small, self-contained pipeline to extract centroids from the provided FITS images, plate-solve them using the repository LSPR star lists, and compute a preliminary orbit for comet 67P using the classical Gauss method.

Goals and scope
- Produce a reproducible, one-command analysis for the four provided FITS files (Feb 2022).
- Keep the repository's original algorithms (LSPR, Gauss) and adapt them to run non-interactively where appropriate.
- Store all generated artifacts under `data/` and `results/` so the repository root stays clean.

High-level layout
- `centroid_batch.py` — automated flux-weighted centroid extraction (fast, non-interactive). Produces `data/comet_67P_centroids.txt`.
- `LSPR_fix.py` — small Python-3 compatible set of functions derived from the original `LSPR.py` so we can plate-solve programmatically.
- `run_lsps_py3_wrapper.py` — wrapper that runs plate-solving for each centroid and writes `data/comet_67P_astrometry.txt`.
- `analyze_comet_67p.py` — main pipeline that reads `data/` and writes all outputs under `results/` (combo results, plots, representative summary, final `67P_results.txt`).
- `run_pipeline.py` — helper that moves generated files into `data/` and `results/`, re-runs analysis, and writes `results/RESULTS.md`.
- `requirements.txt` — minimal third-party packages.

What you should know before running
- The four FITS files you provided are treated as the only raw observation inputs. Those FITS should remain in the repository root (this repo assumes that original placement).
- Automated centroiding is fast but can be less accurate than manual centroiding. For best scientific results, run the original `centroidprogram.py` interactively and then plate-solve.
- The Gauss method is sensitive to observation geometry. With just four observations over ~8 days the preliminary solution may differ substantially from published orbital elements; differential-correction and more observations reduce that error.

Quick start — reproduce everything (PowerShell)
1) Create and activate a venv (optional but recommended) and install deps:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

2) Run the whole pipeline (this runs centroiding, plate-solve, Gauss combos, and writes results):

```powershell
python .\comet_67P_analysis\run_pipeline.py
```

Exact PowerShell commands (copy/paste)
If you've cleared `data/` and `results/` and want the exact commands I used (run these from the repository root or directly from the `comet_67P_analysis` folder):

```powershell
# 1) create + activate venv (only once)
python -m venv .venv
. .\.venv\Scripts\Activate.ps1
pip install -r .\comet_67P_analysis\requirements.txt

# 2) regenerate centroids from the FITS in the analysis folder (preferred)
python .\comet_67P_analysis\centroid_batch.py --fits-dir .\comet_67P_analysis --out "comet_67P_centroids.txt"

# 3) move the centroid CSV into the pipeline data/ folder
Move-Item -Path .\comet_67P_centroids.txt -Destination .\comet_67P_analysis\data\comet_67P_centroids.txt

# 4) run automated plate-solve (creates data/comet_67P_astrometry.txt)
python .\comet_67P_analysis\run_lsps_py3_wrapper.py

# 5) run full pipeline (centroids -> plate-solve -> Gauss -> plots -> results)
python .\comet_67P_analysis\run_pipeline.py

# 6) quick checks (list data and results, preview RESULTS.md and combos)
Get-ChildItem .\comet_67P_analysis\data
Get-ChildItem .\comet_67P_analysis\results
Get-Content .\comet_67P_analysis\results\RESULTS.md -TotalCount 200
Get-Content .\comet_67P_analysis\results\67P_results_combos.txt -Raw
```

Notes:
- If the plate-solve step fails, check `comet_67P_analysis\data\*.lspr.out.txt` for the LSPR output and errors.
- If Gauss fails for all combos, open `comet_67P_analysis\data\comet_67P_astrometry.txt` and verify that the RA/DEC lines correspond to the comet centroid (not the field center). If they look wrong, run the interactive `centroidprogram.py` to vet centroids and repeat step 3.

What that single command does
When you run `python .\comet_67P_analysis\run_pipeline.py` it will perform the full, automated workflow:

- Move any previously-generated temporary files into `data/` and `results/` so the repository root remains clean.
- Run centroid extraction (if `data/comet_67P_centroids.txt` isn't present) and write centroids to `data/`.
- Plate-solve each centroid using the star lists and write `data/comet_67P_astrometry.txt`.
- Run the Gauss method on every 3-observation combination, saving combo outputs to `results/67P_results_combos.txt` and a representative `results/67P_results_summary.txt`.
- Generate plots (combo orbit PNGs and a main `results/67P_orbit.png`).
- Write a human-ready `results/RESULTS.md` (summary + embedded images).

How to regenerate a specific artifact (single steps)
- Regenerate centroids only (fast):

```powershell
python .\centroid_batch.py --fits-dir ".." --out comet_67P_centroids.txt
move ..\comet_67P_centroids.txt .\data\
```

- Regenerate plate-solved astrometry (after centroids are ready):

```powershell
python .\run_lsps_py3_wrapper.py
# writes data/comet_67P_astrometry.txt and per-image .lspr.out.txt files
```

- Re-run only the Gauss analysis and save outputs to results/ (keeps data/ intact):

```powershell
python .\analyze_comet_67p.py
# writes results/67P_results_combos.txt, results/67P_results_summary.txt, results/67P_orbit.png, results/67P_results.txt
```

Reproducibility notes (commands and file locations)
- Centroids: `data/comet_67P_centroids.txt` — columns: filename, x_pix, y_pix, flux_sum, sigma_x, sigma_y.
- Plate-solved astrometry: `data/comet_67P_astrometry.txt` — columns: filename, RA(h:m:s), DEC(d:m:s), lspr_output_file.
- Observations used by Gauss: `data/67P_observations.txt` — date time RA DEC SunX SunY SunZ.
- Gauss outputs: `results/67P_results_combos.txt` and `results/67P_results_summary.txt` and `results/67P_results.txt`.
- Visualizations: `results/67P_orbit.png` and `results/67P_orbit_combo_*.png`.

Interpretation hints (what to look for in the results)
- Semi-major axis and eccentricity: preliminary orbit will often differ from JPL/IAU catalog values when using a small number of observations.
- Inclination and node: small-number geometry may bias these angles; compare with known values (provided in this README) to assess reliability.
- RMS/repeatability: check `results/67P_results_combos.txt` to see how different 3-observation combinations change the derived elements — large spread indicates the solution is geometry-limited.

If you want better accuracy (recommended)
1. Open each FITS with `centroidprogram.py` and pick/verify the comet centroid interactively.
2. Re-run `run_lsps_py3_wrapper.py` to produce `data/comet_67P_astrometry.txt`.
3. Run `analyze_comet_67p.py` to obtain a cleaned Gauss solution.
4. Ask me to add differential correction (least-squares) and jackknife uncertainty — I can wire those to produce uncertainties and improved fits.

Troubleshooting
- If `analyze_comet_67p.py` errors on iteration: the Gauss loop has safety checks and will abort if the geometry is degenerate. Check `data/comet_67P_astrometry.txt` for bad coordinates.
- If plate-solve fails: ensure the star lists (`LSPRstars1.txt`, `LSPRstars2.txt`) are present in the parent repo directory.

Contact / next steps
- If you want, I can:
  - add differential correction + jackknife and append the outputs to `results/RESULTS.md`;
  - help you vet centroids interactively using the original `centroidprogram.py`;
  - or package the `comet_67P_analysis` folder as a clean zip for upload.

Thank you — the pipeline is designed to be reproducible and to keep your repo root uncluttered. If you'd like more scientifically rigorous steps (interactive centroiding, plate-solve diagnostics, or LSQ refinements), tell me which and I will implement them next.
# Comet 67P Orbital Determination

Python 3 implementation of the Gauss Method for determining orbital elements of comet 67P/Churyumov-Gerasimenko from archival FITS observations.

## Installation

```bash
pip install numpy astropy matplotlib
```

## Quick Start

```bash
python analyze_comet_67p.py
```

This single command will:
1. Extract observation data from 4 FITS files (Feb 2022)
2. Compute orbital elements using the Gauss method
3. Generate 3D orbit visualization
4. Save results to text file

## Output Files

- **`67P_observations.txt`** - Extracted observation data (date, time, RA, DEC, Sun-Earth vector)
- **`67P_results.txt`** - Computed orbital elements with position/velocity vectors
- **`67P_orbit.png`** - 3D visualization of orbital trajectory

## Orbital Elements Computed

- **Semi-major axis (a)** - Size of the orbit in AU
- **Eccentricity (e)** - Shape of orbit (0=circle, <1=ellipse)
- **Inclination (i)** - Tilt of orbital plane in degrees
- **Longitude of ascending node (Ω)** - Where orbit crosses ecliptic
- **Argument of perihelion (ω)** - Orientation within orbital plane
- **Mean anomaly (M)** - Position along orbit at epoch

## Known Values for Comet 67P

For validation, the known orbital elements are:
- Semi-major axis: 3.463 AU
- Eccentricity: 0.641
- Inclination: 7.04°
- Orbital period: ~6.45 years

## FITS Files Required

You can place the four FITS either in the repository root (parent directory `../`) or directly inside the `comet_67P_analysis/` folder. The pipeline will prefer FITS located inside `comet_67P_analysis/` if present and fall back to the parent directory otherwise.

Expected filenames (the pipeline's star-list mapping relies on these exact names):
- `20220217191741-313-RA.wcs.proc.fits`
- `20220220190400-656-RA.wcs.proc.fits`
- `20220221190037-229-RA.wcs.proc.fits`
- `20220225192324-727-RA.wcs.proc.fits`

Caveats:
- The plate-solver wrapper (`run_lsps_py3_wrapper.py`) uses a hard-coded mapping from these filenames to the local `LSPRstars1.txt` / `LSPRstars2.txt` star lists located in the repository root. If you rename the FITS files, update `run_lsps_py3_wrapper.py`'s `star_map` accordingly.
- Make sure the star lists (`LSPRstars1.txt`, `LSPRstars2.txt`) remain available in the repository root; the wrapper looks for them there.

## Method

Uses the classical **Gauss Method** for orbit determination:
1. Extracts 3 observations (first, middle, last)
2. Computes direction vectors from Earth to comet
3. Solves 8th-degree polynomial for range
4. Iteratively refines using Taylor series
5. Converts position/velocity to orbital elements

## Notes

- Observation data uses estimated comet positions from ephemeris
- For best accuracy, precise astrometric measurements from FITS images would be needed
- Runtime: ~1-2 seconds on modern hardware

## References

- Original implementation: github.com/jkim117/Method-of-Gauss-Orbital-Determination
- Gauss method for orbital determination (classical celestial mechanics)
- Modernized from Python 2.7 to Python 3 (2025)

## Troubleshooting

**Import errors**: Run `pip install numpy astropy matplotlib`

**FITS files not found**: Ensure FITS files are in parent directory

**Unrealistic results**: Verify FITS file locations and observation data quality

---

**Complete pipeline in one command:**
```bash
python analyze_comet_67p.py
```

Results saved to `67P_results.txt` and visualization to `67P_orbit.png`
