# Replication status and provenance

This is a **table-faithful numerical replication attempt**, not a verified
reproduction of the published frequency-response figures. The reduced
four-state equations are implemented literally as Eq. (10), including the
one-sided electromagnetic term. The underlying source is:

B. Cadiou, C. Stephan, A. Renoult, G. Michon, *Damping adjustment of a nonlinear
vibration absorber using an electro–magnetomechanical coupling*, Journal of
Sound and Vibration 518 (2022), 116508.
https://doi.org/10.1016/j.jsv.2021.116508

## Inputs and evidence

| Input | Value / source | Status |
|---|---|---|
| Primary mass, NES mass | 64 kg, 1.4 kg; Table 1 | Reported |
| Primary damping, NES mechanical damping | 97, 1 N s/m; Table 1 | Reported |
| Primary stiffness | 91,324 N/m; Table 1 | Reported; inconsistent with plotted reference resonance |
| NES linear, cubic stiffness | 500 N/m, 8.44e6 N/m³; Table 1 | Reported |
| Coil resistance | 11.6 ohms; Table 2 | Reported |
| External resistance | infinity for open, zero for closed | Section 4.2 definitions |
| Excitation | 5 N baseline; 25 N open/closed | Figures 11, 16, 17 |
| Theta polynomial | Even powers through degree 8; Eq. (11) | Coefficients reconstructed from Figure 8 |
| Sweep | 4–7 Hz in 0.05 Hz increments | Our discrete dwell protocol, extended past paper's 4–6 Hz range |
| Initial conditions | Rest at each sweep's first point; final state passed forward | Our protocol; not author-supplied |
| Duration | 160 periods per point, final 80 analyzed | Our numerical choice; convergence diagnostics supplied |

No numeric values for a1–a5 were found in the supplied paper. A web search
identified the author's thesis, but a download attempt returned HTTP 503;
no coefficient values from that unavailable source were used.

## Coupling reconstruction

`data/figure8_digitized.csv` holds 837 points recovered from the red
approximation in the embedded 1000×773 raster image on PDF page 6. The
extraction excludes the legend and uses the median red pixel row in each
column. Axes are calibrated at x=128/970 pixels for -40/+40 mm, and
image rows 16/649 for 20/0 V s/m. Exact calibration and source PDF SHA-256
are saved in `data/coupling_fit.json`.

The fit uses q=y/0.04 and an even degree-8 least-squares polynomial to avoid
ill-conditioning from fitting metre powers directly. Ascending SI coefficients
are stored in the JSON file; odd coefficients are exactly zero. Fit residual
RMSE is approximately 0.0165 V s/m and maximum residual 0.0698 V s/m.
These residuals are not uncertainty bounds on the authors' actual coefficients.
Raster line thickness and axis calibration introduce additional uncertainty.
The central coupling is approximately 12.564 V s/m.

The fit is supported only for **-0.04 <= y <= 0.04 m**. Every simulation CSV
reports maximum relative displacement over the entire trajectory and a range
flag. Treat closed-circuit runs outside that interval as unsupported polynomial
extrapolation; open-circuit dynamics do not use theta. No clipping is applied.

To reconstruct it again from your own copy (optional tooling, not needed for
normal simulation):

```bash
python -m pip install -r requirements-digitize.txt
python scripts/digitize_coupling.py /path/to/paper.pdf
```

The source PDF and its figures are not redistributed in this repository.
The coupling chart is newly drawn from derived numerical points.

## A source inconsistency blocks exact agreement

For a locked NES, the two masses move together, so the linear reference is

`A = Fe / sqrt((k1 - (m1+m2)*omega**2)**2 + (c1*omega)**2)`.

The damped displacement resonance is
`sqrt(k1/(m1+m2) - c1**2/(2*(m1+m2)**2)) / (2*pi)`.
The published table gives approximately **5.945 Hz**. Figure 17's locked
reference peaks near **4.85 Hz**, and the paper describes approximately
8 mm at 25 N. The listed parameters instead give a peak near 6.90 mm.
Both the table and plotted axis were visually checked against the PDF.
Removing the NES mass would not resolve this discrepancy.

Consequently, shifting the frequency axis or changing stiffness until the
curves overlap would hide an unresolved input discrepancy. No such tuning
is performed. Confirmation of the actual mass/stiffness for the plotted
experiments (or original author code) is needed for a defensible exact match.
The Figure 17 caption also says Rinf, whereas its legend and Section 4.2
identify that case as Rmin; this implementation follows the legend and text.

## Running and reading the results

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python reproduce.py --workers 4
```

The five sweeps are a 5 N open-circuit baseline, and 25 N open/closed circuits
in both directions. Directions start independently from rest. Each frequency
is held for an integer number of forcing periods, preserving forcing phase
when passing the final state to the next point. This is a stepped dwell
sweep, **not** the experimental continuous chirp; the paper does not provide
a sweep-rate protocol sufficient to reproduce its experimental trajectories.

- `results/frequency_response.png`: newly computed sampled peak curves and
  analytic locked reference. The paper's approximate reference frequency is
  marked to expose the mismatch, not used to calibrate parameters.
- `results/*N_*.csv`: frequency, sampled absolute peak, fundamental amplitude,
  RMS, change between the two late-window peaks, maximum relative displacement,
  and polynomial range flag.
- `results/coupling_fit.png`: digitized red curve and reconstructed polynomial.
- `results/time_response.png` and its CSV: 25 N closed-circuit response at
  5.95 Hz from rest. The plot shows the complete trajectory; the CSV retains
  the final five forcing periods at the original output resolution.
- `results/convergence.json`: separate output-resolution/tolerance refinement
  and duration-extension comparisons at 4.85, 5.95 and 6.2 Hz, both circuits.
- `results/summary.json`: exact settings, library versions, maxima and flags.
- `results/RESULTS.md`: concise numerical findings from the completed run.

`peak_x1_m` is sampled max(abs(x1)); `fundamental_x1_m` is the amplitude of
the Fourier component at the forcing frequency; neither should be assumed
equal to the paper's averaged fixed-point amplitude in modulated motion.
The final half of each dwell is analyzed. A late-window change over 1% is a
diagnostic, not a proof of instability; agreement below 1% is not proof that
all slow modulation has decayed. Convergence checks start from rest and do
not certify every continuation point or alternative attractor.

Tighter settings can be run with `--cycles 320 --samples 160 --points 121`.
This replaces the results directory outputs; preserve earlier CSVs if needed.
The script does not claim to find all attractors, unstable fixed points,
Floquet multipliers, or analytical continuation branches.

## What remains scientifically unresolved

1. Resolve the table/figure resonance inconsistency with original author data.
2. Replace the digitized fit with author coefficients if available.
3. Obtain the experimental sweep rate and processing convention before making
   quantitative experimental-error claims.
4. If full analytical branches are required, implement and validate harmonic
   balance/continuation separately; direct integration cannot recover unstable
   branches.

No email has been sent to the authors. These missing inputs cannot be supplied
by additional numerical precision alone.
