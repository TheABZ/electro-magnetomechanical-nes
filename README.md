# Electro-magnetomechanical NES

Foundational Python time-domain solver for the reduced 2-DOF equations in
Cadiou et al., *Damping adjustment of a nonlinear vibration absorber using an
electro–magnetomechanical coupling*, Journal of Sound and Vibration 518 (2022),
116508. [Paper DOI](https://doi.org/10.1016/j.jsv.2021.116508).

This implements Eq. (10) as supplied. It is a starting point for replication,
not a validated reproduction of the paper's figures or analytical branches.

## Quick start

Use Python 3.10 or newer:

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python nes.py
python -m unittest discover -s tests -v
```

`nes.py` plots a demonstration response and prints the sampled primary peak
for t >= 20 s. Its mechanical values come from Table 1 and Rc from Table 2;
forcing, external resistance, and theta coefficients are illustrative choices.
Replace these before comparing with a particular experimental case.

## State-space formulation

Let `z = [x1, v1, x2, v2]`, `y = x1 - x2`, `u = v1 - v2`,
`S = k2l*y + k2c*y**3`, and `c_em = theta(y)**2 / (Rc + Ri)`. Then:

```text
dx1/dt = v1
dv1/dt = [Fe*sin(omega*t) - c1*v1 - k1*x1 - c2*u - S] / m1
dx2/dt = v2
dv2/dt = [(c2 + c_em)*u + S] / m2
```

The electromagnetic term occurs only in the second equation, exactly as in
the requested model and the paper's Eq. (10). Adding an equal and opposite
reaction to the primary mass would change this model; it is not done here.
The reduced model neglects coil inductance (paper Eq. (8)); a model retaining
current dynamics would require a fifth state.

Use SI units: masses in kg, displacements in m, time in s, damping in N s/m,
linear stiffness in N/m, cubic stiffness in N/m³, resistance in ohms, force
in N, and **omega in rad/s** (`omega = 2*pi*f_hz`). Theta has units N/A,
equivalently V s/m, so polynomial coefficients have units dependent on power.

`theta_coeffs` is in **ascending power order**. For the paper's Eq. (11),
`theta(y) = a1*y**8 + a2*y**6 + a3*y**4 + a4*y**2 + a5`, use:

```python
theta_coeffs = (a5, 0., a4, 0., a3, 0., a2, 0., a1)
```

Do not insert coefficients of theta squared or silently use millimetres.
Validate the fit over the simulated displacement range. `Ri=float('inf')`
represents open circuit; `Ri=0` represents a short across the external load
with the finite coil resistance retained.

## Integration and frequency sweeps

`simulate(parameters, time, initial_state)` returns a SciPy OdeResult.
`result.y[0]` is x1; `result.y[2]` is x2. Initial conditions are defined at
`time[0]`. Failed integrations raise RuntimeError rather than yielding an
apparently valid amplitude. Solver method, tolerances and maximum step size
are keyword arguments; see the
[SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html).

Example independent runs from rest (values below are illustrative):

```python
from dataclasses import replace
import numpy as np
from nes import Parameters, simulate, peak_amplitude

p = Parameters(m1=64., m2=1.4, c1=97., c2=1., k1=91324.,
               k2l=500., k2c=8.44e6, Rc=11.6, Ri=10., Fe=1.,
               omega=37., theta_coeffs=(1., 0., 100.))
omegas = np.linspace(30., 45., 16)
amplitudes = []
for omega in omegas:
    period = 2*np.pi / omega
    # Run 200 forcing periods; retain the final 100 periods.
    time = np.linspace(0., 200*period, 200*100 + 1)
    result = simulate(replace(p, omega=omega), time)
    amplitudes.append(peak_amplitude(result.t, result.y[0],
                                     transient_time=100*period))
np.savetxt('frequency_response.csv', np.column_stack([omegas, amplitudes]),
           delimiter=',', header='omega_rad_s,peak_x1_m', comments='')
```

Nonlinear responses can depend on initial conditions and sweep direction.
For continuation, pass the previous final state as the next initial state
and preserve the intended forcing phase (e.g. use whole forcing periods
when restarting the time origin). Compare upward and downward sweeps.
Time integration samples accessible attractors; it does not compute unstable
branches or reproduce the paper's analytical continuation automatically.

## Interpreting and checking amplitudes

`peak_amplitude` computes `max(abs(x1))` on retained samples. It is neither
RMS nor half peak-to-peak. By default, it discards the first half of elapsed
time. `transient_time` sets an absolute cutoff and overrides that fraction.
The cutoff does **not** automatically establish steady state.

Compare successive late-time windows, extend the simulation, and check that
amplitudes stabilize to your chosen tolerance. For strongly modulated or
chaotic responses, report the observation-window length and check longer
windows. Refine output sampling and solver tolerances separately: sampled
peaks can underestimate extrema between samples. The solver uses adaptive
internal steps; t_eval alone is not an internal step limit. The default
max_step here is capped by output spacing and 1/50 of the forcing period.

Tests check the original equation residuals, the exact uncoupled oscillator,
the analytical linear forced-response amplitude, polynomial ordering,
open-circuit behavior, transient selection, and invalid inputs.
