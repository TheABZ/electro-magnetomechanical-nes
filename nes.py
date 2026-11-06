"""Reduced 2-DOF NES model; state order is [x1, v1, x2, v2].

Implements the supplied equations (Cadiou et al., Eq. 10) literally:
electromagnetic damping acts in the secondary equation only.
All inputs use SI units, and forcing frequency is angular frequency (rad/s).
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Parameters:
    """Mechanical constants, circuit resistance, forcing, and theta polynomial.

    theta_coeffs uses ascending powers: (b0, b1, b2, ...) represents
    theta(y) = b0 + b1*y + b2*y**2 + ... with y in metres.
    Specify coefficients for theta itself, not theta squared.
    Ri=np.inf represents an open circuit; all other values must be finite.
    """
    m1: float
    m2: float
    c1: float
    c2: float
    k1: float
    k2l: float
    k2c: float
    Rc: float
    Ri: float
    Fe: float
    omega: float
    theta_coeffs: tuple[float, ...]

    def __post_init__(self):
        for name in ('m1', 'm2', 'c1', 'c2', 'k1', 'k2l', 'k2c',
                     'Rc', 'Fe', 'omega'):
            if not np.isfinite(getattr(self, name)):
                raise ValueError(f'{name} must be finite')
        if self.m1 <= 0 or self.m2 <= 0:
            raise ValueError('Masses must be positive')
        if min(self.c1, self.c2, self.Rc, self.omega) < 0:
            raise ValueError('Damping, Rc and omega must be nonnegative')
        if np.isnan(self.Ri) or self.Ri < 0 or self.Rc + self.Ri <= 0:
            raise ValueError('Ri must be nonnegative and Rc + Ri positive')
        coeffs = np.asarray(self.theta_coeffs, dtype=float)
        if coeffs.ndim != 1 or coeffs.size == 0 or not np.all(np.isfinite(coeffs)):
            raise ValueError('theta_coeffs must be a nonempty finite 1-D sequence')
        object.__setattr__(self, 'theta_coeffs', tuple(coeffs))


def theta(relative_displacement, coefficients):
    """Evaluate theta(y) using coefficients ordered from constant upwards."""
    return np.polynomial.polynomial.polyval(relative_displacement, coefficients)


def rhs(t, state, p):
    """Return [v1, a1, v2, a2] for the supplied two second-order ODEs."""
    x1, v1, x2, v2 = state
    y = x1 - x2
    relative_velocity = v1 - v2
    spring_force = p.k2l * y + p.k2c * y**3
    mechanical_force = p.c2 * relative_velocity + spring_force
    # Avoid evaluating a potentially large polynomial for an open circuit.
    c_em = 0.0 if np.isinf(p.Ri) else theta(y, p.theta_coeffs)**2 / (p.Rc + p.Ri)
    a1 = (p.Fe * np.sin(p.omega * t) - p.c1 * v1
          - p.k1 * x1 - mechanical_force) / p.m1
    # (x2 - x1)^3 = -y^3, so both spring terms enter a2 positively.
    a2 = (mechanical_force + c_em * relative_velocity) / p.m2
    return np.array([v1, a1, v2, a2])


def _time_array(time):
    time = np.asarray(time, dtype=float)
    if (time.ndim != 1 or time.size < 2 or not np.all(np.isfinite(time))
            or not np.all(np.diff(time) > 0)):
        raise ValueError('time must contain at least two finite, increasing times')
    return time


def simulate(p, time, initial_state=(0., 0., 0., 0.), *, method='DOP853',
             rtol=1e-8, atol=1e-10, max_step=None):
    """Integrate at requested times and return SciPy's OdeResult.

    Initial conditions apply at time[0]. Rows of result.y are x1, v1, x2, v2.
    By default, cap internal steps at the smallest output interval and 1/50
    of the forcing period. Sampling and solver accuracy still need convergence
    checks. Use method='Radau' or 'BDF' if the system becomes stiff.
    """
    time = _time_array(time)
    initial_state = np.asarray(initial_state, dtype=float)
    if initial_state.shape != (4,) or not np.all(np.isfinite(initial_state)):
        raise ValueError('initial_state must contain four finite values')
    if max_step is None:
        max_step = float(np.min(np.diff(time)))
        if p.omega > 0:
            max_step = min(max_step, 2 * np.pi / (50 * p.omega))
    if np.isnan(max_step) or max_step <= 0:
        raise ValueError('max_step must be positive')
    result = solve_ivp(rhs, (time[0], time[-1]), initial_state, args=(p,),
                       t_eval=time, method=method, rtol=rtol, atol=atol,
                       max_step=max_step)
    if not result.success:
        raise RuntimeError(f'Integration failed: {result.message}')
    if not np.all(np.isfinite(result.y)):
        raise RuntimeError('Integration returned nonfinite state values')
    return result


def peak_amplitude(time, x1, *, transient_fraction=0.5, transient_time=None):
    """Return sampled max(abs(x1)) after a user-selected transient cutoff.

    transient_time is an absolute time and overrides transient_fraction.
    The fraction removes that fraction of elapsed duration, not sample count.
    This is an amplitude relative to zero, not RMS or half peak-to-peak.
    It does not establish that steady state has been reached. Refine sampling
    and compare longer windows, especially for modulated/chaotic responses.
    """
    time = _time_array(time)
    x1 = np.asarray(x1, dtype=float)
    if x1.shape != time.shape or not np.all(np.isfinite(x1)):
        raise ValueError('x1 must be finite and have the same shape as time')
    if transient_time is None:
        if not np.isfinite(transient_fraction) or not 0 <= transient_fraction < 1:
            raise ValueError('transient_fraction must be in [0, 1)')
        cutoff = time[0] + transient_fraction * (time[-1] - time[0])
    else:
        cutoff = float(transient_time)
        if not np.isfinite(cutoff) or not time[0] <= cutoff < time[-1]:
            raise ValueError('transient_time must lie in [time[0], time[-1])')
    retained = x1[time >= cutoff]
    if retained.size < 2:
        raise ValueError('At least two samples must remain after the cutoff')
    return float(np.max(np.abs(retained)))


def main():
    """Run a demonstration, not a calibrated reproduction of a paper figure."""
    import matplotlib.pyplot as plt

    # Mechanical constants from Table 1 and Rc from Table 2.
    # Fe, Ri, omega and theta below are demonstration choices, NOT a fitted law.
    p = Parameters(m1=64., m2=1.4, c1=97., c2=1., k1=91324.,
                   k2l=500., k2c=8.44e6, Rc=11.6, Ri=10., Fe=1.,
                   omega=np.sqrt(91324. / 64.), theta_coeffs=(1., 0., 100.))
    time = np.linspace(0., 30., 15001)
    result = simulate(p, time)
    amplitude = peak_amplitude(result.t, result.y[0], transient_time=20.)
    print(f'Demo sampled peak |x1| for t >= 20 s: {amplitude:.6g} m')
    fig, ax = plt.subplots()
    ax.plot(result.t, result.y[0], label='Primary mass x1')
    ax.plot(result.t, result.y[2], label='NES x2', alpha=0.7)
    ax.axvline(20., color='black', linestyle='--', label='Transient cutoff')
    ax.set(xlabel='Time (s)', ylabel='Displacement (m)',
           title='Demonstration parameters — not a calibrated paper result')
    ax.legend()
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
