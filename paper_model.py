"""Traceable Table 1/2 parameters and reconstructed Figure 8 coupling."""
import json
from pathlib import Path

import numpy as np
from nes import Parameters

ROOT = Path(__file__).resolve().parent


def coupling_fit():
    return json.loads((ROOT/'data/coupling_fit.json').read_text())


def paper_parameters(force=25., frequency_hz=6., circuit='closed'):
    if circuit not in ('open', 'closed'):
        raise ValueError("circuit must be 'open' or 'closed'")
    return Parameters(m1=64., m2=1.4, c1=97., c2=1., k1=91324.,
                      k2l=500., k2c=8.44e6, Rc=11.6,
                      Ri=np.inf if circuit == 'open' else 0., Fe=force,
                      omega=2*np.pi*frequency_hz,
                      theta_coeffs=tuple(coupling_fit()['theta_coeffs']))


def locked_amplitude(frequency_hz, p):
    """Locked NES moves with the LO, so the total moving mass is m1+m2."""
    omega = 2*np.pi*np.asarray(frequency_hz)
    return abs(p.Fe)/np.hypot(p.k1-(p.m1+p.m2)*omega**2, p.c1*omega)


def locked_peak(p):
    mass = p.m1+p.m2
    omega_peak = np.sqrt(p.k1/mass-p.c1**2/(2*mass**2))
    return float(omega_peak/(2*np.pi)), float(locked_amplitude(omega_peak/(2*np.pi), p))


def calibrated_parameters(force=25., frequency_hz=4.85, circuit='closed', target_hz=4.85):
    """Infer k1 from the approximate plotted locked displacement-resonance peak.

    Masses and damping remain reported values. This is a calibration hypothesis,
    not a verified correction or an author-supplied parameter.
    """
    from dataclasses import replace
    if not np.isfinite(target_hz) or target_hz <= 0:
        raise ValueError('target_hz must be positive and finite')
    p = paper_parameters(force, frequency_hz, circuit)
    mass = p.m1+p.m2
    k1 = mass*(2*np.pi*target_hz)**2 + p.c1**2/(2*mass)
    return replace(p, k1=k1)
