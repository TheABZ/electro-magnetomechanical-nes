"""Reproducible table-faithful time-domain study; see REPLICATION.md.

python reproduce.py --workers 4
Uses continuation within each sweep, and independent workers across sweeps.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy

from nes import simulate, peak_amplitude, theta
from paper_model import coupling_fit, paper_parameters, locked_amplitude, locked_peak

OUT = Path(__file__).resolve().parent/'results'
FIELDS = ['frequency_hz', 'peak_x1_m', 'fundamental_x1_m', 'rms_x1_m',
          'late_window_peak_change', 'max_relative_displacement_m', 'fit_range_ok']


def run_point(p, initial, cycles=160, samples=80, rtol=1e-8):
    period = 2*np.pi/p.omega
    time = np.linspace(0., cycles*period, cycles*samples+1)
    sol = simulate(p, time, initial, rtol=rtol, atol=rtol/100,
                   max_step=period/40)
    start = cycles*samples//2
    x = sol.y[0, start:]
    t = sol.t[start:]
    # Integer-period window: Fourier amplitude is distinct from sampled peak.
    integral = getattr(np, 'trapezoid', np.trapz if hasattr(np, 'trapz') else None)
    s = 2/(t[-1]-t[0])*integral(x*np.sin(p.omega*t), t)
    c = 2/(t[-1]-t[0])*integral(x*np.cos(p.omega*t), t)
    half = len(x)//2
    a, b = np.max(np.abs(x[:half])), np.max(np.abs(x[half:]))
    relative = float(np.max(np.abs(sol.y[0]-sol.y[2])))
    metrics = [p.omega/(2*np.pi), peak_amplitude(sol.t, sol.y[0]),
               float(np.hypot(s, c)), float(np.sqrt(integral(x*x, t)/(t[-1]-t[0]))),
               float(abs(a-b)/max(a, b, 1e-15)), relative, int(relative <= .04)]
    return metrics, sol


def sweep(task):
    force, circuit, direction, cycles, samples, count = task
    frequencies = np.linspace(4., 7., count)
    if direction == 'down':
        frequencies = frequencies[::-1]
    state = np.zeros(4)
    rows = []
    for f in frequencies:
        row, sol = run_point(paper_parameters(force, f, circuit), state, cycles, samples)
        rows.append(row)
        # Each dwell is an integer number of forcing periods: phase resets to zero.
        state = sol.y[:, -1]
    name = f'{force:g}N_{circuit}_{direction}'
    np.savetxt(OUT/f'{name}.csv', rows, delimiter=',', header=','.join(FIELDS), comments='')
    return name, rows


def checks(cycles, samples):
    rows = []
    for circuit in ('open', 'closed'):
        for f in (4.85, 5.95, 6.2):
            p = paper_parameters(25., f, circuit)
            base, sol = run_point(p, np.zeros(4), cycles, samples)
            refined, _ = run_point(p, np.zeros(4), cycles, samples*2, 1e-10)
            extended, _ = run_point(p, np.zeros(4), cycles*2, samples, 1e-8)
            rows.append(dict(circuit=circuit, frequency_hz=f,
                             peak_m=base[1], refined_peak_m=refined[1], extended_peak_m=extended[1],
                             resolution_relative_change=abs(refined[1]-base[1])/base[1],
                             duration_relative_change=abs(extended[1]-base[1])/base[1],
                             late_window_relative_change=base[4], fit_range_ok=bool(base[-1])))
            if circuit == 'closed' and f == 5.95:
                np.savetxt(OUT/'time_response_last5cycles_25N_closed_5p95Hz.csv',
                           np.column_stack([sol.t, sol.y.T])[-5*samples-1:],
                           delimiter=',', fmt='%.10e',
                           header='time_s,x1_m,v1_m_per_s,x2_m,v2_m_per_s', comments='')
                fig, ax = plt.subplots(figsize=(9, 4))
                ax.plot(sol.t, sol.y[0]*1000, label='Primary')
                ax.plot(sol.t, sol.y[2]*1000, label='NES', alpha=.6)
                ax.axvline(sol.t[-1]/2, color='black', ls='--', label='Analysis cutoff')
                ax.set(xlabel='Time (s)', ylabel='Displacement (mm)', title='25 N, closed circuit, 5.95 Hz; Table 1 parameters')
                ax.legend(); fig.tight_layout(); fig.savefig(OUT/'time_response.png', dpi=160); plt.close(fig)
    (OUT/'convergence.json').write_text(json.dumps(rows, indent=2)+'\n')
    return rows


def plot(results):
    points = np.loadtxt(Path(__file__).resolve().parent/'data/figure8_digitized.csv', delimiter=',', skiprows=1)
    y = np.linspace(-.04, .04, 1001)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(points[:, 0]*1000, points[:, 1], '.', ms=2, label='Digitized Figure 8 red curve')
    ax.plot(y*1000, theta(y, coupling_fit()['theta_coeffs']), label='Reconstructed even degree-8 fit')
    ax.set(xlabel='Relative displacement (mm)', ylabel='Coupling (V s/m)', title='Reconstructed coupling, valid over ±40 mm')
    ax.legend(); fig.tight_layout(); fig.savefig(OUT/'coupling_fit.png', dpi=160); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    f = np.linspace(4, 7, 1201)
    for ax, force in zip(axes, (5, 25)):
        ax.plot(f, locked_amplitude(f, paper_parameters(force))*1000, 'k--', label='Locked NES (analytic)')
        for name, rows in results:
            if name.startswith(f'{force}N_'):
                a = np.asarray(rows)
                ax.plot(a[:, 0], a[:, 1]*1000, '.-', ms=3, label=name.split('N_')[1].replace('_', ' '))
        ax.axvline(4.85, color='gray', alpha=.5, ls=':', label='Paper plotted peak ≈4.85 Hz')
        ax.set(xlabel='Frequency (Hz)', ylabel='Sampled peak |x1| (mm)', title=f'{force} N: published Table 1 constants')
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Time-domain comparison — parameter mismatch prevents exact figure replication')
    fig.tight_layout(); fig.savefig(OUT/'frequency_response.png', dpi=160); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--cycles', type=int, default=160)
    ap.add_argument('--samples', type=int, default=80)
    ap.add_argument('--points', type=int, default=61)
    args = ap.parse_args()
    if args.cycles < 8 or args.cycles % 4 or args.samples < 20 or args.points < 3 or args.workers < 1:
        ap.error('cycles must be divisible by 4 and >=8; samples>=20, points>=3, workers>=1')
    OUT.mkdir(exist_ok=True)
    tasks = [(5, 'open', 'up', args.cycles, args.samples, args.points)]
    tasks += [(25, circuit, direction, args.cycles, args.samples, args.points)
              for circuit in ('open', 'closed') for direction in ('up', 'down')]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(sweep, tasks))
    plot(results)
    validation = checks(args.cycles, args.samples)
    summary = dict(settings=vars(args), python=platform.python_version(), numpy=np.__version__,
                   scipy=scipy.__version__, matplotlib=matplotlib.__version__,
                   source='DOI 10.1016/j.jsv.2021.116508; Tables 1/2, digitized Figure 8',
                   locked_25N_peak_hz_m=locked_peak(paper_parameters()),
                   paper_reference='Figures 16/17: locked peak about 8 mm near 4.85 Hz; coupled about 3.8 mm (analytical)',
                   comparisons={}, convergence_checks=validation)
    for name, rows in results:
        a = np.asarray(rows); peak = a[np.argmax(a[:, 1])]
        summary['comparisons'][name] = dict(max_peak_m=float(peak[1]), frequency_hz=float(peak[0]),
            points_with_window_change_over_1pct=int(np.sum(a[:, 4]>.01)),
            points_outside_fit_range=int(np.sum(a[:, -1]==0)), max_relative_displacement_m=float(np.max(a[:, 5])))
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    from scripts.summarize_results import write_report
    write_report(summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
