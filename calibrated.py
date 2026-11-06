"""Run a stiffness-calibrated comparison without replacing published inputs.

python calibrated.py --workers 4
Outputs are separate from the original study in results/calibrated/.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy

from paper_model import paper_parameters, calibrated_parameters, locked_amplitude, locked_peak
from reproduce import run_point, FIELDS

OUT = Path(__file__).resolve().parent/'results/calibrated'


def sweep(task):
    force, circuit, direction, target, cycles, samples, count = task
    frequencies = np.linspace(4., 6., count)
    if direction == 'down':
        frequencies = frequencies[::-1]
    rows, state = [], np.zeros(4)
    for f in frequencies:
        row, sol = run_point(calibrated_parameters(force, f, circuit, target), state, cycles, samples)
        rows.append(row)
        state = sol.y[:, -1]
    name = f'{force:g}N_{circuit}_{direction}'
    np.savetxt(OUT/f'{name}.csv', rows, delimiter=',', header=','.join(FIELDS), comments='')
    return name, rows


def make_plots(results, target):
    f = np.linspace(4, 6.5, 1201)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(f, locked_amplitude(f, paper_parameters())*1000, label='Published stiffness')
    ax.plot(f, locked_amplitude(f, calibrated_parameters(target_hz=target))*1000, label='Calibrated stiffness')
    ax.axvline(target, color='gray', ls=':', label='Calibration target ≈4.85 Hz')
    ax.set(xlabel='Frequency (Hz)', ylabel='Locked amplitude (mm)', title='25 N locked reference before and after calibration')
    ax.legend(); ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(OUT/'baseline_comparison.png', dpi=160); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    f = np.linspace(4, 6, 1001)
    for ax, force in zip(axes, (5, 25)):
        ax.plot(f, locked_amplitude(f, calibrated_parameters(force, target_hz=target))*1000, 'k--', label='Locked reference')
        for name, rows in results:
            if name.startswith(f'{force}N_'):
                a = np.asarray(rows)
                ax.plot(a[:, 0], a[:, 1]*1000, '.-', ms=3, label=name.split('N_')[1].replace('_', ' '))
        ax.set(xlabel='Frequency (Hz)', ylabel='Sampled peak |x1| (mm)', title=f'{force} N with inferred stiffness')
        ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.suptitle('Calibrated model: masses unchanged; coupling reconstructed from Figure 8')
    fig.tight_layout(); fig.savefig(OUT/'frequency_response.png', dpi=160); plt.close(fig)


def validate(target, cycles, samples):
    checks = []
    for circuit in ('open', 'closed'):
        for f in (4.7, 4.85, 5.1):
            p = calibrated_parameters(25., f, circuit, target)
            base, sol = run_point(p, np.zeros(4), cycles, samples)
            fine, _ = run_point(p, np.zeros(4), cycles, samples*2, 1e-10)
            long, _ = run_point(p, np.zeros(4), cycles*2, samples)
            checks.append(dict(circuit=circuit, frequency_hz=f, peak_m=base[1],
                               resolution_change=abs(fine[1]-base[1])/base[1],
                               duration_change=abs(long[1]-base[1])/base[1],
                               late_window_change=base[4], fit_range_ok=bool(base[-1])))
            if circuit == 'closed' and f == 4.85:
                fig, axes = plt.subplots(2, 1, figsize=(8, 5.5))
                for ax in axes:
                    ax.plot(sol.t, sol.y[0]*1000, label='Primary mass')
                    ax.plot(sol.t, sol.y[2]*1000, alpha=.7, label='NES')
                    ax.set(xlabel='Time (s)', ylabel='Displacement (mm)')
                axes[0].axvline(sol.t[-1]/2, color='gray', ls='--')
                axes[0].set_title('25 N closed circuit at 4.85 Hz from rest')
                axes[0].legend(); axes[1].set_xlim(sol.t[-1]-3/f, sol.t[-1])
                axes[1].set_title('Final three forcing periods')
                fig.tight_layout(); fig.savefig(OUT/'time_response.png', dpi=160); plt.close(fig)
    return checks


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--target-hz', type=float, default=4.85)
    ap.add_argument('--cycles', type=int, default=160)
    ap.add_argument('--samples', type=int, default=80)
    ap.add_argument('--points', type=int, default=41)
    args = ap.parse_args()
    if args.cycles < 8 or args.cycles % 4 or args.samples < 20 or args.points < 3 or args.workers < 1:
        ap.error('Require cycles divisible by 4 and >=8, samples>=20, points>=3, workers>=1')
    p = calibrated_parameters(target_hz=args.target_hz)
    OUT.mkdir(parents=True, exist_ok=True)
    cases = [(5, 'open', 'up')] + [(25, c, d) for c in ('open', 'closed') for d in ('up', 'down')]
    tasks = [(*case, args.target_hz, args.cycles, args.samples, args.points) for case in cases]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(sweep, tasks))
    make_plots(results, args.target_hz)
    checks = validate(args.target_hz, args.cycles, args.samples)
    summary = dict(status='Inferred stiffness calibration; not a confirmed correction to the paper',
                   settings=vars(args), calibrated_k1_N_per_m=p.k1,
                   published_k1_N_per_m=paper_parameters().k1,
                   locked_peak_hz_m=locked_peak(p),
                   sensitivity_k1_N_per_m={str(f): calibrated_parameters(target_hz=f).k1 for f in (4.80, 4.85, 4.90)},
                   python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                   comparisons={}, checks=checks)
    for name, rows in results:
        a = np.asarray(rows); peak=a[np.argmax(a[:, 1])]
        summary['comparisons'][name] = dict(max_peak_m=float(peak[1]), frequency_hz=float(peak[0]),
            window_flags=int(np.sum(a[:, 4]>.01)), range_flags=int(np.sum(a[:, -1]==0)))
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
