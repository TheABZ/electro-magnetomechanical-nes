"""Generate a concise report from results/summary.json; called by reproduce.py."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def write_report(summary):
    p = summary['locked_25N_peak_hz_m']
    lines = ['# Numerical findings', '',
             'These are simulations using the reported table values and a digitized coupling fit. '
             'They are **not a validated reproduction of the published figures**.', '',
             f'Analytic locked-NES reference at 25 N: {p[1]*1000:.4f} mm at {p[0]:.4f} Hz. '
             'The paper describes about 8 mm and plots the peak near 4.85 Hz.', '',
             '## Sweeps', '',
             '| Run | Maximum sampled primary peak (mm) | Frequency (Hz) | Points with late-window peak change >1% | Points outside coupling-fit range |',
             '|---|---:|---:|---:|---:|']
    for name, row in summary['comparisons'].items():
        lines.append(f"| {name} | {row['max_peak_m']*1000:.4f} | {row['frequency_hz']:.2f} | "
                     f"{row['points_with_window_change_over_1pct']} | {row['points_outside_fit_range']} |")
    lines += ['', 'The maxima are over the sampled frequency grid and observation windows, not '
              'optimized resonance peaks. Peak amplitude and forcing-frequency Fourier amplitude '
              'are both saved in the CSVs. The open-circuit curve depends on continuation direction.', '',
              '## Selected convergence checks (independent runs from rest)', '',
              '| Circuit | Frequency (Hz) | Peak (mm) | Resolution/tolerance change (%) | Doubled-duration change (%) |',
              '|---|---:|---:|---:|---:|']
    for row in summary['convergence_checks']:
        lines.append(f"| {row['circuit']} | {row['frequency_hz']:.2f} | {row['peak_m']*1000:.4f} | "
                     f"{row['resolution_relative_change']*100:.5f} | {row['duration_relative_change']*100:.5f} |")
    max_resolution = max(r['resolution_relative_change'] for r in summary['convergence_checks'])
    max_duration = max(r['duration_relative_change'] for r in summary['convergence_checks'])
    lines += ['', f'Maximum selected resolution/tolerance change: {max_resolution*100:.5f}%. '
              f'Maximum selected doubled-duration change: {max_duration*100:.5f}%.', '',
              'These checks quantify sensitivity at six points; they do not certify every sweep '
              'point or establish a unique attractor. Larger changes can reflect long transients '
              'or modulation as well as numerical sensitivity.', '',
              '## Figures', '', '![Frequency response](frequency_response.png)', '',
              '![Coupling reconstruction](coupling_fit.png)', '',
              '![Time response](time_response.png)', '',
              '## Interpretation', '',
              'The numerical workflow is operational, reproducible, and tested. The reported '
              'mass/stiffness mismatch must be resolved before claiming agreement with Figures '
              '11, 16 or 17. No parameter was tuned to conceal that discrepancy. The reconstructed '
              'coupling coefficients also remain approximate rather than author supplied.', '',
              'See [REPLICATION.md](../REPLICATION.md) for methods, units, provenance, '
              'polynomial validity limits, and remaining scientific questions.', '']
    (ROOT/'results/RESULTS.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    write_report(json.loads((ROOT/'results/summary.json').read_text()))
