"""Checks for paper parameter provenance and reconstructed coupling."""
import unittest
from dataclasses import replace

import numpy as np
from numpy.testing import assert_allclose

from nes import theta
from paper_model import coupling_fit, paper_parameters, locked_amplitude, locked_peak
from reproduce import run_point


class PaperTests(unittest.TestCase):
    def test_calibrated_resonance_and_unchanged_parameters(self):
        from paper_model import calibrated_parameters
        from dataclasses import asdict
        original, calibrated = paper_parameters(), calibrated_parameters(frequency_hz=6.)
        self.assertAlmostEqual(locked_peak(calibrated)[0], 4.85, places=12)
        self.assertAlmostEqual(calibrated.k1, 60804.45, delta=1.)
        for key, value in asdict(original).items():
            if key != 'k1':
                self.assertEqual(value, asdict(calibrated)[key])
        self.assertEqual(original.k1, 91324.)
        for invalid in (0., -1., float('nan')):
            with self.assertRaises(ValueError):
                calibrated_parameters(target_hz=invalid)

    def test_reported_values_and_circuits(self):
        p = paper_parameters()
        self.assertEqual((p.m1, p.m2, p.k1, p.Rc, p.Ri), (64., 1.4, 91324., 11.6, 0.))
        self.assertTrue(np.isinf(paper_parameters(circuit='open').Ri))
        with self.assertRaises(ValueError):
            paper_parameters(circuit='typo')

    def test_even_polynomial_reconstructs_digitized_curve(self):
        from paper_model import ROOT
        points = np.loadtxt(ROOT/'data/figure8_digitized.csv', delimiter=',', skiprows=1)
        coeffs = coupling_fit()['theta_coeffs']
        assert_allclose(theta(points[:, 0], coeffs), theta(-points[:, 0], coeffs))
        rmse = np.sqrt(np.mean((theta(points[:, 0], coeffs)-points[:, 1])**2))
        self.assertLess(rmse, .025)
        self.assertGreater(theta(0., coeffs), 12.)
        self.assertLess(theta(0., coeffs), 13.)

    def test_locked_baseline_has_total_mass(self):
        p = paper_parameters()
        f, peak = locked_peak(p)
        self.assertAlmostEqual(f, 5.945, delta=.01)
        self.assertGreater(peak, locked_amplitude(f-.1, p))
        self.assertGreater(peak, locked_amplitude(f+.1, p))
        self.assertEqual(locked_amplitude(0, p), p.Fe/p.k1)

    def test_fourier_metric_against_linear_analytic_response(self):
        p = replace(paper_parameters(force=1., frequency_hz=5.),
                    c2=0., k2l=0., k2c=0., theta_coeffs=(0.,))
        # Start exactly on the steady forced orbit to isolate Fourier integration.
        d = p.k1-p.m1*p.omega**2
        denominator = d*d+(p.c1*p.omega)**2
        initial = [-p.Fe*p.c1*p.omega/denominator, p.omega*p.Fe*d/denominator, 0., 0.]
        metrics, _ = run_point(p, initial, cycles=8, samples=80)
        expected = p.Fe/np.hypot(d, p.c1*p.omega)
        self.assertAlmostEqual(metrics[2], expected, delta=expected*1e-6)
        self.assertAlmostEqual(metrics[3], expected/np.sqrt(2), delta=expected*1e-6)


if __name__ == '__main__':
    unittest.main()
