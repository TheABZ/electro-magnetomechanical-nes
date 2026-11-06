"""Physics and numerical regression checks; run with unittest discovery."""
import unittest
from dataclasses import replace

import numpy as np
from numpy.testing import assert_allclose

from nes import Parameters, peak_amplitude, rhs, simulate, theta


class ModelTests(unittest.TestCase):
    def setUp(self):
        self.p = Parameters(2., 0.5, 0.3, 0.2, 4., 0.7, 1.2,
                            2., 3., 0.8, 1.7, (0.4, -0.2, 0.1))

    def test_original_equation_residuals(self):
        p = self.p
        t = 0.37
        x1, v1, x2, v2 = state = [0.2, -0.3, -0.1, 0.4]
        _, a1, _, a2 = rhs(t, state, p)
        em = (0.4 - 0.2*(x1-x2) + 0.1*(x1-x2)**2)**2 / 5
        r1 = (p.m1*a1 + p.c1*v1 + p.k1*x1 + p.c2*(v1-v2)
              + p.k2l*(x1-x2) + p.k2c*(x1-x2)**3
              - p.Fe*np.sin(p.omega*t))
        r2 = (p.m2*a2 + (p.c2+em)*(v2-v1) + p.k2l*(x2-x1)
              + p.k2c*(x2-x1)**3)
        assert_allclose([r1, r2], 0, atol=1e-14)

    def test_uncoupled_exact_solution(self):
        p = replace(self.p, c1=0., c2=0., k2l=0., k2c=0., Fe=0.,
                    theta_coeffs=(0.,))
        t = np.linspace(0, 8, 401)
        sol = simulate(p, t, (1., 0., 0., 0.))
        assert_allclose(sol.y[0], np.cos(np.sqrt(p.k1/p.m1)*t), atol=1e-8)
        assert_allclose(sol.y[2:], 0, atol=1e-12)

    def test_linear_forced_amplitude(self):
        p = replace(self.p, c1=1., c2=0., k2l=0., k2c=0., theta_coeffs=(0.,))
        t = np.linspace(0, 100, 20001)
        sol = simulate(p, t)
        expected = abs(p.Fe) / np.hypot(p.k1-p.m1*p.omega**2, p.c1*p.omega)
        self.assertAlmostEqual(peak_amplitude(t, sol.y[0], transient_time=80),
                               expected, delta=expected*1e-5)

    def test_polynomial_and_open_circuit(self):
        self.assertEqual(theta(2., (3., 4., 5.)), 31.)
        assert_allclose(rhs(0.2, [0.2, 0.3, 0.1, -0.5], replace(self.p, Ri=np.inf)),
                        rhs(0.2, [0.2, 0.3, 0.1, -0.5], replace(self.p, theta_coeffs=(0.,))))

    def test_peak_cutoff_and_absolute_value(self):
        self.assertEqual(peak_amplitude([0, 1, 2, 3], [100, 50, -3, 2]), 3)
        self.assertEqual(peak_amplitude([0, 1, 2, 3], [100, 50, -3, 2],
                                        transient_time=1), 50)
        # Fraction applies to elapsed time, even with nonuniform samples.
        self.assertEqual(peak_amplitude([10, 10.1, 10.2, 18, 20],
                                       [99, 99, 99, -2, 1]), 2)

    def test_invalid_inputs(self):
        for t in ([0, 0], [1, 0], [0, np.nan], [0]):
            with self.assertRaises(ValueError):
                simulate(self.p, t)
        with self.assertRaises(ValueError):
            replace(self.p, m1=0)
        with self.assertRaises(ValueError):
            replace(self.p, Rc=0, Ri=0)
        with self.assertRaises(ValueError):
            peak_amplitude([0, 1], [0, 1], transient_time=1)


if __name__ == '__main__':
    unittest.main()
