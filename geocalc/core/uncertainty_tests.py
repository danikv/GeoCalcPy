import unittest
import uncertainty
import numpy.testing


class UncertaintyConversionTests(unittest.TestCase):
    def test_convert_uncertainty_matrix_to_ellipse(self):
        matrix = [[2, 2], [2, 4]]
        major, minor, angle = uncertainty.convert_uncertainty_matrix_to_ellipse(matrix)
        self.assertAlmostEqual(angle, 58.28252558, places=1)
        self.assertAlmostEqual(major, 6.864736833, places=1)
        self.assertAlmostEqual(minor, 2.622096146, places=1)

    def test_convert_uncertainty_ellipse_to_matrix(self):
        matrix = uncertainty.convert_uncertainty_ellipse_to_matrix(6.864736833, 2.622096146, 58.28252558)
        numpy.testing.assert_array_almost_equal(matrix, [[2, 2], [2, 4]])
