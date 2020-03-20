import unittest

from geocalc.core import velocity


class MyTestCase(unittest.TestCase):
    def test_cartesian2directional(self):
        vx, vy, vz = 21, 12.2, 100
        speed, course, vz = velocity.cartesian2directional(vx, vy, vz)
        self.assertAlmostEqual(speed, 24.3, places=1)
        self.assertAlmostEqual(course, 30, places=0)
        self.assertAlmostEqual(vz, 100, places=1)

    def test_cartesian2directional2(self):
        vx, vy, vz = 14.6, -8.62, 100
        speed, course, vz = velocity.cartesian2directional(vx, vy, vz)
        self.assertAlmostEqual(speed, 17, places=1)
        self.assertAlmostEqual(course, -30.6, places=1)
        self.assertAlmostEqual(vz, 100, places=1)

    def test_cartesian2directional_vx_is_zero_vy_is_negative(self):
        vx, vy, vz = 0, -8.62, 100
        speed, course, vz = velocity.cartesian2directional(vx, vy, vz)
        self.assertEqual(speed, 8.62)
        self.assertEqual(course, -90)
        self.assertEqual(vz, 100)

    def test_cartesian2directional_vx_is_zero_vy_is_positive(self):
        vx, vy, vz = 0, 8.62, 100
        speed, course, vz = velocity.cartesian2directional(vx, vy, vz)
        self.assertEqual(speed, 8.62)
        self.assertEqual(course, 90)
        self.assertEqual(vz, 100)

    def test_directional2cartesian(self):
        speed, course, vz = 24.3, 30, 100
        vx, vy, vz = velocity.directional2cartesian(speed, course, vz)
        self.assertAlmostEqual(vx, 21, places=1)
        self.assertAlmostEqual(vy, 12.2, places=0)
        self.assertAlmostEqual(vz, 100, places=1)


if __name__ == '__main__':
    unittest.main()
