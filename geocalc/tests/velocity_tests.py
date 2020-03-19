import unittest

from geocalc.core import velocity


class MyTestCase(unittest.TestCase):
    def test_cartesian2directional(self):
        vx, vy, vz = 21, 12.2, 100
        speed, course, vz = velocity.cartesian2directional(vx, vy, vz)
        self.assertAlmostEqual(speed, 24.3, places=1)
        self.assertAlmostEqual(course, 30, places=0)
        self.assertAlmostEqual(vz, 100, places=1)


if __name__ == '__main__':
    unittest.main()
