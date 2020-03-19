import unittest

from geocalc.core import position


class CoordinatesConverterTests(unittest.TestCase):
    def test_lla2ecef(self):
        lat, lon, alt = 30, 31, 100
        ecef = position.lla2ecef(lat, lon, alt)
        x, y, z = ecef[0], ecef[1], ecef[2]
        self.assertAlmostEqual(x, 4738715.05395, places=4)
        self.assertAlmostEqual(y, 2847307.26071, places=4)
        self.assertAlmostEqual(z, 3170423.73533, places=4)

    def test_ecef2lla(self):
        x, y, z = 4738715.05, 2847307.26, 3170423.73
        lat, lon, alt = position.ecef2lla([x, y, z])
        self.assertAlmostEqual(lat, 30, places=1)
        self.assertAlmostEqual(lon, 31, places=1)
        self.assertAlmostEqual(alt, 100, places=1)

    def test_ecef2ned(self):
        lat_ref, lon_ref, alt_ref = 29, 30, 0
        x, y, z = 4738715.05, 2847307.26, 3170423.73
        ned = position.ecef2ned([x, y, z], lat_ref, lon_ref, alt_ref)
        north, east, down = ned[0], ned[1], ned[2]
        self.assertAlmostEqual(north, 93129.29, places=1)
        self.assertAlmostEqual(east, 96482.89, places=1)
        self.assertAlmostEqual(down, -6371513.44, places=1)

    def test_ned2ecef(self):
        lat_ref, lon_ref, alt_ref = 29, 30, 0
        north, east, down = 93129.29, 96482.89, -6371513.44
        ecef = position.ned2ecef([north, east, down], lat_ref, lon_ref, alt_ref)
        x, y, z = ecef[0], ecef[1], ecef[2]
        self.assertAlmostEqual(x, 4738715.05, places=1)
        self.assertAlmostEqual(y, 2847307.26, places=1)
        self.assertAlmostEqual(z, 3170423.73, places=1)

    def test_ned2lla(self):
        lat_ref, lon_ref, alt_ref = 29, 30, 0
        north, east, down = 93129.29, 96482.89, -6371513.44
        lat, lon, alt = position.ned2lla([north, east, down], lat_ref, lon_ref, alt_ref)
        self.assertAlmostEqual(lat, 30, places=1)
        self.assertAlmostEqual(lon, 31, places=1)
        self.assertAlmostEqual(alt, 100, places=1)

    def test_lla2ned(self):
        lat_ref, lon_ref, alt_ref = 29, 30, 0
        lat, lon, alt = 30, 31, 100
        ned = position.lla2ned(lat, lon, alt, lat_ref, lon_ref, alt_ref)
        north, east, down = ned[0], ned[1], ned[2]
        self.assertAlmostEqual(north, 93129.29, places=1)
        self.assertAlmostEqual(east, 96482.89, places=1)
        self.assertAlmostEqual(down, -6371513.44, places=1)

    def test_ecef2polar(self):
        x, y, z = 4738715.05395, 2847307.26071, 3170423.73533
        x_ref, y_ref, z_ref = 4834879.70745, 2791419.10059, 3073901.20054
        r, bearing, azimuth = position.ecef2polar([x, y, z], [x_ref, y_ref, z_ref])
        self.assertAlmostEqual(r, 147267.53459829, places=7)
        self.assertAlmostEqual(bearing, 300.16393546, places=7)
        self.assertAlmostEqual(azimuth, 49.04826534, places=7)

    def test_polar2ecef(self):
        r, bearing, azimuth = 147267.53459829, 300.16393546, 49.04826534
        x_ref, y_ref, z_ref = 4834879.70745, 2791419.10059, 3073901.20054
        ecef = position.polar2ecef(r, bearing, azimuth, [x_ref, y_ref, z_ref])
        x, y, z = ecef[0], ecef[1], ecef[2]
        self.assertAlmostEqual(x, 4738715.05395, places=4)
        self.assertAlmostEqual(y, 2847307.26071, places=4)
        self.assertAlmostEqual(z, 3170423.73533, places=4)


if __name__ == '__main__':
    unittest.main()
