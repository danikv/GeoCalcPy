"""
WGS 84 four defining parameters and several commonly used derived parameters.
All parameters are stored with the exact number of significant digits provided
by the WGS 84 rublished report.  Available parameters:

Parameters
----------

a : semi-major axis [m]
f : flattenning
omega_E : angular velocity of the Earth [rad/s]
GM : earth's gravitational constant [m^3/s^2]
     (note: for GPS applications, use GM_GPS)

_b : semi-minor axis [m]
_ecc : first eccentricity
_ecc_sqrd : first eccentricity squared

Copyright (c) 2014 NavPy Developers. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in
LICENSE.txt

References
----------
.. [1] NIMA Technical Report TR8350.2, "Department of Defense World Geodetic
       System 1984, Its Definition and Relationships With Local Geodetic Systems"

       http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
       Accessed on Nov. 19, 2013 

"""
import numpy as np

# Table 3.1: WGS 84  Four Defining Parameters
a = 6378137.0  # Semi-major Axis [m]
f = 1. / 298.257223563  # Flattening
omega_E = 7292115.0e-11  # Angular velocity of the Earth [rad/s]
omega_E_GPS = 7292115.1467e-11  # Angular velocity of the Earth [rad/s]
# According to ICD-GPS-200

GM = 3986004.418e8  # Earth's Gravitational Constant [m^3/s^2]
# (mass of earth's atmosphere included)

GM_GPS = 3986005.0e8  # The WGS 84 GM value recommended for GPS receiver usage
# by the GPS interface control document (ICD-GPS-200)
# differs from the current refined WGS 84 GM value.
#
# Details for this difference can be read in the WGS84
# reference: 3.2.3.2 "Special Considerations for GPS"

# Table 3.3: WGS 84 Ellipsoid Derived Geometric Constants
b = 6356752.3142  # Semi-minor axis [m]

f = 0.0034
e = np.sqrt(((a ** 2) - (b ** 2)) / (a ** 2))
e2 = np.sqrt(((a ** 2) - (b ** 2)) / (b ** 2))
e_sqrd = e ** 2
half_pi = np.pi / 2
half_pi_degree = 90.0
two_pi_degree = 360.0
