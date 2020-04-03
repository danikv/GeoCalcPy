import numpy as np
import math

from geocalc.core import constants
from geocalc.core.utils import input_check_Nx1 as _input_check_Nx1
from geocalc.core.utils import input_check_Nx3 as _input_check_Nx3
from geocalc.core.utils import positive_angle

def radiansChecker(deg, angle_in_radians = False):
    if not angle_in_radians:
        deg = np.radians(deg)
    return deg

def lla2ecef(lat, lon, alt, angle_in_radians = False):
    lat, lon = radiansChecker([lat, lon], angle_in_radians)

    cos_lat = math.cos(lat)
    sin_lat = math.sin(lat)
    cos_lon = math.cos(lon)
    sin_lon = math.sin(lon)

    R = math.pow(constants.a,2)/math.sqrt(math.pow(constants.a,2)*math.pow(cos_lat,2)+
                                    math.pow(constants.b,2)*math.pow(sin_lat,2))
    X = (R+alt)*cos_lat*cos_lon
    Y = (R+alt)*cos_lat*sin_lon
    Z = (math.pow(constants.b,2)*R/math.pow(constants.a,2)+alt)*sin_lat
    return X, Y, Z


def ecef2lla(ecef):
    """
    Calculate the Latitude, Longitude and Altitude of a point located on earth 
    given the ECEF Coordinates.
    
    References
    ----------
    .. [1] Jekeli, C.,"Inertial Navigation Systems With Geodetic
       Applications", Walter de Gruyter, New York, 2001, pp. 24
    
    Parameters
    ----------
    ecef : {(N,3)} array like input of ECEF coordinate in X, Y, and Z column, unit is meters
    latlon_unit : {('deg','rad')} specifies the output latitude and longitude unit
    
    Returns
    -------
    lat : {(N,)} array like latitude in unit specified by latlon_unit
    lon : {(N,)} array like longitude in unit specified by latlon_unit
    alt : {(N,)} array like altitude in meters
    """
    ecef, N = _input_check_Nx3(ecef)
    if (N > 1):
        x = ecef[:, 0];
        y = ecef[:, 1];
        z = ecef[:, 2]
    else:
        x = ecef[0];
        y = ecef[1];
        z = ecef[2]

    p = np.sqrt(x ** 2 + y ** 2)

    theta = np.arctan2((z * constants.a), (p * constants.b))

    lon = np.arctan2(y, x)

    lat = np.arctan2((z + (constants.e2 ** 2) * constants.b * (np.sin(theta) ** 3)),
                     ((p - (constants.e_sqrd) * constants.a * (np.cos(theta) ** 3))))
    N = constants.a / (np.sqrt(1 - ((constants.e_sqrd) * (np.sin(lat) ** 2))))

    m = (p / np.cos(lat))
    alt = m - N

    lat = np.rad2deg(lat)
    lon = np.rad2deg(lon)

    return lat, lon, alt


def lla2ned(lat, lon, alt, lat_ref, lon_ref, alt_ref, latlon_unit='deg', angle_in_radians = False):
    """
    Convert Latitude, Longitude, Altitude to its resolution in the NED
    coordinate. The center of the NED coordiante is given by lat_ref, lon_ref,
    and alt_ref.
    
    For example, this can be used to convert GPS data to a local NED frame.
    
    Parameters
    ----------
    lat : {(N,)} array like latitude, unit specified by latlon_unit, default in deg
    lon : {(N,)} array like longitude, unit specified by latlon_unit, default in deg
    alt : {(N,)} array like altitude, unit specified by alt_unit, default in m
    
    lat_ref : Reference latitude, unit specified by latlon_unit, default in deg
    lon_ref : Reference longitude, unit specified by latlon_unit, default in deg
    alt : Reference altitude, unit specified by alt_unit, default in m
    
    Returns
    -------
    ned : {(N,3)} array like ecef position, unit is the same as alt_unit        
    """
    lat, lon, lat_ref, lon_ref = radiansChecker([lat, lon, lat_ref, lon_ref], angle_in_radians)
    ecef = lla2ecef(lat, lon, alt, True)
    ned = ecef2ned(ecef, lat_ref, lon_ref, alt_ref, True)
    return ned


def ned2lla(ned, lat_ref, lon_ref, alt_ref, angle_in_radians = False):
    """
    Calculate the Latitude, Longitude and Altitude of points given by NED coordinates
    where NED origin given by lat_ref, lon_ref, and alt_ref.
    Parameters
    ----------
    ned : {(N,3)} array like input of NED coordinate in N, E, and D column, unit is meters
    lat_ref : Reference latitude, unit specified by latlon_unit, default in deg
    lon_ref : Reference longitude, unit specified by latlon_unit, default in deg
    alt_ref : Reference altitude, unit specified by alt_unit, default in m
    latlon_unit : {('deg','rad')} specifies the output latitude and longitude unit
    
    Returns
    -------
    lat : {(N,)} array like latitude in unit specified by latlon_unit
    lon : {(N,)} array like longitude in unit specified by latlon_unit
    alt : {(N,)} array like altitude in meters
    Note
    ----
    This method is a wrapper on ned2ecef (add ecef of NED-origin) and ecef2lla.
    """

    lat_ref, lon_ref = radiansChecker([lat_ref, lon_ref], angle_in_radians)
    ecef = ned2ecef(ned, lat_ref, lon_ref, alt_ref, True)
    lla = ecef2lla(ecef) #In deg?
    return lla


def ned2ecef(ned, lat_ref, lon_ref, alt_ref, angle_in_radians = False):
    """
    Transform a vector resolved in NED (origin given by lat_ref, lon_ref, and alt_ref)
    coordinates to its ECEF representation. 
    Parameters
    ----------
    ned : {(N,3)} input array, units of meters
    lat_ref : Reference latitude, unit specified by latlon_unit, default in deg
    lon_ref : Reference longitude, unit specified by latlon_unit, default in deg
    alt_ref : Reference altitude, unit specified by alt_unit, default in m
    
    Returns
    -------
    ecef : {(N,3)} array like ned vector, in the ECEF frame, units of meters
    Notes
    -----
    The NED vector is treated as a relative vector, and hence the ECEF representation
    returned is NOT converted into an absolute coordinate.  This means that the 
    magnitude of `ned` and `ecef` will be the same (bar numerical differences).
    
    Examples
    --------
    >>> import navpy
    >>> ned = [0, 0, 1]
    >>> lat_ref, lon_ref, alt_ref = 45.0, -93.0, 250.0 # deg, meters
    >>> ecef = navpy.ned2ecef(ned, lat_ref, lon_ref, alt_ref)
    >>> print("NED:", ned)
    >>> print("ECEF:", ecef)
    >>> print("Notice that 'down' is not same as 'ecef-z' coordinate.")
    """
    lat_ref, N1 = _input_check_Nx1(lat_ref)
    lon_ref, N2 = _input_check_Nx1(lon_ref)
    alt_ref, N3 = _input_check_Nx1(alt_ref)
    lat_ref, lon_ref = radiansChecker([lat_ref, lon_ref], angle_in_radians)

    if ((N1 != 1) or (N2 != 1) or (N3 != 1)):
        raise ValueError('Reference Location can only be 1')

    ned, N = _input_check_Nx3(ned)

    ned = ned.T

    C = np.zeros((3, 3))

    C[0, 0] = -np.sin(lat_ref) * np.cos(lon_ref)
    C[0, 1] = -np.sin(lat_ref) * np.sin(lon_ref)
    C[0, 2] = np.cos(lat_ref)

    C[1, 0] = -np.sin(lon_ref)
    C[1, 1] = np.cos(lon_ref)
    C[1, 2] = 0

    C[2, 0] = -np.cos(lat_ref) * np.cos(lon_ref)
    C[2, 1] = -np.cos(lat_ref) * np.sin(lon_ref)
    C[2, 2] = -np.sin(lat_ref)

    # C defines transoformation: ned = C * ecef.  Hence used transpose.
    ecef = np.dot(C.T, ned)
    ecef = ecef.T

    if (N == 1):
        ecef = ecef.reshape(3)

    return ecef


def ecef2ned(ecef, lat_ref, lon_ref, alt_ref, angle_in_radians = False):
    """
    Transform a vector resolved in ECEF coordinate to its resolution in the NED
    coordinate. The center of the NED coordiante is given by lat_ref, lon_ref,
    and alt_ref.
    
    Parameters
    ----------
    ecef : {(N,3)} input vector expressed in the ECEF frame
    lat_ref : Reference latitude, unit specified by latlon_unit, default in deg
    lon_ref : Reference longitude, unit specified by latlon_unit, default in deg
    alt : Reference altitude, unit specified by alt_unit, default in m
    
    Returns
    -------
    ned : {(N,3)} array like ecef position, unit is the same as alt_unit
    """

    lat_ref, N1 = _input_check_Nx1(lat_ref)
    lon_ref, N2 = _input_check_Nx1(lon_ref)
    alt_ref, N3 = _input_check_Nx1(alt_ref)
    lat_ref, lon_ref = radiansChecker([lat_ref, lon_ref], angle_in_radians)

    if ((N1 != 1) or (N2 != 1) or (N3 != 1)):
        raise ValueError('Reference Location can only be 1')

    ecef, N = _input_check_Nx3(ecef)

    ecef = ecef.T

    C = np.zeros((3, 3))

    C[0, 0] = -np.sin(lat_ref) * np.cos(lon_ref)
    C[0, 1] = -np.sin(lat_ref) * np.sin(lon_ref)
    C[0, 2] = np.cos(lat_ref)

    C[1, 0] = -np.sin(lon_ref)
    C[1, 1] = np.cos(lon_ref)
    C[1, 2] = 0

    C[2, 0] = -np.cos(lat_ref) * np.cos(lon_ref)
    C[2, 1] = -np.cos(lat_ref) * np.sin(lon_ref)
    C[2, 2] = -np.sin(lat_ref)

    ned = np.dot(C, ecef)
    ned = ned.T

    if (N == 1):
        ned = ned.reshape(3)

    return ned


def ecef2polar(ecef, ecef_ref):
    """
    Transform a vector resolved in ECEF coordinate to its resolution in the Polar
    coordinate.
    
    Parameters
    ----------
    ecef : {(N,3)} input vector expressed in the ECEF frame
    ecef_ref : {(N,3)} reference vector expressed in the ECEF frame
    
    Returns
    -------
    polar : range, bearing, elevation
    """
    ecef = np.array(ecef) - np.array(ecef_ref)
    ecef, N = _input_check_Nx3(ecef)
    if (N > 1):
        x = ecef[:, 0];
        y = ecef[:, 1];
        z = ecef[:, 2]
    else:
        x = ecef[0];
        y = ecef[1];
        z = ecef[2]

    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    elevation = np.rad2deg(np.arctan2(np.sqrt(x ** 2 + y ** 2), z))
    bearing = positive_angle(np.rad2deg(constants.half_pi - np.arctan2(y, x)))
    return r, bearing, elevation


def lla2polar(lat, lon, alt, lat_ref, lon_ref, alt_ref, angle_in_radians = False):
    """
    Transform a vector resolved in LLA coordinate to its resolution in the Polar
    coordinate.
    
    Parameters
    ----------
    lla : lat, long, alt input coordiantes
    lla : lat, long, alt reference input coordinates
    
    Returns
    -------
    polar : range, bearing, elevation
    """

    lat, lon, lat_ref, lon_ref = radiansChecker([lat, lon, lat_ref, lon_ref], angle_in_radians)
    ecef = lla2ecef(lat, lon, alt, True)
    ecef_ref = lla2ecef(lat_ref, lon_ref, alt_ref, True)
    return ecef2polar(ecef, ecef_ref)


def polar2ecef(r, bearing, elevation, ecef_ref, angle_in_radians = False):
    """
    Transform a vector resolved in ECEF coordinate to its resolution in the Polar
    coordinate.
    
    Parameters
    ----------
    range, bearing, elevation : {(N,3)} input vector expressed in the Polar coordinates
    ecef_ref : {(N,3)} reference vector expressed in the ECEF frame
    
    Returns
    -------
    polar : range, bearing, elevation
    """
    ecef, N = _input_check_Nx3(ecef_ref)
    if (N > 1):
        x = ecef[:, 0];
        y = ecef[:, 1];
        z = ecef[:, 2]
    else:
        x = ecef[0];
        y = ecef[1];
        z = ecef[2]

    bearing, elevation = radiansChecker([bearing, elevation], angle_in_radians)
    bearing, elevation = constants.half_pi - bearing, elevation
    x = r * np.cos(bearing) * np.sin(elevation)
    y = r * np.sin(bearing) * np.sin(elevation)
    z = r * np.cos(elevation)
    return np.array([x, y, z]) + np.array(ecef_ref)


def polar2lla(r, bearing, elevation, lat_ref, lon_ref, alt_ref, angle_in_radians = False):
    """
    Transform a vector resolved in ECEF coordinate to its resolution in the Polar
    coordinate.
    
    Parameters
    ----------
    range, bearing, elevation : {(N,3)} input vector expressed in the Polar coordinates
    lat, lon, alt reference : {(N,3)} reference vector expressed in the LLA coordinates
    
    Returns
    -------
    polar : range, bearing, elevation
    """

    bearing, elevation, lat_ref, lon_ref = radiansChecker([bearing, elevation, lat_ref, lon_ref], angle_in_radians)
    ecef_ref = lla2ecef(lat_ref, lon_ref, alt_ref, True)
    ecef = polar2ecef(r, bearing, elevation, ecef_ref, True)
    return ecef2lla(ecef)


def enu2polar(x_east, y_north, z_up):
    projected_range = math.sqrt(x_east * x_east + y_north * y_north)
    range = math.sqrt(projected_range*projected_range + z_up * z_up)
    bearing = math.acos(y_north / projected_range)
    elevation = math.atan(z_up, projected_range)
    return range, bearing, elevation


def polar2enu(range, bearing, elevation, angle_in_radians = False):
    bearing, elevation = radiansChecker([bearing, elevation], angle_in_radians)
    projected_range = range * math.cos(elevation)
    x_east = projected_range * math.sin(bearing)
    y_north = projected_range * math.cos(bearing)
    z_up = range * math.sin(elevation)
    return x_east, y_north, z_up


def polar_height2elevation(ownship_lat, ownship_long, ownship_height, range, bearing, elevation, angle_in_radians = False):
    ownship_lat, ownship_long, bearing, elevation = radiansChecker([ownship_lat, ownship_long, bearing, elevation], angle_in_radians)
    ENU_x, ENU_y, ENU_z = polar2enu(range, bearing, elevation, True)
    ECEF_x, ECEF_y, ECEF_z = enu2ecef(ENU_x, ENU_y, ENU_z, ownship_lat, ownship_long, ownship_height, True)
    lat, lon, alt = ecef2lla(ECEF_x, ECEF_y, ECEF_z) #In deg?
    return alt


def enu2ecef(x, y, z, lat_ENU, long_ENU, h_ENU, angle_in_radians = False):
    lat_ENU, long_ENU = radiansChecker([lat_ENU, long_ENU], angle_in_radians)
    transMat = [[-math.sin(long_ENU), -math.sin(lat_ENU)*math.cos(long_ENU), math.cos(lat_ENU)*math.cos(long_ENU)],
                [math.cos(long_ENU), -math.sin(lat_ENU)*math.sin(long_ENU), math.cos(lat_ENU)*math.sin(long_ENU)],
                [0, math.cos(lat_ENU), math.sin(lat_ENU)]]
    transCoordinates = np.dot(transMat, [x,y,z])
    # print(transCoordinates)
    relativeCoordinates = lla2ecef(lat_ENU, long_ENU, h_ENU, True)
    # print(relativeCoordinates)
    Coordinates_new = np.add(transCoordinates, relativeCoordinates)
    return Coordinates_new
