import math
import numpy as np


def convert_uncertainty_matrix_to_ellipse(uncertainty, sigma=3):
    """
    Convert Uncertainty matrix to ellipse

    Parameters
    ----------
    uncertainty : {(2,2)} array like unceratinty matrix, must be symmetric positive definite
    sigma : float, the sigma value of the ellipse, responsible for size of the ellipse

    Returns
    -------
    major : float stating the major radius of the ellipse
    minor : float stating the minor radius of the ellipse
    angle : float stating the rotation angle of the ellipse w.r.t the X axis, in degrees
    """
    if (not isinstance(uncertainty, list)) or np.array(uncertainty).shape != (2, 2):
        raise ValueError('Not a 2 x 2 array')

    uncertainty = np.array(uncertainty)

    if not np.allclose(uncertainty, uncertainty.T, rtol=1e-5, atol=1e-8):
        raise ValueError('Matrix is not symmetric')

    eig_vals, eig_vecs = np.linalg.eigh(uncertainty)

    if not np.all(eig_vals > 0):
        raise ValueError('Matrix is not positive definite')

    max_eig_vec = eig_vecs[:, 1]
    angle = np.degrees(np.arctan2(max_eig_vec[1], max_eig_vec[0]))
    minor, major = sigma * np.sqrt(eig_vals)

    return major, minor, angle


def convert_uncertainty_ellipse_to_matrix(major, minor, angle, sigma=3):
    """
    Convert Uncertainty ellipse to matrix

    Parameters
    ----------
    major : float, major axis of the ellipse
    minor : float, minor axis of the ellipse
    angle : float stating the rotation angle of the ellipse w.r.t the X axis, in degrees
    sigma : float, the sigma value of the ellipse, responsible for size of the ellipse

    Returns
    -------
    uncertainty : {(2,2)} array like unceratinty matrix, must be symmetric positive definite
    """
    if not (0 < angle < 360):
        raise ValueError('Angle is not valid')

    if minor < 0 or major < 0 or major < minor:
        raise ValueError('Illegal major or minor axis')

    major /= sigma
    minor /= sigma
    angle_in_rad = math.radians(angle)
    var_x = (major ** 2) * (math.cos(angle_in_rad) ** 2) + (minor ** 2) * (math.sin(angle_in_rad) ** 2)
    var_y = (major ** 2) * (math.sin(angle_in_rad) ** 2) + (minor ** 2) * (math.cos(angle_in_rad) ** 2)
    covar_xy = (major ** 2 - minor ** 2) * math.sin(angle_in_rad) * math.cos(angle_in_rad)

    return np.array([[var_x, covar_xy], [covar_xy, var_y]])


def uncertainty3_spherical2cartesian(uncertainty):
    P = 0
    return P


def uncertainty3_cartesian2spherical(uncertainty):
    P = 0
    return P


def J3_spherical2cartesian(theta, r, phi):
    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    J3_11 = r*cos_theta*cos_phi
    J3_12 = sin_theta*cos_phi
    J3_13 = -r*sin_theta*sin_phi
    J3_21 = -r*sin_theta*cos_phi
    J3_22 = cos_theta*cos_phi
    J3_23 = -r*cos_theta*sin_phi
    J3_31 = 0
    J3_32 = sin_phi
    J3_33 = r*cos_phi

    J3 = [[J3_11,J3_12,J3_13], [J3_21,J3_22,J3_23], [J3_31,J3_32, J3_33]]
    return J3


def J3_cartesian2spherical(x,y,z):
    R = math.sqrt(math.pow(x,2)+math.pow(y,2)+math.pow(z,2))
    r = math.sqrt(math.pow(x,2)+math.pow(y,2))

    J3_11 = y/math.pow(r,2)
    J3_12 = x/math.pow(r,2)
    J3_13 = 0
    J3_21 = x/R
    J3_22 = y/R
    J3_23 = z/R
    J3_31 = (x*z)/(math.pow(R,2)*r)
    J3_32 = (y*z)/(math.pow(R,2)*r)
    J3_33 = r/math.pow(R,2)

    J3 = [[J3_11,J3_12,J3_13], [J3_21,J3_22,J3_23], [J3_31,J3_32, J3_33]]
    return J3
