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
    return


def uncertainty3_cartesian2spherical(uncertainty):
    return


def J3_spherical2cartesian(theta, r, phi):
    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    J11 = r*cos_theta*cos_phi
    J12 = sin_theta*cos_phi
    J13 = -r*sin_theta*sin_phi
    J21 = -r*sin_theta*cos_phi
    J22 = cos_theta*cos_phi
    J23 = -r*cos_theta*sin_phi
    J31 = 0
    J32 = sin_phi
    J33 = r*cos_phi
    return


def J3_cartesian2spherical(uncertainty):
    J11 = 0
    J12 = 0
    J13 = 0
    J21 = 0
    J22 = 0
    J23 = 0
    J31 = 0
    J32 = 0
    J33 = 0
    return