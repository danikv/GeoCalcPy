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


def uncertainty3_spherical2cartesian(theta, r, phi, uncertainty):
    J = J3_spherical2cartesian(theta, r, phi)
    P = np.matmul(np.transpose(J), uncertainty)
    P = np.matmul(P, J)
    return P


def uncertainty3_cartesian2spherical(x, y, z, uncertainty):
    J = J3_cartesian2spherical(x, y, z)
    P = np.matmul(np.transpose(J), uncertainty)
    P = np.matmul(P, J)
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


def J3_cartesian2spherical(x, y, z):
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


def uncertainty6_spherical2cartesian(theta, r, phi, range_rate, bearing_rate, elevation_rate, uncertainty):
    J = J6_spherical2cartesian(theta, r, phi, range_rate, bearing_rate, elevation_rate)
    P = np.matmul(np.transpose(J), uncertainty)
    P = np.matmul(P, J)
    return P


def uncertainty6_cartesian2spherical(x, y, z, vx, vy, vz, uncertainty):
    J = J6_cartesian2spherical(x, y, z, vx, vy, vz)
    P = np.matmul(np.transpose(J), uncertainty)
    P = np.matmul(P, J)
    return P


def J6_spherical2cartesian(theta, r, phi, r_dot, theta_dot, phi_dot):
    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    
    J6 = np.zeros((6, 6))
    J6[1][1] = r*cos_theta*cos_phi
    J6[1][3] = sin_theta*cos_phi
    J6[1][5] = -r*sin_theta*sin_phi

    J6[2][1] = r_dot*cos_theta*cos_phi - theta_dot*r*sin_theta*cos_phi - phi_dot*r*cos_theta*sin_phi
    J6[2][2] = r*cos_theta*cos_phi
    J6[2][3] = -theta_dot*sin_theta*cos_phi - phi_dot*sin_theta*sin_phi
    J6[2][4] = sin_theta*cos_phi
    J6[2][5] = -r_dot*sin_theta*sin_phi - theta_dot*r*cos_theta*sin_phi - phi_dot*r*sin_theta*cos_phi
    J6[2][6] = -r*sin_theta*sin_phi

    J6[3][1] = -r*sin_theta*cos_phi
    J6[3][3] = sin_theta*cos_phi
    J6[3][5] = -r*cos_theta*sin_phi

    J6[4][1] = -r_dot*sin_theta*cos_phi - theta_dot*r*cos_theta*cos_phi + phi_dot*r*cos_theta*sin_phi
    J6[4][2] = -r*sin_theta*cos_phi
    J6[4][3] = -theta_dot*sin_theta*cos_phi - phi_dot*sin_theta*sin_phi
    J6[4][4] = cos_theta*cos_phi
    J6[4][5] = -r_dot*cos_theta*sin_phi + theta_dot*r*sin_theta*sin_phi - phi_dot*r*cos_theta*cos_phi
    J6[4][6] = -r*cos_theta*sin_phi

    J6[5][3] = sin_phi
    J6[5][5] = r*cos_phi

    J6[6][3] = phi_dot*cos_phi
    J6[6][4] = sin_phi
    J6[6][5] = r_dot*cos_phi - phi_dot*r*sin_phi
    J6[6][6] = r*cos_phi

    return J6


def J6_cartesian2spherical(x, y, z, vx, vy, vz):
    R = math.sqrt(math.pow(x,2)+math.pow(y,2)+math.pow(z,2))
    r = math.sqrt(math.pow(x,2)+math.pow(y,2))

    J6 = np.zeros((6, 6))
    J6[1][1] = y/(r**2)
    J6[1][3] = -x/(r**2)

    J6[2][1] = ((x**2-y**2)*vx-2*x*y*vy)/(r**4)
    J6[2][2] = y/(r**2)
    J6[2][3] = ((x**2-y**2)*vy-2*x*y*vx)/(r**4)
    J6[2][4] = -x/(r**2)

    J6[3][1] = x/R
    J6[3][3] = y/R
    J6[3][5] = z/R

    J6[4][1] = (vx*(y**2+z**2)-x*(y*vy+z*vz))/(R**3)
    J6[4][2] = x/R
    J6[4][3] = (vy*(x**2+z**2)-y*(x*vx+z*vz))/(R**3)
    J6[4][4] = y/R
    J6[4][5] = (vz*(y**2+x**2)-z*(y*vy+x*vx))/(R**3)
    J6[4][6] = z/R

    J6[5][1] = (x*z)/((R**2)*r)
    J6[5][3] = (y*z)/((R**2)*r)
    J6[5][5] = r/(R**2)

    J6[6][1] = -(vx*z*(-2*(x**4)-(x**2)*(y**2)+y**4)+vy*x*y*z*(-3*(x**2)-3*(y**2)-(z**2))+vz*x*(-(x**2)*(z**2)-(y**2)*(z**2)+(x**4)+2*(x**2)*(y**2)+(y**4)))/((R**4)*(r**3))
    J6[6][2] = (x*z)/((R**2)*r)
    J6[6][3] = -(vx*x*y*z*(-3*(x**2)-3*(y**2)-(z**2))+vy*z*((x**4)-(x**2)*(y**2)-2*(y**4)+(x**2)*(z**2))+vz*y*(-(x**2)*(z**2)-(y**2)*(z**2)+(x**4)+2*(x**2)*(y**2)+(y**4)))/((R**4)*(r**3))
    J6[6][4] = (y*z)/((R**2)*r)
    J6[6][5] = -(vx*((x**3)+x*(y**2)-x*(z**2))+vy*y*((x**2)+(y**2)-(z**2))+vz*z*2*((x**2)+(y**2)))/((R**4)*r)
    J6[6][6] = r/(R**2)

    return J6
