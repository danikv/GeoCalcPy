import math


def cartesian_to_directional(vx, vy, vz):
    speed = math.sqrt(vx * vx + vy * vy)
    course = math.degrees(math.atan2(vy, vx))
    return speed, course, vz


def directional_to_cartesian(speed, course, vz):
    course_in_rad = math.radians(course)
    vx = speed * math.cos(course_in_rad)
    vy = speed * math.sin(course_in_rad)
    return vx, vy, vz


def enu2spherical(x, y, z, vx, vy, vz):
    R = math.sqrt(math.pow(x,2)+math.pow(y,2)+math.pow(z,2))
    r = math.sqrt(math.pow(x,2)+math.pow(y,2))
    range_rate = (x*vx+y*vy+z*vz)/R
    bearing_rate = (vx*y-x*vy)/math.pow(r,2)
    elevation_rate = (vz*math.pow(r,2)-z*(x*vx+y*vy))/(r*math.pow(R,2))
    return range_rate, bearing_rate, elevation_rate


def spherical2enu(range, bearing, elevation, range_rate, bearing_rate, elevation_rate):
    r = range
    rdot = range_rate
    bdot = bearing_rate
    edot = elevation_rate
    cos_b = math.cos(bearing)
    sin_b = math.sin(bearing)
    cos_e = math.cos(elevation)
    sin_e = math.sin(elevation)

    vx = rdot*sin_b*cos_e + r*bdot*cos_b*cos_e - r*edot*sin_b*sin_e
    vy = rdot*cos_b*cos_e - r*bdot*sin_b*cos_e - r*edot*cos_b*sin_e
    vz = rdot*sin_e + r*edot*cos_e
    return vx, vy, vz

#range_rate, bearing_rate, elevation_rate = enu2spherical(1, 1, 1, 1, 1, 1)
#print(range_rate)
#print(bearing_rate)
#print(elevation_rate)

vx, vy, vz = spherical2enu(1, 0, math.pi/2, 1, 0, 0)
print(vx)
print(1-vy)
print(vz)