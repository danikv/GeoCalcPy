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


def cartesian_to_polar(vx, vy, vz):
    range_rate = 1
    bearing_rate = 1
    elevation_rate = 1
    return range_rate, bearing_rate, elevation_rate


def polar_to_ned_velocity(range, bearing, elevation, range_rate, bearing_rate, elevation_rate):
    vx = 1
    vy = 1
    vz = 1
    return vx, vy, vz
