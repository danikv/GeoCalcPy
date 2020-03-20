import math


def cartesian2directional(vx, vy, vz):
    speed = math.sqrt(vx * vx + vy * vy)
    course = math.degrees(math.atan2(vy, vx))
    return speed, course, vz


def directional2cartesian(speed, course, vz):
    course_in_rad = math.radians(course)
    vx = speed * math.cos(course_in_rad)
    vy = speed * math.sin(course_in_rad)
    return vx, vy, vz
