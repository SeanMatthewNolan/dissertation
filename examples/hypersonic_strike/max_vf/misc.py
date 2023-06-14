# from math import sin, cos, asin, sqrt, atan2, pi
import pickle

from numpy import pi, sin, cos, sqrt
from numpy import arcsin as asin
from numpy import arctan2 as atan2

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import HostAxes
from mpl_toolkits.axisartist import GridHelperCurveLinear


def wrap_ang(ang):
    """
    Wraps angles to fall betweem -pi and pi
    :param ang:
    :return:
    """
    return (ang + pi) % (2 * pi) - pi


def hav(theta):
    """
    Haversine function
    :param theta:
    :return: haversine(theta)
    """
    return (1 - cos(theta))/2


def ahav(x):
    """
    Inverse haversine function
    :param x:
    :return:
    """
    return 2*asin(sqrt(x))


def circle_ang_dist(lat1, long1, lat2, long2):
    """
    Calculate angular distance between two lat-lon coordinates
    :param lat1: latitude for first point
    :param long1: longitude for first point
    :param lat2: latitude for second point
    :param long2: longitude for second point
    :return: angular distance
    """
    hav_theta = hav(lat2 - lat1) + cos(lat1)*cos(lat2)*hav(long2 - long1)
    return ahav(hav_theta)


def calc_bearing(lat1, long1, lat2, long2):
    """
    Calculate bearing angle from first coorginates to second coordinates
    :param lat1: latitude for first point
    :param long1: longitude for first point
    :param lat2: latitude for second point
    :param long2: longitude for second point
    :return: bearing angle
    """
    dlong = long2 - long1
    x = cos(lat2) * sin(dlong)
    y = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlong)
    return atan2(x, y)


def curvelinear(fig):
    """Polar projection, but in a rectangular box."""
    # see demo_curvelinear_grid.py for details
    tr = Affine2D().scale(np.pi / 180., 1.) + PolarAxes.PolarTransform()

    extreme_finder = angle_helper.ExtremeFinderCycle(
            20, 20, lon_cycle=360, lat_cycle=None, lon_minmax=None, lat_minmax=(0, np.inf))

    grid_locator1 = angle_helper.LocatorDMS(12)

    tick_formatter1 = angle_helper.FormatterDMS()

    grid_helper = GridHelperCurveLinear(
            tr, extreme_finder=extreme_finder, grid_locator1=grid_locator1, tick_formatter1=tick_formatter1)

    ax1 = fig.add_subplot(axes_class=HostAxes, grid_helper=grid_helper)

    # Now creates floating axis

    # # floating axis whose first coordinate (theta) is fixed at 60
    # ax1.axis["lat"] = axis = ax1.new_floating_axis(0, 0)
    # axis.label.set_text(r"$\theta = 60^{\circ}$")
    # axis.label.set_visible(True)
    #
    # # floating axis whose second coordinate (r) is fixed at 6
    # ax1.axis["lon"] = axis = ax1.new_floating_axis(1, 6)
    # axis.label.set_text(r"$r = 6$")

    ax1.set_aspect(1.)
    # ax1.set_xlim(-5, 12)
    # ax1.set_ylim(-5, 10)

    ax1.grid(True)

    return ax1


def load(filename):
    with open(filename, 'rb') as file:
        sols = pickle.load(file)

    return sols
