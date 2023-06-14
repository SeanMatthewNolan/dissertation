from itertools import product

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from giuseppe.guess_generation.gauss_newton import gauss_newton


# https://homepage.univie.ac.at/herwig.hauser/gallery.html

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 8
rcParams['mathtext.fontset'] = 'stix'


def sphere(arr):
    return np.expand_dims(arr[0]**2 + arr[1]**2 + arr[2] ** 2 - 1**2, axis=0)


def torus(arr):
    x, y, z = arr
    return np.array([
        (x**2 + y**2 + z**2 + 1 - 0.5**2) ** 2 - 4 * (x**2 + y**2)
    ])


def heart(arr):
    x, y, z = arr
    return np.array([
        (x**2 + 9/4 * y**2 + z**2 - 1) ** 3 - x ** 2 * z ** 3 - 9 / 89 * y ** 2 * z ** 3
    ])


def paraboloid(arr):
    x, y, z = arr
    return np.array([
        z ** 2 - x ** 2 - y ** 2 - 1
    ])


def subway(arr):
    x, y, z = arr
    return np.array([
        x**2 * y**2 + (z**2 - 1)**3
    ])


def planes(arr):
    return np.array([
        arr[0] * arr[1] * arr[2]
    ])


def tangle(arr):
    x, y, z = arr
    a, b, c = 0.0, -5.0, 11.8
    return np.expand_dims(np.array(x**4+y**4+z**4+a*(x**2+y**2+z**2)**2+b*(x**2+y**2+z**2)+c), axis=0)


def crixxi(arr):
    x, y, z = arr
    return np.array([
        (y**2 + z**2 - 1)**2 +(x**2 + y**2 - 1)**3
    ])


def set_lims(arr, pad_ratio=0.05):
    a, b = arr[0], arr[-1]
    d = b - a
    pad = d * pad_ratio
    return a - pad, b + pad


_side = np.linspace(-2, 2, 11)
x_list, y_list, z_list = [], [], []
for _x, _y, _z in product(_side, _side, _side):
    x_list.append(_x), y_list.append(_y), z_list.append(_z)

states = np.array((x_list, y_list, z_list))
# states = np.array((x_list, y_list, z_list)) + np.random.random((11 ** 3,)) / 100


def make_plot(points, filename):
    fig1 = plt.figure(figsize=(3, 3))
    ax111 = fig1.add_subplot(111, projection='3d')
    ax111.set_proj_type('ortho')
    ax111.set_aspect('equal')
    ax111.set(
            xlim=set_lims(states[0, :]),
            ylim=set_lims(states[1, :]),
            zlim=set_lims(states[2, :]),
    )

    ax111.plot(points[0, :], points[1, :], points[2, :], '.')

    ax111.set_xlabel(r'$x$')
    ax111.set_ylabel(r'$y$')
    ax111.set_zlabel(r'$z$')

    fig1.subplots_adjust(
            top=0.95,
            bottom=0.1,
            left=0.1,
            right=0.9,
    )

    plt.show()

    fig1.savefig(f'images/{filename}.eps')


make_plot(states, 'grid')

for constraint, name in [
    (sphere, 'sphere'), (torus, 'torus'), (heart, 'heart'), (planes, 'planes'), (tangle, 'tangle'),
    (paraboloid, 'paraboloid'), (subway, 'subway'), (crixxi, 'crixxi')
]:
    proj_states = np.empty_like(states)

    for idx, states_i in enumerate(states.T):
        try:
            proj_states[:, idx] = gauss_newton(
                    constraint, states_i, use_line_search=True, max_steps=51, abs_tol=1e-8, rel_tol=1e-8, verbose=True)
        except RuntimeError:
            #  proj_states[:, idx] = states_i
            proj_states[:, idx] = np.nan
            print(states_i)

    # proj_states = np.array(proj_states).T
    make_plot(proj_states, name)
