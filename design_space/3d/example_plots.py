from math import sin, cos, pi, sqrt
import itertools

import matplotlib.pyplot as plt
import numpy as np

import giuseppe as gp
from matplotlib import rcParams, gridspec, rc

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

cmap = plt.colormaps['viridis']

fig = plt.figure(figsize=(6.5, 3.5))
fig.suptitle('3D Sampling of Continuation Direction')

gs = gridspec.GridSpec(1, 2)

ax_rand = fig.add_subplot(gs[0, 0], projection='3d')
ax_spiral = fig.add_subplot(gs[0, 1], projection='3d')

np.random.seed(1927)

n_variables = 3
n_samples = 250

d = np.random.normal(size=(3, n_samples))
d /= np.linalg.norm(d, axis=0)

ax_rand.plot(d[0, :], d[1, :], d[2, :], '*')
ax_rand.set_title('Random Sampling')

# n_turns = 8
# n_samples_circle = 15
#
# points = []
# for n in range(n_turns * n_samples_circle):
#     points.append(
#             (
#                 -cos(n / (n_turns * n_samples_circle) * pi),
#                 sin(n / n_samples_circle * 2 * pi) * sin(n / (n_turns * n_samples_circle) * pi),
#                 cos(n / n_samples_circle * 2 * pi) * sin(n / (n_turns * n_samples_circle) * pi),
#             )
#     )
#
# spiral_points = np.array(points).T

n_spiral = 250

golden_angle = pi * (3 - sqrt(5))

points = []
for n in range(n_spiral):
    theta = (golden_angle * n)
    z = (1 - 1 / n_spiral) * (1 - 2 * n / (n_spiral - 1))
    radius = sqrt(1 - z**2)
    points.append((radius * cos(theta), radius * sin(theta), z))


spiral_points = np.array(points).T

ax_spiral.plot(spiral_points[0, :], spiral_points[1, :], spiral_points[2, :], '*', color='C1')
ax_spiral.set_title('Spiral Sampling')

plt.show()
