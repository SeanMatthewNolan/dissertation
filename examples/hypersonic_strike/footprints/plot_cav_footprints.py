import os
import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec

import giuseppe.data_classes.solution_set
from misc import circle_ang_dist, calc_bearing, curvelinear, load

os.chdir(os.path.dirname(__file__))

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

SKIP = 1
# RE = 20_902_900 / 5280
RE = 20_902_900 / 6076.118
rad2deg = 180 / 3.141592653589793

cav_h_sols = giuseppe.data_classes.load_sol_set('cav_h_footprint.bin')
cav_l_sols = giuseppe.data_classes.load_sol_set('cav_l_footprint.bin')

fig = plt.figure(figsize=(6.5, 5))
ax = curvelinear(fig)

vfss = []
for sols, label in ((cav_h_sols, 'CAV-H'), (cav_l_sols, 'CAV-L')):
    vfs = []
    footprint_theta = []
    footprint_phi = []

    for sol in sols[::SKIP]:
        theta = sol.x[1, :]
        phi = sol.x[2, :]

        footprint_theta.append(theta[-1])
        footprint_phi.append(phi[-1])
        vfs.append(sol.x[3, -1])

    footprint_theta = np.array(footprint_theta)
    footprint_phi = np.array(footprint_phi)

    r = RE * circle_ang_dist(0., 0., footprint_theta, footprint_phi)
    bear = calc_bearing(0., 0., footprint_theta, footprint_phi)

    ax.plot(r * np.cos(bear), r * np.sin(bear), label=label)
    vfss.append(vfs)


ax.set_title('Reachability Footprints for CAVs')
ax.set_xlabel(r'Angle to Glide Origin [deg]')
ax.set_ylabel(r'Range to Glide Origin [NM]')
ax.legend()

ax.autoscale()

fig.subplots_adjust(
        top=0.925,
        bottom=0.125,
        left=0.125,
        right=0.975,
        hspace=0.4,
        wspace=0.4
)

plt.show()
