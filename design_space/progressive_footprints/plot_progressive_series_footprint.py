import os

import giuseppe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, rcParams
from matplotlib.colors import Normalize

from misc import circle_ang_dist, calc_bearing, curvelinear, load

os.chdir(os.path.dirname(__file__))

SKIP = 1
RE = 20_902_900 / 6076.118

rad2deg = 180/3.141592653589793

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

fig = plt.figure(figsize=(4, 4.5))
ax = curvelinear(fig)

data = {
    0: load('footprint_t0.bin'),
    50: load('footprint_t50.bin'),
    100: load('footprint_t100.bin'),
    150: load('footprint_t150.bin'),
    200: load('footprint_t200.bin'),
    250: load('footprint_t250.bin'),
    300: load('footprint_t300.bin'),
    350: load('footprint_t350.bin'),
    400: load('footprint_t400.bin'),
    450: load('footprint_t450.bin'),
    500: load('footprint_t500.bin'),
    550: load('footprint_t550.bin'),
    600: load('footprint_t600.bin'),
    650: load('footprint_t650.bin'),
    700: load('footprint_t700.bin'),
    750: load('footprint_t750.bin'),
    800: load('footprint_t800.bin'),
    850: load('footprint_t850.bin'),
    900: load('footprint_t900.bin'),
    950: load('footprint_t950.bin'),
    1000: load('footprint_t1000.bin'),
    1050: load('footprint_t1050.bin'),
    1100: load('footprint_t1100.bin'),
    1150: load('footprint_t1150.bin'),
    1200: load('footprint_t1200.bin'),
    1250: load('footprint_t1250.bin'),
    1300: load('footprint_t1300.bin'),
    1350: load('footprint_t1350.bin'),
    1400: load('footprint_t1400.bin'),
    1450: load('footprint_t1450.bin'),
    1500: load('footprint_t1500.bin'),
    1550: load('footprint_t1550.bin'),
}

ref_sol = giuseppe.load_sol('ref_traj.bin')

ref_phi = ref_sol.x[1, :]
ref_theta = ref_sol.x[2, :]

ref_range = RE * circle_ang_dist(0., 0., ref_phi, ref_theta)
ref_bear = calc_bearing(0., 0., ref_phi, ref_theta)

ax.plot(ref_range * np.cos(ref_bear), ref_range * np.sin(ref_bear), color='C1', linestyle='-', label='Reference')

ax.legend()

cmap = plt.colormaps['viridis']

for ti, sols in data.items():
    label = f'{ti}'
    color = cmap(ti/1600)

    footprint_theta = []
    footprint_phi = []

    for sol in sols[::SKIP]:
        theta = sol.x[1, :]
        phi = sol.x[2, :]

        footprint_theta.append(theta[-1])
        footprint_phi.append(phi[-1])

    footprint_theta = np.array(footprint_theta)
    footprint_phi = np.array(footprint_phi)

    r = RE * circle_ang_dist(0., 0., footprint_theta, footprint_phi)
    bear = calc_bearing(0., 0., footprint_theta, footprint_phi)

    ax.plot(r * np.cos(bear), r * np.sin(bear), label=label, color=color)

# ax.plot(0, 0, label='Initial Point', marker='+', linestyle='')

ax.set_title('Evolution of Footprints Over Time')
ax.set_xlabel(r'Angle to Glide Origin')
ax.set_ylabel(r'Range to Glide Origin [NM]')

ax.autoscale()

# ax.legend(title='Time [s]', loc='upper left', bbox_to_anchor=(1, 1), fontsize=10, title_fontsize=10)
fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, ref_sol.t[-1]), cmap=cmap), ax=ax, orientation='horizontal',
             aspect=50, label=r'Initial Time from Reference [sec]')

fig.tight_layout()

plt.show()
