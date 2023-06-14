import os.path
import re
import random

import numpy as np
from matplotlib import rcParams, rc, gridspec
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, LogNorm, SymLogNorm
from scipy.integrate import simpson

import giuseppe

from misc import circle_ang_dist, calc_bearing, curvelinear


os.chdir(os.path.dirname(__file__))

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

SKIP = 1
RE = 20_902_900 / 6076.118

rad2deg = 180/3.141592653589793

gs = gridspec.GridSpec(2, 2)

fig = plt.figure(figsize=(6.5, 6.5))

ax_h = curvelinear(fig, gs[0, 0])
ax_alpha = curvelinear(fig, gs[0, 1])
ax_beta = curvelinear(fig, gs[1, 0])
ax_h_uu = curvelinear(fig, gs[1, 1])

cmap_h = plt.colormaps['viridis']
cmap_alpha = plt.colormaps['plasma']
cmap_beta = plt.colormaps['cividis']
cmap_h_uu = plt.colormaps['magma']

directory = 'star_cav_h_45_00'
files = [f'./{directory}/{file}'
         for file in os.listdir(f'./{directory}') if re.match(r'(ray)(\d*)(.bin)', file)]

data = [giuseppe.load_sol_set(file) for file in files[::]]

h_penalties = []
alpha_penalties = []
beta_penalties = []
h_uus = []

all_sols = []
for sols in data:
    all_sols += sols

    h_penalties.append([])
    alpha_penalties.append([])
    beta_penalties.append([])
    h_uus.append([])
    for sol in sols:
        h_penalties[-1].append(simpson(sol.aux['Penalty: h'], sol.t))
        alpha_penalties[-1].append(simpson(sol.aux['Penalty: alpha'], sol.t))
        beta_penalties[-1].append(simpson(sol.aux['Penalty: beta'], sol.t))
        h_uus[-1].append(min([min(np.linalg.eig(h_uu)[0]) for h_uu in sol.aux['H_uu']]))

norm_h = LogNorm(np.min(np.concatenate(h_penalties)), np.max(np.concatenate(h_penalties)))
norm_alpha = LogNorm(np.min(np.concatenate(alpha_penalties)), np.max(np.concatenate(alpha_penalties)))
norm_beta = SymLogNorm(1e-5, vmin=np.min(np.concatenate(beta_penalties)), vmax=np.max(np.concatenate(beta_penalties)))
norm_h_uu = LogNorm(np.min(np.concatenate(h_uus)), np.max(np.concatenate(h_uus)))

for sols, h_p, alpha_p, beta_p, h_uu_min in zip(data, h_penalties, alpha_penalties, beta_penalties, h_uus):
    star_theta = []
    star_phi = []

    for sol in sols[::SKIP]:
        theta = sol.x[2, :]
        phi = sol.x[1, :]
        # norms.append(np.min(sol.aux['H_uu'][:, 1, 1]))

        star_theta.append(theta[-1])
        star_phi.append(phi[-1])

    star_theta = np.array(star_theta)
    star_phi = np.array(star_phi)

    r = RE * circle_ang_dist(0., 0., star_phi, star_theta)
    bear = calc_bearing(0., 0., star_phi, star_theta)

    x = r * np.cos(bear)
    y = r * np.sin(bear)

    points = np.array([x, y]).T.reshape(-1, 1, 2)

    segments = np.concatenate([points[:-1, :], points[1:, :]], axis=1)

    lc_h = LineCollection(segments, cmap=cmap_h, norm=norm_h)
    lc_h.set_array(h_p)
    ax_h.add_collection(lc_h)

    lc_alpha = LineCollection(segments, cmap=cmap_alpha, norm=norm_alpha)
    lc_alpha.set_array(alpha_p)
    ax_alpha.add_collection(lc_alpha)

    lc_beta = LineCollection(segments, cmap=cmap_beta, norm=norm_beta)
    lc_beta.set_array(beta_p)
    ax_beta.add_collection(lc_beta)

    lc_h_uu = LineCollection(segments, cmap=cmap_h_uu, norm=norm_h_uu)
    lc_h_uu.set_array(h_uu_min)
    ax_h_uu.add_collection(lc_h_uu)

fig.colorbar(plt.cm.ScalarMappable(norm=norm_h, cmap=cmap_h), label=r'Penalty $\int P_{h}\; dt$',
             ax=ax_h, orientation='horizontal')
fig.colorbar(plt.cm.ScalarMappable(norm=norm_alpha, cmap=cmap_alpha), label=r'Penalty $\int P_{\alpha}\; dt$',
             ax=ax_alpha, orientation='horizontal')
fig.colorbar(plt.cm.ScalarMappable(norm=norm_beta, cmap=cmap_beta), label=r'Penalty $\int{P_{\beta}}\; dt$',
             ax=ax_beta, orientation='horizontal')
fig.colorbar(plt.cm.ScalarMappable(norm=norm_h_uu, cmap=cmap_h_uu), label=r'Min. Eigenvalue of $H_{\bm{uu}}$',
             ax=ax_h_uu, orientation='horizontal')

ax_h.set_title('Altitude Constraint\nPenalty')
ax_h.set_xlabel(r'Azimuth [deg]')
ax_h.set_ylabel(r'Range [NM]')
ax_h.autoscale()

ax_alpha.set_title('AOA Constraint\nPenalty')
ax_alpha.set_xlabel(r'Azimuth [deg]')
ax_alpha.set_ylabel(r'Range [NM]')
ax_alpha.autoscale()

ax_beta.set_title('Bank Constraint\nPenalty')
ax_beta.set_xlabel(r'Azimuth [deg]')
ax_beta.set_ylabel(r'Range [NM]')
ax_beta.autoscale()

ax_h_uu.set_title(r'$H_{\bm{uu}}$'+'\nMin. Eigenvalue')
ax_h_uu.set_xlabel(r'Azimuth [deg]')
ax_h_uu.set_ylabel(r'Range [NM]')
ax_h_uu.autoscale()

fig.subplots_adjust(
        top=0.92,
        bottom=0.046,
        left=0.08,
        right=0.977,
        hspace=0.268,
        wspace=0.194
)

plt.show()
