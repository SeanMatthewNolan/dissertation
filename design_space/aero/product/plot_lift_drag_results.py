import os.path
import re
import pickle
import random

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib import rc, rcParams, gridspec


os.chdir(os.path.dirname(__file__))

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

SKIP = 100
RE = 20_902_900

rad2deg = 180/3.141592653589793

fig1 = plt.figure(figsize=(6.5, 7.5))
gs = gridspec.GridSpec(3, 2, height_ratios=[2, 1, 0.1])

ax_heat = fig1.add_subplot(gs[0, :])
ax_1 = fig1.add_subplot(gs[1, 0])
ax_2 = fig1.add_subplot(gs[1, 1])
cax = fig1.add_subplot(gs[2, :])

terminal_v = []
lift_per = []
drag_per = []

t_to_plot = []
phi_to_plot = []
theta_to_plot = []
h_to_plot = []
alpha_to_plot = []
bank_to_plot = []
v_to_plot = []

terminal_vs = []
lift_pers = []
drag_pers = []
with open('set_phi_060_theta_000.data', 'rb') as f:
    data = pickle.load(f)

for sols in data['inter_node_sets']:
    terminal_v = []
    lift_per = []
    drag_per = []

    for sol in sols:
        t = sol.t

        h = sol.x[0, :]
        phi = sol.x[1, :] * rad2deg
        theta = sol.x[2, :] * rad2deg

        v = sol.x[3, :]

        alpha = sol.u[0, :] * rad2deg

        terminal_v.append(v[-1])
        lift_per.append((1 + sol.k[22]) * 100)
        drag_per.append((1 + sol.k[23]) * 100)

    terminal_vs.append(terminal_v)
    lift_pers.append(lift_per)
    drag_pers.append(drag_per)

terminal_vs_nodes = []
lift_pers_nodes = []
drag_pers_nodes = []
for sol in data['node_sols']:
    t = sol.t

    h = sol.x[0, :]
    phi = sol.x[1, :] * rad2deg
    theta = sol.x[2, :] * rad2deg

    v = sol.x[3, :]

    alpha = sol.u[0, :] * rad2deg
    bank = sol.u[1, :] * rad2deg

    t_to_plot.append(t)
    phi_to_plot.append(phi)
    theta_to_plot.append(theta)
    h_to_plot.append(h)
    alpha_to_plot.append(alpha)
    bank_to_plot.append(bank)
    v_to_plot.append(v)

    terminal_vs_nodes.append(v[-1])
    lift_pers_nodes.append((1 + sol.k[22]) * 100)
    drag_pers_nodes.append((1 + sol.k[23]) * 100)


all_terminal_vs = np.concatenate(terminal_vs)
max_v = max(all_terminal_vs)
min_v = min(all_terminal_vs)
cmap = plt.colormaps['viridis']
norm = Normalize(min_v, max_v)

def generate_line_collection(_x, _y, _norms):
    points = np.array([_x, _y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1, :], points[1:, :]], axis=1)
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(_norms)
    return lc


for _lift_per, _drag_per, _terminal_v in zip(lift_pers, drag_pers, terminal_vs):
    ax_heat.add_collection(generate_line_collection(_drag_per, _lift_per,  _terminal_v))


ax_heat.scatter(drag_pers_nodes, lift_pers_nodes, marker='x', c=norm(terminal_vs_nodes))

ax_heat.autoscale()

for t, phi, h, alpha, v in zip(t_to_plot, phi_to_plot, h_to_plot, alpha_to_plot, v_to_plot):
    ax_1.plot(phi, h / 1e5, color=cmap(v[-1] / max_v), linewidth=1)
    ax_2.plot(t, alpha, color=cmap(v[-1] / max_v), linewidth=1)
fig1.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, max_v), cmap=cmap),
              cax=cax, label='Terminal Velocity [ft/s]', orientation='horizontal')

ax_heat.set_xlabel(r'Percent of Nominal Drag, $k_{\text{D, Mult.}}$ [\%]')
ax_heat.set_ylabel(r'Percent of Nominal Lift,, $k_{\text{L, Mult.}}$ [\%]')

# ax221.set_title('Maximum Terminal Energy Sample Space')
ax_1.set_xlabel(r'Downrange $\phi$ [deg]')
ax_1.set_ylabel(r'Altitude $h$ [$\times 10^5$ ft]')

ax_2.set_xlabel(r'Time $t$ [s]')
ax_2.set_ylabel(r'Angle of Attack $\alpha$ [deg]')

fig1.suptitle('Sparse Product Space Continuation:\nVariation of Terminal Velocity with Aerodynamic Functions')

fig1.subplots_adjust(
        top=0.918,
        bottom=0.075,
        left=0.101,
        right=0.977,
        hspace=0.355,
        wspace=0.243
)

plt.show()
