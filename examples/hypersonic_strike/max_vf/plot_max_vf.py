import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec


def haversine(_ang):
    return (1 - np.cos(_ang)) / 2


def ahav(x):
    return 2*np.arcsin(np.sqrt(x))


def compute_range(downrange, crossrange):
    return ahav(haversine(crossrange) + np.cos(crossrange) * haversine(downrange))


rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

cav_h_set = gp.load_sol_set('cav_h.bin')
cav_l_set = gp.load_sol_set('cav_l.bin')

cav_h_sol = cav_h_set[-1]
cav_l_sol = cav_l_set[-1]

fig = plt.figure(figsize=(6.5, 7.5))
title = fig.suptitle('Nominal CAV Maximum Terminal Energy Trajectories')

gs = gridspec.GridSpec(3, 2)

ax_path = fig.add_subplot(gs[0, :])
ax_ground = fig.add_subplot(gs[1, 1])

ax_hv = fig.add_subplot(gs[1, 0])

ax_alpha = fig.add_subplot(gs[2, 0])
ax_beta = fig.add_subplot(gs[2, 1])

ax_path.plot(compute_range(cav_h_sol.x[1, :], cav_h_sol.x[2, :]) * 20_902_900 / 6076.11, cav_h_sol.x[0, :] / 1e5,
             label=r'CAV-H')
ax_path.plot(compute_range(cav_l_sol.x[1, :], cav_l_sol.x[2, :]) * 20_902_900 / 6076.11, cav_l_sol.x[0, :] / 1e5,
             label=r'CAV-L')
ax_path.set_xlabel(r'Range to Origin [NM]')
ax_path.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

fig.legend(loc='lower center', ncols=3)

ax_ground.plot(cav_h_sol.x[1, :] * 180 / np.pi, cav_h_sol.x[2, :] * 180 / np.pi,
               label=r'CAV-H')
ax_ground.plot(cav_l_sol.x[1, :] * 180 / np.pi, cav_l_sol.x[2, :] * 180 / np.pi,
               label=r'CAV-L')
ax_ground.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_ground.set_ylabel(r'Crossrange, $\theta$ [deg]')

ax_hv.plot(cav_h_sol.x[3, :] / 1e3, cav_h_sol.x[0, :] / 1e5, label=r'CAV-H')
ax_hv.plot(cav_l_sol.x[3, :] / 1e3, cav_l_sol.x[0, :] / 1e5, label=r'CAV-L')
ax_hv.set_xlabel(r'Velocity, $v$ [$\times10^3$ ft/s]')
ax_hv.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

ax_alpha.plot(cav_h_sol.t, cav_h_sol.u[0, :] * 180 / np.pi, label=r'CAV-H')
ax_alpha.plot(cav_l_sol.t, cav_l_sol.u[0, :] * 180 / np.pi, label=r'CAV-L')
ax_alpha.set_xlabel(r'Time $t$ [s]')
ax_alpha.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')

ax_beta.plot(cav_h_sol.t, cav_h_sol.u[1, :] * 180 / np.pi, label=r'CAV-H')
ax_beta.plot(cav_l_sol.t, cav_l_sol.u[1, :] * 180 / np.pi, label=r'CAV-L')
ax_beta.set_xlabel(r'Time $t$ [s]')
ax_beta.set_ylabel(r'Bank Angle, $\beta$ [deg]')

fig.subplots_adjust(
        top=0.925,
        bottom=0.125,
        left=0.125,
        right=0.975,
        hspace=0.4,
        wspace=0.4
)

plt.show()
