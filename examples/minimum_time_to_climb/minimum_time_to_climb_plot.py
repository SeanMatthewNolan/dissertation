import pickle

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import numpy as np

from giuseppe.utils.examples import Atmosphere1976

from minimum_time_to_climb import S, adiff_dual
from lookup_tables import thrust_table_bspline, eta_table_bspline_expanded, CLalpha_table_bspline_expanded,\
    CD0_table_bspline_expanded, temp_table_bspline, dens_table_bspline

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

SMALL_FIGSIZE = (6.5, 3)
MED_FIGSIZE = (6.5, 5)
LARGE_FIGSIZE = (6.5, 7.5)
T_LAB = 'Time [sec]'

DATA = 0

if DATA == 0:
    with open('sol_set.data', 'rb') as file:
        sols = pickle.load(file)
        sol = sols[-1]
elif DATA == 1:
    with open('seed.data', 'rb') as file:
        sol = pickle.load(file)
elif DATA == 2:
    with open('guess.data', 'rb') as file:
        sol = pickle.load(file)

r2d = 180 / np.pi
d2r = np.pi / 180

h = sol.x[0, :]
V = sol.x[1, :]
alpha = sol.u[0, :]
alpha_hat = alpha * r2d

atm = Atmosphere1976(use_metric=False)

# T = np.asarray([atm.temperature(alt) for alt in h])
# rho = np.asarray([atm.density(alt) for alt in h])
T = np.asarray(temp_table_bspline(h)).flatten()
rho = np.asarray(dens_table_bspline(h)).flatten()

a = np.sqrt(atm.specific_heat_ratio * atm.gas_constant * T)

M = V/a
Qdyn = 0.5 * rho * V**2

thrust = np.asarray(thrust_table_bspline(np.vstack((M.T, h.T)))).flatten()
eta = np.asarray(eta_table_bspline_expanded(M)).flatten()
CLalpha = np.asarray(CLalpha_table_bspline_expanded(M)).flatten()
CD0 = np.asarray(CD0_table_bspline_expanded(M)).flatten()

CD = CD0 + eta * CLalpha * alpha_hat**2
CL = CLalpha * alpha_hat

LoD = CL / CD

drag = 0.5 * CD * S * rho * V**2
lift = 0.5 * CL * S * rho * V**2

# FIGURE 1 (STATES)
fig1 = plt.figure(figsize=LARGE_FIGSIZE)
title = fig1.suptitle('Minimum Time to Climb Problem')

# Alt. vs. Time
ax1 = fig1.add_subplot(321)
ax1.plot(sol.t, sol.x[0, :] / 1_000)
xlabel_1 = ax1.set_xlabel(T_LAB)
ylabel_1 = ax1.set_ylabel('Altitude, $h$ [1000 ft]')
ax1.grid()

# Velocity vs. Time
ax2 = fig1.add_subplot(322)
ax2.plot(sol.t, sol.x[1, :] / 100)
xlabel_2 = ax2.set_xlabel(T_LAB)
ylabel_2 = ax2.set_ylabel('Velocity, $v$ [100 ft/s]')
ax2.grid()

# FPA vs. Time
ax3 = fig1.add_subplot(323)
ax3.plot(sol.t, sol.x[2, :] * r2d)
xlabel_3 = ax3.set_xlabel(T_LAB)
ylabel_3 = ax3.set_ylabel(r'FPA, $\gamma$ [deg]')
ax3.grid()

# Weight vs. Time
ax4 = fig1.add_subplot(324)
ax4.plot(sol.t, sol.x[3, :] / 10_000)
xlabel_4 = ax4.set_xlabel(T_LAB)
ylabel_4 = ax4.set_ylabel(r'Weight, $w$  [10,000 lb]')
ax4.grid()

# AoA vs. Time
ax5 = fig1.add_subplot(325)
ax5.plot(sol.t, sol.u[0, :] * r2d)
xlabel_5 = ax5.set_xlabel(T_LAB)
ylabel_5 = ax5.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')
ax5.grid()

# Alt. Vs. Velocity
ax6 = fig1.add_subplot(326)
ax6.plot(sol.x[1, :] / 100, sol.x[0, :] / 1_000)
xlabel_6 = ax6.set_xlabel('Velocity, $v$ [100 ft/s]')
ylabel_6 = ax6.set_ylabel('Altitude, $h$ [1000 ft]')
ax6.grid()

fig1.tight_layout()

plt.show()
