import os
import pickle

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import gridspec, rcParams

import giuseppe.utils.compilation
from giuseppe.numeric_solvers import SciPySolver
from giuseppe.guess_generation import auto_propagate_guess
from giuseppe.problems.input import StrInputProb
from giuseppe.problems.symbolic import SymDual, PenaltyConstraintHandler
from giuseppe.utils import Timer

giuseppe.utils.compilation.JIT_COMPILE = False

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 8
rcParams['mathtext.fontset'] = 'stix'

os.chdir(os.path.dirname(__file__))  # Set directory to current location

ocp = StrInputProb()

ocp.set_independent('t')

ocp.add_expression('r', 're + h')
ocp.add_expression('g', 'mu / r**2')
ocp.add_expression('rho', 'rho_0 * exp(-h / h_ref)')
ocp.add_expression('dyn_pres', '1 / 2 * rho * v ** 2 ')
ocp.add_expression('lift', 'c_l * s_ref * dyn_pres')
ocp.add_expression('drag', 'c_d * s_ref * dyn_pres')
ocp.add_expression('c_l', 'a_0 + a_1 * alpha_hat')
ocp.add_expression('c_d', 'b_0 + b_1 * alpha_hat + b_2 * alpha_hat**2')
ocp.add_expression('alpha_hat', 'alpha * 180 / pi')

ocp.add_state('h', 'v * sin(gamma)')
ocp.add_state('phi', 'v * cos(gamma) * sin(psi) / (r * cos(theta))')
ocp.add_state('theta', 'v * cos(gamma) * cos(psi) / r')
ocp.add_state('v', '-drag / m - g * sin(gamma)')
ocp.add_state('gamma', 'lift * cos(beta) / (m * v) + cos(gamma) * (v / r - g / v)')
ocp.add_state('psi', 'lift * sin(beta)/(m * v * cos(gamma)) + v * cos(gamma) * sin(psi) * sin(theta)/(r * cos(theta))')

ocp.add_control('alpha')
ocp.add_control('beta')

ocp.add_constant('rho_0', 0.002378)
ocp.add_constant('h_ref', 23_800)
ocp.add_constant('re', 20_902_900)
ocp.add_constant('m', 203_000 / 32.174)
ocp.add_constant('mu', 0.14076539e17)

ocp.add_constant('a_0', -0.20704)
ocp.add_constant('a_1', 0.029244)
ocp.add_constant('b_0', 0.07854)
ocp.add_constant('b_1', -0.61592e-2)
ocp.add_constant('b_2', 0.621408e-3)
ocp.add_constant('s_ref', 2690)

ocp.add_constant('xi', 0)

ocp.add_constant('eps_alpha', 1e-5)
ocp.add_constant('alpha_min', -80 / 180 * 3.1419)
ocp.add_constant('alpha_max', 80 / 180 * 3.1419)

ocp.add_constant('eps_beta', 1e-5)
ocp.add_constant('beta_min', -85 / 180 * 3.1419)
ocp.add_constant('beta_max', 85 / 180 * 3.1419)

ocp.add_constant('h_0', 260_000)
ocp.add_constant('phi_0', 0)
ocp.add_constant('theta_0', 0)
ocp.add_constant('v_0', 25_600)
ocp.add_constant('gamma_0', -1 / 180 * np.pi)
ocp.add_constant('psi_0', np.pi / 2)

ocp.add_constant('h_f', 80_000)
ocp.add_constant('v_f', 2_500)
ocp.add_constant('gamma_f', -5 / 180 * np.pi)

ocp.set_cost('0', '0', '-phi * cos(xi) - theta  * sin(xi)')

# ocp.add_constraint('initial', 't')
# ocp.add_constraint('initial', 'h - h_0')
# ocp.add_constraint('initial', 'phi - phi_0')
# ocp.add_constraint('initial', 'theta - theta_0')
# ocp.add_constraint('initial', 'v - v_0')
# ocp.add_constraint('initial', 'gamma - gamma_0')
# ocp.add_constraint('initial', 'psi - psi_0')
#
# ocp.add_constraint('terminal', 'h - h_f')
# ocp.add_constraint('terminal', 'v - v_f')
# ocp.add_constraint('terminal', 'gamma - gamma_f')

ocp.add_constraint('initial', 't')
ocp.add_constraint('initial', 'h - h_0')
ocp.add_constraint('initial', 'phi - phi_0')
ocp.add_constraint('initial', 'theta - theta_0')
ocp.add_constraint('initial', 'v - v_0')
ocp.add_constraint('initial', 'gamma - gamma_0')
ocp.add_constraint('initial', 'psi - psi_0')

ocp.add_constraint('terminal', 'h - h_f')
ocp.add_constraint('terminal', 'v - v_f')
ocp.add_constraint('terminal', 'gamma - gamma_f')

ocp.add_inequality_constraint('path', 'alpha', lower_limit='alpha_min', upper_limit='alpha_max',
                              regularizer=PenaltyConstraintHandler('eps_alpha', method='sec'))
ocp.add_inequality_constraint('path', 'beta', lower_limit='beta_min', upper_limit='beta_max',
                              regularizer=PenaltyConstraintHandler('eps_beta', method='sec'))

comp_dual = SymDual(ocp, control_method='differential')
solver = SciPySolver(comp_dual)

INCLUDE_INIT = False

t_span = 100
init_guess = giuseppe.guess_generation.initialize_guess(comp_dual, t_span=t_span)
# guess = auto_propagate_guess(comp_dual, control=(15/180*3.14159, 0), max_step=20,
#                              t_span=100, abs_tol=1e-8, rel_tol=1e-8, verbose=True, default_value=-1)
guess = auto_propagate_guess(comp_dual, control=(15/180*3.14159, 0), t_span=t_span, max_step=20,
                             condition_adjoints=True, verbose=True)
sol = solver.solve(guess)

fig = plt.figure(figsize=(6, 8))
fig.suptitle('Comparison Between Guess and Solution', fontsize=10)

gs = gridspec.GridSpec(5, 3, height_ratios=[1, 1, 1, 1, 1], width_ratios=[1, 1, 1])

ax_h   = fig.add_subplot(gs[0, 0])
ax_phi = fig.add_subplot(gs[0, 1])
ax_the = fig.add_subplot(gs[0, 2])
ax_v   = fig.add_subplot(gs[1, 0])
ax_gam = fig.add_subplot(gs[1, 1])
ax_psi = fig.add_subplot(gs[1, 2])

ax_lam_h   = fig.add_subplot(gs[2, 0])
ax_lam_phi = fig.add_subplot(gs[2, 1])
ax_lam_the = fig.add_subplot(gs[2, 2])
ax_lam_v   = fig.add_subplot(gs[3, 0])
ax_lam_gam = fig.add_subplot(gs[3, 1])
ax_lam_psi = fig.add_subplot(gs[3, 2])

ax_alpha = fig.add_subplot(gs[4, 0])
ax_beta  = fig.add_subplot(gs[4, 1])
ax_ham   = fig.add_subplot(gs[4, 2])

linewidth = 2

ax_h.plot(guess.t, guess.x[0, :], label='Guess',    linewidth=linewidth)
ax_h.plot(sol.t,   sol.x[0, :],   label='Solution', linewidth=linewidth)
ax_h.set_xlabel(r'Time $t$, [sec]')
ax_h.set_ylabel(r'Altitude $h$, [ft]')

ax_phi.plot(guess.t, guess.x[1, :] * 180 / 3.14159 , label='Guess', linewidth=linewidth)
ax_phi.plot(sol.t, sol.x[1, :] * 180 / 3.14159 , label='Solution', linewidth=linewidth)
ax_phi.set_xlabel(r'Time $t$, [sec]')
ax_phi.set_ylabel(r'Downrange $\phi$, [deg]')

ax_the.plot(guess.t, guess.x[2, :] * 180 / 3.14159 , label='Guess', linewidth=linewidth)
ax_the.plot(sol.t, sol.x[2, :] * 180 / 3.14159 , label='Solution', linewidth=linewidth)
ax_the.set_xlabel(r'Time $t$, [sec]')
ax_the.set_ylabel(r'Crossrange $\theta$, [deg]')
ax_the.set_ylim((-0.1, 0.1))

ax_v.plot(guess.t, guess.x[3, :], label='Guess',    linewidth=linewidth)
ax_v.plot(sol.t,   sol.x[3, :],   label='Solution', linewidth=linewidth)
ax_v.set_xlabel(r'Time $t$, [sec]')
ax_v.set_ylabel(r'Velocity $v$, [ft/s]')

ax_gam.plot(guess.t, guess.x[4, :] * 180 / 3.14159 , label='Guess', linewidth=linewidth)
ax_gam.plot(sol.t, sol.x[4, :] * 180 / 3.14159 , label='Solution', linewidth=linewidth)
ax_gam.set_xlabel(r'Time $t$, [sec]')
ax_gam.set_ylabel(r'FPA $\gamma$, [deg]')

ax_psi.plot(guess.t, guess.x[5, :] * 180 / 3.14159 , label='Guess', linewidth=linewidth)
ax_psi.plot(sol.t, sol.x[5, :] * 180 / 3.14159 , label='Solution', linewidth=linewidth)
ax_psi.set_xlabel(r'Time $t$, [sec]')
ax_psi.set_ylabel(r'Heading $\psi$, [deg]')
ax_psi.set_ylim((80, 100))

ax_lam_h.plot(guess.t, guess.lam[0, :], label='Guess', linewidth=linewidth)
ax_lam_h.plot(sol.t, sol.lam[0, :],   label='Solution', linewidth=linewidth)
if INCLUDE_INIT:
    ax_lam_h.plot(init_guess.t, init_guess.lam[0, :], ':', label='Initialized Value', linewidth=linewidth)
ax_lam_h.set_xlabel(r'Time $t$, [sec]')
ax_lam_h.set_ylabel(r'$\lambda_h$, [1/ft]')

ax_lam_phi.plot(guess.t, guess.lam[1, :], label='Guess', linewidth=linewidth)
ax_lam_phi.plot(sol.t, sol.lam[1, :], label='Solution', linewidth=linewidth)
if INCLUDE_INIT:
    ax_lam_phi.plot(init_guess.t, init_guess.lam[1, :], ':', label='Initialized Value', linewidth=linewidth)
ax_lam_phi.set_xlabel(r'Time $t$, [sec]')
ax_lam_phi.set_ylabel(r'$\lambda_\phi$, [1]')
if not INCLUDE_INIT:
    ax_lam_phi.set_ylim((-2, 0))

ax_lam_the.plot(guess.t, guess.lam[2, :], label='Guess', linewidth=linewidth)
ax_lam_the.plot(sol.t, sol.lam[2, :], label='Solution', linewidth=linewidth)
if INCLUDE_INIT:
    ax_lam_the.plot(init_guess.t, init_guess.lam[2, :], ':', label='Initialized Value', linewidth=linewidth)
ax_lam_the.set_xlabel(r'Time $t$, [sec]')
ax_lam_the.set_ylabel(r'$\lambda_\theta$, [1]')
if not INCLUDE_INIT:
    ax_lam_the.set_ylim((-0.001, 0.001))

ax_lam_v.plot(guess.t, guess.lam[3, :], label='Guess',    linewidth=linewidth)
ax_lam_v.plot(sol.t,   sol.lam[3, :],   label='Solution', linewidth=linewidth)
if INCLUDE_INIT:
    ax_lam_v.plot(init_guess.t, init_guess.lam[3, :], ':', label='Initialized Value', linewidth=linewidth)
ax_lam_v.set_xlabel(r'Time $t$, [sec]')
ax_lam_v.set_ylabel(r'$\lambda_v$, [s/ft]')

ax_lam_gam.plot(guess.t, guess.lam[4, :], label='Guess', linewidth=linewidth)
ax_lam_gam.plot(sol.t, sol.lam[4, :], label='Solution', linewidth=linewidth)
if INCLUDE_INIT:
    ax_lam_gam.plot(init_guess.t, init_guess.lam[4, :], ':', label='Initialized Value', linewidth=linewidth)
ax_lam_gam.set_xlabel(r'Time $t$, [sec]')
ax_lam_gam.set_ylabel(r'$\lambda_\gamma$, [1]')

ax_lam_psi.plot(guess.t, guess.lam[5, :], label='Guess', linewidth=linewidth)
ax_lam_psi.plot(sol.t, sol.lam[5, :], label='Solution', linewidth=linewidth)
if INCLUDE_INIT:
    ax_lam_psi.plot(init_guess.t, init_guess.lam[5, :], ':', label='Initialized Value', linewidth=linewidth)
ax_lam_psi.set_xlabel(r'Time $t$, [sec]')
ax_lam_psi.set_ylabel(r'$\lambda_\psi$, [1]')
if not INCLUDE_INIT:
    ax_lam_psi.set_ylim((-0.001, 0.001))

ax_alpha.plot(guess.t, guess.u[0, :] * 180 / 3.14159 , label='Guess', linewidth=linewidth)
ax_alpha.plot(sol.t, sol.u[0, :] * 180 / 3.14159 , label='Solution', linewidth=linewidth)
ax_alpha.set_xlabel(r'Time $t$, [sec]')
ax_alpha.set_ylabel(r'AoA $\alpha$, [deg]')

ax_beta.plot(guess.t, guess.u[1, :] * 180 / 3.14159 , label='Guess', linewidth=linewidth)
ax_beta.plot(sol.t, sol.u[1, :] * 180 / 3.14159 , label='Solution', linewidth=linewidth)
ax_beta.set_xlabel(r'Time $t$, [sec]')
ax_beta.set_ylabel(r'Bank $\beta$, [deg]')
ax_beta.set_ylim((-1, 1))

ham_guess = np.array([
    comp_dual.compute_hamiltonian(_t, _x, _lam, _u, guess.p, guess.k)
    for _t, _x, _lam, _u in zip(guess.t, guess.x.T, guess.lam.T, guess.u.T)
])

ham_sol = np.array([
    comp_dual.compute_hamiltonian(_t, _x, _lam, _u, sol.p, sol.k)
    for _t, _x, _lam, _u in zip(sol.t, sol.x.T, sol.lam.T, sol.u.T)
])

ax_ham.plot(guess.t, ham_guess , label='Guess', linewidth=linewidth)
ax_ham.plot(sol.t, ham_sol, label='Solution', linewidth=linewidth)
ax_ham.set_xlabel(r'Time $t$, [sec]')
ax_ham.set_ylabel(r'Hamiltonian $H$, [rad]')
# ax_ham.set_ylim((-90, 90))

handles, labels = ax_lam_h.get_legend_handles_labels()
fig.legend(handles, labels, ncols=3, loc='lower center')

# fig.tight_layout()

fig.subplots_adjust(
        top=0.929,
        bottom=0.094,
        left=0.149,
        right=0.971,
        hspace=0.691,
        wspace=0.737
)

plt.show()
