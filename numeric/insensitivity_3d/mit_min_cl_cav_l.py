import os
from math import cos, sin

import numpy as np

from giuseppe.continuation import ContinuationHandler
from giuseppe.guess_generation import auto_propagate_guess
from giuseppe.problems.input import StrInputProb
from giuseppe.problems.symbolic import SymDual
from giuseppe.problems.symbolic.regularization import PenaltyConstraintHandler
from giuseppe.numeric_solvers.bvp import SciPySolver

os.chdir(os.path.dirname(__file__))

RE = 20_902_900
G0 = 32.17404856

ocp = StrInputProb()

ocp.set_independent('t')

ocp.add_expression('r', 're + h')
ocp.add_expression('g', 'mu / r**2')
ocp.add_expression('rho', 'rho_0 * exp(-h / h_ref)')
ocp.add_expression('dyn_pres', '1 / 2 * rho * v ** 2 ')
ocp.add_expression('lift', 'c_l * s_ref * dyn_pres')
ocp.add_expression('drag', 'c_d * s_ref * dyn_pres')
ocp.add_expression('c_l', 'c_l1 * alpha + c_l2 * exp(c_l3 * v / a) + c_l0')
ocp.add_expression('c_d', 'c_d1 * alpha **2 + c_d2 * exp(c_d3 * v / a) + c_d0')

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
ocp.add_constant('mu', 0.14076539e17)

ocp.add_constant('c_l0', -0.1177)
ocp.add_constant('c_l1', 0.0366*180/3.14159)
ocp.add_constant('c_l2', 0.2664)
ocp.add_constant('c_l3', -0.20)

ocp.add_constant('c_d0', 0.0442)
ocp.add_constant('c_d1', 7.53e-4*180**2/3.14159**2)
ocp.add_constant('c_d2', 0.27)
ocp.add_constant('c_d3', -0.407)

# T = g0 * h_ref / R
temp = G0 * 23_800 / 1716
# a = sqrt(gamma * R * T)
a = (1.4 * 1716 * temp)**0.5

ocp.add_constant('a', a)
ocp.add_constant('s_ref', 500/144)
ocp.add_constant('m', 1800 / G0)

ocp.add_constant('angle_scale', 1)
ocp.add_constant('c_l_scale', 0.5)
ocp.add_constant('k_mit', 100)

ocp.add_constant('eps_h', 1e-4)
ocp.add_constant('h_max', 200_000)

ocp.add_constant('eps_alpha', 1e-4)
ocp.add_constant('alpha_max', 25 / 180 * 3.1419)

ocp.add_constant('eps_beta', 1e-4)
ocp.add_constant('beta_max', 89 / 180 * 3.1419)

ocp.add_constant('h_0', 165_000)
ocp.add_constant('phi_0', 0)
ocp.add_constant('theta_0', 0)
ocp.add_constant('v_0', 20_000)
ocp.add_constant('gamma_0', -1 / 180 * np.pi)
ocp.add_constant('psi_0', np.pi / 2)

ocp.add_constant('h_f', 10_000)
ocp.add_constant('phi_f', 0.01)
ocp.add_constant('theta_f', 0.0)

ocp.set_cost('0', '(lift) ** 2 + k_mit / 2 * beta ** 2', '0')
# ocp.set_cost('0', 'k_mit / 2 * beta ** 2', 't')

ocp.add_constraint('initial', 't')
ocp.add_constraint('initial', '(h - h_0)/h_0')
ocp.add_constraint('initial', '(phi - phi_0)/angle_scale')
ocp.add_constraint('initial', '(theta - theta_0)/angle_scale')
ocp.add_constraint('initial', '(v - v_0)/v_0')
ocp.add_constraint('initial', '(gamma - gamma_0)/angle_scale')
ocp.add_constraint('initial', '(psi - psi_0)/angle_scale')

ocp.add_constraint('terminal', '(h - h_f)/h_0')
ocp.add_constraint('terminal', '(phi - phi_f)/angle_scale')
ocp.add_constraint('terminal', '(theta - theta_f)/angle_scale')

# ocp.add_inequality_constraint('path', 'h', lower_limit='-h_max', upper_limit='h_max',
#                               regularizer=PenaltyConstraintHandler('eps_h'))
# ocp.add_inequality_constraint('path', 'alpha', lower_limit='-alpha_max', upper_limit='alpha_max',
#                               regularizer=PenaltyConstraintHandler('eps_alpha'))
# ocp.add_inequality_constraint('path', 'beta', lower_limit='-beta_max', upper_limit='beta_max',
#                               regularizer=PenaltyConstraintHandler('eps_beta'))

sym_ocp = SymDual(ocp, control_method='differential')
num_solver = SciPySolver(sym_ocp, verbose=2)

sym_ocp.add_post_process('compute_hamiltonian')
sym_ocp.add_post_process('compute_huu')

guess = auto_propagate_guess(sym_ocp, control=(1.5 / 180 * 3.14159, 0 / 180 * 3.14159), t_span=25)

cont = ContinuationHandler(num_solver, guess)
cont.add_linear_series(100, {'h_f': 0, 'phi_f': 40 / 180 * np.pi})
cont.add_linear_series(100, {'theta_f': 10 / 180 * np.pi})
# cont.add_linear_series(100, {'k_mit': 0})
sol_set = cont.run_continuation()

sol_set.save('min_cl.bin')
