import os

import numpy as np

from scipy.interpolate import PchipInterpolator
from giuseppe import ContinuationHandler, StrInputProb, SymDual, auto_propagate_guess, \
    SymPenaltyConstraintHandler, SciPySolver, load_sol

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

ocp.add_constant('c_l0', -0.2317)
ocp.add_constant('c_l1', 0.0513*180/3.14159)
ocp.add_constant('c_l2', 0.2945)
ocp.add_constant('c_l3', -0.1028)

ocp.add_constant('c_d0', 0.024)
ocp.add_constant('c_d1', 7.24e-4*180**2/3.14159**2)
ocp.add_constant('c_d2', 0.406)
ocp.add_constant('c_d3', -0.323)

# T = g0 * h_ref / R
temp = G0 * 23_800 / 1716
# a = sqrt(gamma * R * T)
a = (1.4 * 1716 * temp)**0.5

ocp.add_constant('a', a)
ocp.add_constant('s_ref', 750 / 144)
ocp.add_constant('m', 2_000 / G0)

ocp.add_constant('eps_alpha', 1e-5)
ocp.add_constant('alpha_max', 30 / 180 * 3.1419)

ocp.add_constant('eps_beta', 1e-5)
ocp.add_constant('beta_max', 180 / 180 * 3.1419)

ocp.add_constant('xi', 0)

ocp.add_constant('h_0', 30_000)  # 165_000
ocp.add_constant('phi_0', 0)
ocp.add_constant('theta_0', 0)
ocp.add_constant('v_0', 900)   # 20_500
ocp.add_constant('gamma_0', -1 / 180 * np.pi)
ocp.add_constant('psi_0', np.pi / 2)

ocp.add_constant('h_f', 10_000)

ocp.set_cost('0', '0', '-phi * cos(xi) - theta  * sin(xi)')

ocp.add_constraint('initial', 't')
ocp.add_constraint('initial', 'h - h_0')
ocp.add_constraint('initial', 'phi - phi_0')
ocp.add_constraint('initial', 'theta - theta_0')
ocp.add_constraint('initial', 'v - v_0')
ocp.add_constraint('initial', 'gamma - gamma_0')
ocp.add_constraint('initial', 'psi - psi_0')

ocp.add_constraint('terminal', 'h - h_f')

ocp.add_inequality_constraint('path', 'alpha', lower_limit='-alpha_max', upper_limit='alpha_max',
                              regularizer=SymPenaltyConstraintHandler('eps_alpha'))
ocp.add_inequality_constraint('path', 'beta', lower_limit='-beta_max', upper_limit='beta_max',
                              regularizer=SymPenaltyConstraintHandler('eps_beta'))

comp_prob = SymDual(ocp, control_method='differential')
num_solver = SciPySolver(comp_prob)

guess = auto_propagate_guess(comp_prob, control=(5 / 180 * 3.14159, 0), t_span=100)

cont = ContinuationHandler(num_solver, guess)
cont.add_linear_series(100, {'h_0': 165_000, 'v_0': 20_000})
cont.add_linear_series(100, {'phi_0': 0, 'theta_0': 0, 'gamma_0': -5 / 180 * np.pi, 'psi_0': np.pi/2})
cont.add_linear_series(100, {'h_f': 0})
sol_set = cont.run_continuation()

root = sol_set[-1]

ref_traj = load_sol('ref_traj.bin')
ref_state_inter = PchipInterpolator(ref_traj.t, ref_traj.x, axis=1)

for ti in np.arange(0, ref_traj.t[-1], 50):
    print(f'{ti}/{ref_traj.t[-1]}')
    ref_states = ref_state_inter(ti)
    cont_states = ContinuationHandler(num_solver, root)
    cont_states.add_linear_series(
            100,
            {'h_0': ref_states[0], 'phi_0': ref_states[1], 'theta_0': ref_states[2],
             'v_0': ref_states[3], 'gamma_0': ref_states[4], 'psi_0': ref_states[5]}
    )
    cont_states.add_linear_series(10, {'xi': ref_states[5] - np.pi / 2})
    sol_set_states = cont_states.run_continuation()
    root = sol_set_states[-1]

    cont_left = ContinuationHandler(num_solver, root)
    cont_left.add_linear_series(200, {'xi': 200 / 180 * np.pi - (ref_states[5] - np.pi / 2)})
    sol_set_left = cont_left.run_continuation()

    cont_footprint = ContinuationHandler(num_solver, sol_set_left[-1])
    cont_footprint.add_linear_series(400, {'xi': -200 / 180 * np.pi - (ref_states[5] - np.pi / 2)})
    sol_set_footprint = cont_footprint.run_continuation()

    sol_set_footprint.save(f'footprint_t{round(ti)}.bin')
