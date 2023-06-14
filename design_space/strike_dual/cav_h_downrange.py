import os

import numpy as np

from giuseppe import ContinuationHandler, InteractiveGuessGenerator, StrInputProb, SymDual, auto_propagate_guess, \
    SymPenaltyConstraintHandler, SciPySolver

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
ocp.add_state('phi', 'v * cos(gamma) / r')
ocp.add_state('v', '-drag / m - g * sin(gamma)')
ocp.add_state('gamma', 'lift / (m * v) + cos(gamma) * (v / r - g / v)')

ocp.add_control('alpha')

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

ocp.add_constant('eps_h', 1e-6)
ocp.add_constant('h_max', 210_000)

ocp.add_constant('eps_alpha', 1e-5)
ocp.add_constant('alpha_max', 30 / 180 * 3.1419)

ocp.add_constant('h_0', 30_000)  # 165_000
ocp.add_constant('phi_0', 0)
ocp.add_constant('v_0', 900)   # 12_500
ocp.add_constant('gamma_0', -1 / 180 * np.pi)

ocp.add_constant('h_f', 10_000)
ocp.add_constant('v_f', 900)

ocp.set_cost('0', '0', '-phi')

ocp.add_constraint('initial', 't')
ocp.add_constraint('initial', 'h - h_0')
ocp.add_constraint('initial', 'phi - phi_0')
ocp.add_constraint('initial', 'v - v_0')
ocp.add_constraint('initial', 'gamma - gamma_0')

ocp.add_constraint('terminal', 'h - h_f')

ocp.add_inequality_constraint('path', 'h', lower_limit='-h_max', upper_limit='h_max',
                              regularizer=SymPenaltyConstraintHandler('eps_h'))
ocp.add_inequality_constraint('path', 'alpha', lower_limit='-alpha_max', upper_limit='alpha_max',
                              regularizer=SymPenaltyConstraintHandler('eps_alpha'))

sym_ocp = SymDual(ocp, control_method='differential')
comp_prob = sym_ocp.compile()
num_solver = SciPySolver(comp_prob)

guess = auto_propagate_guess(comp_prob, control=(7.5 / 180 * 3.14159,), t_span=100, max_step=5)

cont = ContinuationHandler(num_solver, guess)
cont.add_linear_series(100, {'h_0': 165_000, 'v_0': 20_000})
cont.add_linear_series(100, {'h_f': 0})
cont.add_linear_series(20, {'alpha_max': 20 / 180 * 3.14159})
cont.add_linear_series(50, {'gamma_0': -5 / 180 * 3.14159})
sol_set = cont.run_continuation()

sol_set[-1].save('downrange_unconstrained.bin')
