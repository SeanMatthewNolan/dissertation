import itertools
import os

import sympy
import numpy as np

from giuseppe import ContinuationHandler, auto_propagate_guess, StrInputProb, SciPySolver, SymDual, Timer,\
    SymPenaltyConstraintHandler, propagate_guess

os.chdir(os.path.dirname(__file__))  # Set directory to file location

x, y, theta, v, omega = sympy.symbols('x, y, theta, v, omega')
x_b, y_b = sympy.symbols('x_b, y_b')
sig_v, sig_omega, sig_rho = sympy.symbols('sig_v, sig_omega, sig_rho')

x_vec = sympy.Matrix([[x, y, theta]]).T
f_vec = sympy.Matrix([[v * sympy.cos(theta), v * sympy.sin(theta), omega]]).T

rho = sympy.sqrt((x - x_b) ** 2 + (y - y_b) ** 2)
h = sympy.Matrix([[rho]])

f_mat = f_vec.jacobian(x_vec)
h_mat = h.jacobian(x_vec)

q_mat = sympy.Matrix(
        [[sig_v**2, 0],
         [0,        sig_omega**2]]
)

r_mat = sympy.Matrix(
        [[sig_rho**2]]
)

g_mat = sympy.Matrix(
        [
            [sympy.cos(theta), 0],
            [sympy.sin(theta), 0],
            [0, 1]
        ]
)

p_11, p_12, p_13, p_22, p_23, p_33 = sympy.symbols('p_11, p_12, p_13, p_22, p_23, p_33')

p_mat = sympy.Matrix(
        [
            [p_11, p_12, p_13],
            [p_12, p_22, p_23],
            [p_13, p_23, p_33]
        ]
)

p_dot = f_mat @ p_mat + p_mat @ f_mat.T - p_mat @ h_mat.T @ r_mat.inv() @ h_mat @ p_mat + g_mat @ q_mat @ g_mat.T


input_ocp = StrInputProb()

input_ocp.set_independent('t')

input_ocp.add_state('x', 'v * cos(theta)')
input_ocp.add_state('y', 'v * sin(theta)')
input_ocp.add_state('theta', 'omega')

input_ocp.add_control('omega')

input_ocp.add_constant('v', 10)

input_ocp.add_constant('x_b', 500)
input_ocp.add_constant('y_b', 250)

input_ocp.add_constant('sig_v', 0.3)
input_ocp.add_constant('sig_omega', 0.1)
input_ocp.add_constant('sig_rho', 1)

input_ocp.add_constant('omega_max', 0.1)
input_ocp.add_constant('eps_omega', 1)

input_ocp.add_constant('k_mit', 1e-2)

input_ocp.add_constant('x_0', 0)
input_ocp.add_constant('y_0', 0)
input_ocp.add_constant('theta_0', 0)

input_ocp.add_constant('x_f', 1)
input_ocp.add_constant('y_f', 0)

input_ocp.set_cost('0', 'p_11 + p_22', '0')
# input_ocp.set_cost('0', 'k_mit / 2 * omega**2', 'p_11 + p_22')

input_ocp.add_constraint('initial', 't')
input_ocp.add_constraint('initial', 'x - x_0')
input_ocp.add_constraint('initial', 'y - y_0')
input_ocp.add_constraint('initial', 'theta - theta_0')

input_ocp.add_constraint('terminal', 'x - x_f')
input_ocp.add_constraint('terminal', 'y - y_f')

for i, j in itertools.combinations_with_replacement(range(3), 2):
    input_ocp.add_state(str(p_mat[i, j]), str(p_dot[i, j]))
    input_ocp.add_constraint('initial', str(p_mat[i, j]))


input_ocp.add_inequality_constraint(
        'path', 'omega', lower_limit='-omega_max', upper_limit='omega_max',
        regularizer=SymPenaltyConstraintHandler('eps_omega'))

with Timer('Setup Time:'):
    comp_dual = SymDual(input_ocp, control_method='differential')
    comp_dual.add_post_process('compute_hamiltonian')
    comp_dual.add_post_process('compute_huu')
    solver = SciPySolver(comp_dual)

guess = propagate_guess(comp_dual, control=0 / 180 * 3.14159, t_span=100, initial_states=np.zeros((9,)))

guess.save('nav_straight.bin')
