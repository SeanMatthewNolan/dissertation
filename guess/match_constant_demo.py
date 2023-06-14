import copy
import os
import pickle
import math

import numpy as np

from giuseppe.numeric_solvers import SciPySolver
from giuseppe.guess_generation import auto_propagate_guess, initialize_guess_from_partial_solution
from giuseppe.guess_generation.gauss_newton import match_constants_to_boundary_conditions, gauss_newton
from giuseppe.problems.input import StrInputProb
from giuseppe.problems.symbolic import SymDual, PenaltyConstraintHandler, SymOCP

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

ocp.add_constant('eps_alpha', 1e-4)
ocp.add_constant('alpha_min', -80 / 180 * 3.1419)
ocp.add_constant('alpha_max', 80 / 180 * 3.1419)

ocp.add_constant('eps_beta', 1e-4)
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


comp_ocp = SymOCP(ocp)
comp_dual = SymDual(ocp, control_method='differential')
solver = SciPySolver(comp_dual)
guess_ocp = auto_propagate_guess(comp_ocp, control=(15/180*3.14159, 0), t_span=100, match_constants=False)
guess = auto_propagate_guess(comp_dual, control=(15/180*3.14159, 0), t_span=100, match_constants=False)

kuldge_guess = initialize_guess_from_partial_solution(comp_dual, guess_ocp)


def match_constants_to_dynamics(_guess):
    _t, _x, _u, _p = _guess.t, _guess.x, _guess.u, _guess.p

    _compute_dynamics = comp_dual.compute_dynamics
    _compute_boundary_conditions = comp_dual.compute_boundary_conditions

    _h_arr = np.diff(_t)
    _t_mid = (_t[:-1] + _t[1:]) / 2
    _u_mid = (_u[:, :-1] + _u[:, 1:]) / 2

    def _fitting_function(_k: np.ndarray) -> np.ndarray:
        res_bc = _compute_boundary_conditions(_t, _x, _p, _k)

        _x_dot_nodes = np.array([
            _compute_dynamics(_t_i, _x_i, _u_i, _p, _k)
            for _t_i, _x_i, _u_i in zip(_t, _x.T, _u.T)
        ]).T

        _x_dot_diff = np.diff(_x_dot_nodes)
        _x_mid = (_x[:, :-1] + _x[:, 1:]) / 2 - _h_arr / 8 * np.diff(_x_dot_nodes)

        _x_dot_mid = np.array([
            _compute_dynamics(_t_i, _x_i, _u_i, _p, _k)
            for _t_i, _x_i, _u_i in zip(_t_mid, _x_mid.T, _u_mid.T)
        ]).T

        _delta_x = np.diff(_x)

        res_dyn = _delta_x - _h_arr * ((_x_dot_nodes[:, :-1] + _x_dot_nodes[:, 1:]) / 6 + 2 / 3 * _x_dot_mid)

        return np.concatenate((res_bc, res_dyn.flatten()))

    return gauss_newton(_fitting_function, guess.k, verbose=True, max_steps=500)


out_str = ''
for constant, default, state_matched, dual_matched, dyna_matched \
        in zip(comp_dual.annotations.constants, comp_dual.default_values,
               match_constants_to_boundary_conditions(comp_dual, kuldge_guess, verbose=True).k,
               match_constants_to_boundary_conditions(comp_dual, kuldge_guess, use_adjoint_bcs=True, verbose=True).k,
               match_constants_to_dynamics(kuldge_guess)
               ):
    # out_str += f'{constant:10} & {default:12.6g} & {state_matched:12.6g} & {dual_matched:12.6g} \\\\\n'
    out_str += f' & {default:12.3g} & {state_matched:12.3g} & {dual_matched:12.3g} & {dyna_matched:12.3g} \\\\\n'

print(out_str)
