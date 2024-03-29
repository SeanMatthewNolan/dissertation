import copy
import os
import pickle
import math

import numpy as np

from giuseppe.numeric_solvers import SciPySolver
from giuseppe.continuation import ContinuationHandler
from giuseppe.guess_generation import auto_propagate_guess
from giuseppe.problems.input import StrInputProb
from giuseppe.problems.conversions import convert_dual_to_bvp
from giuseppe.problems.symbolic import SymDual, PenaltyConstraintHandler
from giuseppe.utils import Timer
from giuseppe.guess_generation.propagate_guess import propagate_dual_guess_from_guess, propagate_bvp_guess_from_guess
from giuseppe.guess_generation.gauss_newton import match_constants_to_boundary_conditions

os.chdir(os.path.dirname(__file__))  # Set directory to current location

# t_span = 10
t_span = 100
# t_span = 200

min_var = 0.0001

max_var = 10
# max_var = 100

num_cases = 10 * (math.log10(max_var) - math.log10(min_var)) + 1

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

ocp.add_constant('eps_beta', 1e-7)
ocp.add_constant('beta_min', -90 / 180 * 3.1419)
ocp.add_constant('beta_max', 90 / 180 * 3.1419)

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

ocp.add_inequality_constraint('path', 'beta', lower_limit='beta_min', upper_limit='beta_max',
                              regularizer=PenaltyConstraintHandler('eps_beta', method='sec'))


comp_dual = SymDual(ocp, control_method='differential')
solver = SciPySolver(comp_dual)
guess = auto_propagate_guess(comp_dual, control=(15/180*3.14159, 0), t_span=t_span)

nominal_sol = solver.solve(guess)

print(f'\n\nGiuseppe Converged: {nominal_sol.converged} \n\n')

comp_bvp = convert_dual_to_bvp(comp_dual)
# prop_guess = comp_bvp.preprocess_data(copy.deepcopy(nominal_sol))
prop_guess = copy.deepcopy(nominal_sol)

np.random.seed(1927)

beluga_guesses = []
variations = []
bc_reses = []
beluga_attempts = []
converged = []

scale = np.array([])

for _ in range(100):
    for variation in np.geomspace(min_var, max_var, 51):
        variations.append(variation)

        rand = np.random.normal(size=(6, len(prop_guess.t)))
        prop_guess.lam = nominal_sol.lam + rand / np.linalg.norm(rand) * variation
        beluga_guess = propagate_bvp_guess_from_guess(
                comp_bvp, t_span, prop_guess, abs_tol=1e-4, rel_tol=1e-4)
        beluga_guess = match_constants_to_boundary_conditions(
                comp_dual, beluga_guess, abs_tol=1e-6, rel_tol=1e-8)
        beluga_guesses.append(beluga_guess)

        bc_res = np.linalg.norm(
                np.concatenate((
                    comp_dual.compute_initial_adjoint_boundary_conditions(
                            beluga_guess.t[0], beluga_guess.x[:, 0], beluga_guess.lam[:, 0], beluga_guess.u[:, 0],
                            beluga_guess.p, beluga_guess.nu0, beluga_guess.k
                    ),
                    comp_dual.compute_terminal_adjoint_boundary_conditions(
                            beluga_guess.t[-1], beluga_guess.x[:, -1], beluga_guess.lam[:, -1], beluga_guess.u[:, -1],
                            beluga_guess.p, beluga_guess.nuf, beluga_guess.k
                    ),
                ))
        )

        bc_reses.append(bc_res)

        try:
            beluga_attempt = solver.solve(beluga_guess)
            beluga_attempts.append(beluga_attempt)
            if beluga_attempt.converged:
                print(f'Variation: {variation:.4g}; Res: {bc_res:.4g}; Converged')
                converged.append(True)
            else:
                print(f'Variation: {variation:.4g}; Res: {bc_res:.4g}; Not Converged')
                converged.append(False)
        except (ZeroDivisionError, ValueError):
            print(f'Variation: {variation:.4g}; Res: {bc_res:.4g}; Not Converged')
            converged.append(False)

with open(f'beluga_guess_test_{t_span}s.pickle', 'wb') as file:
    pickle.dump((beluga_guesses, variations, bc_reses, converged), file)
