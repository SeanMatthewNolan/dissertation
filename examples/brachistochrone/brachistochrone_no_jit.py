import os

import numpy as np

import giuseppe.utils.compilation
from giuseppe.continuation import ContinuationHandler
from giuseppe.guess_generation import auto_guess, initialize_guess
from giuseppe.numeric_solvers import SciPySolver
from giuseppe.problems.input import StrInputProb
from giuseppe.problems.symbolic import SymDual
from giuseppe.utils import Timer

giuseppe.utils.compilation.JIT_COMPILE = False

os.chdir(os.path.dirname(__file__))  # Set directory to file location

input_ocp = StrInputProb()

input_ocp.set_independent('t')

input_ocp.add_state('x', 'v*cos(theta)')
input_ocp.add_state('y', 'v*sin(theta)')
input_ocp.add_state('v', '-g*sin(theta)')

input_ocp.add_control('theta')

input_ocp.add_constant('g', 32.2)

input_ocp.add_constant('x_0', 0)
input_ocp.add_constant('y_0', 0)
input_ocp.add_constant('v_0', 1)

input_ocp.add_constant('x_f', 1)
input_ocp.add_constant('y_f', -1)

input_ocp.set_cost('0', '0', 't')

input_ocp.add_constraint('initial', 't')
input_ocp.add_constraint('initial', 'x - x_0')
input_ocp.add_constraint('initial', 'y - y_0')
input_ocp.add_constraint('initial', 'v - v_0')

input_ocp.add_constraint('terminal', 'x - x_f')
input_ocp.add_constraint('terminal', 'y - y_f')

with Timer('Setup Time:'):
    comp_dual = SymDual(input_ocp, control_method='differential')
    solver = SciPySolver(comp_dual)
    guess = auto_guess(comp_dual, u=-15 / 180 * 3.14159)

cont = ContinuationHandler(solver, guess)
cont.add_linear_series(5, {'x_f': 10, 'y_f': -10})
cont.add_linear_series(5, {'v_0': 0})
sol_set = cont.run_continuation()

sol_set.save('sol_set_diff.bin')
