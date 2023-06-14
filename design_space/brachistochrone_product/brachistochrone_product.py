import os

import numpy as np

from giuseppe.continuation import ContinuationHandler
from giuseppe.guess_generation import auto_guess
from giuseppe.problems.input import StrInputProb
from giuseppe.numeric_solvers import SciPySolver
from giuseppe.problems.symbolic import SymDual
from giuseppe.utils import Timer

os.chdir(os.path.dirname(__file__))  # Set directory to file location

ocp = StrInputProb()

ocp.set_independent('t')

ocp.add_state('x', 'v*cos(theta)')
ocp.add_state('y', 'v*sin(theta)')
ocp.add_state('v', '-g*sin(theta)')

ocp.add_control('theta')

ocp.add_constant('g', 32.2)

ocp.add_constant('x_0', 0)
ocp.add_constant('y_0', 0)
ocp.add_constant('v_0', 0)

ocp.add_constant('x_f', 1)
ocp.add_constant('y_f', -1)

ocp.set_cost('0', '1', '0')

ocp.add_constraint('initial', 't')
ocp.add_constraint('initial', 'x - x_0')
ocp.add_constraint('initial', 'y - y_0')
ocp.add_constraint('initial', 'v - v_0')

ocp.add_constraint('terminal', 'x - x_f')
ocp.add_constraint('terminal', 'y - y_f')

with Timer('Setup Time:'):
    comp_dual = SymDual(ocp, control_method='algebraic')
    num_solver = SciPySolver(comp_dual)
    guess = auto_guess(comp_dual, u=-45/180 * 3.14159)

cont = ContinuationHandler(num_solver, guess)
cont.add_linear_series(5, {'x_f': 10, 'y_f': -10})
sol_set = cont.run_continuation()

cont = ContinuationHandler(num_solver, sol_set[-1])
for row in np.linspace(20, 100, 5):
    cont.add_linear_series(10, {'x_f': 100})
    cont.add_linear_series(1, {'y_f': -row + 10})
    cont.add_linear_series(10, {'x_f': 10})
    cont.add_linear_series(1, {'y_f': -row})
cont.add_linear_series(10, {'x_f': 100})

sol_set = cont.run_continuation()

sol_set.save('sol_set_product.bin')
