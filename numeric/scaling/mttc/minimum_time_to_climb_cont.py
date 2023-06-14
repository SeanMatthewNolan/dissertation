import numpy as np
import casadi as ca

import giuseppe

from minimum_time_to_climb_prob import adiff_dual, num_solver

# Guess Generation (overwrites the terminal conditions in order to converge)
guess = giuseppe.guess_generation.auto_propagate_guess(adiff_dual, control=6/180 * 3.14159, t_span=1)
# guess = giuseppe.InteractiveGuessGenerator(adiff_dual, num_solver=num_solver, init_guess=guess).run()

# Continuations (from guess BCs to desired BCs)
cont = giuseppe.continuation.ContinuationHandler(num_solver, guess)
cont.add_linear_series(50, {'h_f': 10_000, 'v_f': 1_000, 'gam_f': 35 * np.pi/180})
cont.add_linear_series(50, {'h_f': 65_600.0, 'v_f': 968.148})
cont.add_linear_series(50, {'gam_f': 0})
sol_set = cont.run_continuation()

# Save Solution
sol_set.save('sol_set.data')
