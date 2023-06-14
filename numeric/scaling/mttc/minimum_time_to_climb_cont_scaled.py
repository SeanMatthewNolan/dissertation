import giuseppe
import numpy as np

from minimum_time_to_climb_prob_scaled import adiff_dual, num_solver

normal_sol_set = giuseppe.load_sol_set('sol_set.data')
sol = normal_sol_set[-1]

max_states = np.max(np.abs(sol.x), axis=1)

# Continuations (from guess BCs to desired BCs)
cont = giuseppe.continuation.ContinuationHandler(num_solver, normal_sol_set[-1])
cont.add_linear_series(100, {
    't_scale': normal_sol_set[-1].t[-1],
#     'h_scale': max_states[0],
#     'v_scale': max_states[1],
#     'gam_scale': max_states[2],
#     'w_scale': max_states[3],
})
sol_set_scaled = cont.run_continuation()

# Save Solution
sol_set_scaled.save('sol_set_scaled.data')

