import copy

import numpy as np
import casadi as ca

import giuseppe
from giuseppe.problems.automatic_differentiation.utils import lambdify_ca

from lookup_tables import thrust_table_bspline, eta_table_bspline_expanded, CLalpha_table_bspline_expanded,\
    CD0_table_bspline_expanded, temp_table_bspline, dens_table_bspline
from giuseppe.utils.examples.atmosphere1976 import Atmosphere1976

from minimum_time_to_climb_prob_scaled import adiff_dual

sol_set_scaled = giuseppe.load_sol_set('sol_set_scaled.data')

# bc_0 = adiff_dual.ca_initial_boundary_conditions
bc_0 = adiff_dual.ca_initial_adjoint_boundary_conditions
psi_0 = lambdify_ca(bc_0)
ja = ca.jacobian(bc_0, adiff_dual.states)
jac_0 = lambdify_ca(ja)

bc_f = adiff_dual.ca_terminal_adjoint_boundary_conditions
psi_f = lambdify_ca(bc_f)
jac_f = lambdify_ca(bc_f.jacobian())

for sol in sol_set_scaled:
    scales = sol.k[-8:-3]
    print(scales)

    # psi_0 = bc_0(sol.t[0], sol.x[:, 0], sol.p, sol.k)
    psi_0 = bc_0(sol.t[0], sol.x[:, 0], sol.lam[:, 0], sol.u[:, 0], sol.p, sol.nu0, sol.k)
    print(psi_0)

    # j_s = jac_0(sol.t[0], sol.x[:, 0], sol.p, sol.k, psi_0)
    j_s = jac_0(sol.t[0], sol.x[:, 0], sol.lam[:, 0], sol.u[:, 0], sol.p, sol.nu0, sol.k, psi_0)
    print(np.linalg.cond(j_s))
