import numpy as np
import casadi as ca

import giuseppe

from lookup_tables import thrust_table_bspline, eta_table_bspline_expanded, CLalpha_table_bspline_expanded,\
    CD0_table_bspline_expanded, temp_table_bspline, dens_table_bspline
from giuseppe.utils.examples.atmosphere1976 import Atmosphere1976


ocp = giuseppe.problems.automatic_differentiation.ADiffInputProb(dtype=ca.MX)

# Independent Variable
t = ca.MX.sym('t', 1)
ocp.set_independent(t)

# Control
alpha = ca.MX.sym('alpha', 1)
ocp.add_control(alpha)

# Immutable Constant Parameters
Isp = 1600.0
g0 = 32.174
S = 530
mu = 1.4076539e16
Re = 20902900

# States
h = ca.MX.sym('h', 1)
v = ca.MX.sym('v', 1)
gam = ca.MX.sym('gam', 1)
w = ca.MX.sym('w', 1)

# a = 1125.33  # ft/s
# M = v / a  # assume Sea level speed of sound (i.e. constant temperature)
# h_ref = 23800
# rho_0 = 0.002378
# rho = rho_0 * ca.exp(-h / h_ref)

# atm = Atmosphere1976(use_metric=False)
# T, __, rho = atm.get_ca_atm_expr(h)
# a = ca.sqrt(atm.specific_heat_ratio * atm.gas_constant * T)
# M = v/a

atm = Atmosphere1976(use_metric=False)
T = temp_table_bspline(h)
rho = dens_table_bspline(h)
a = ca.sqrt(atm.specific_heat_ratio * atm.gas_constant * T)
a_func = ca.Function('a', (h,), (a,), ('h',), ('a',))
M = v / a

# Look-Up Tables
thrust = thrust_table_bspline(ca.vertcat(M, h))
CLalpha = CLalpha_table_bspline_expanded(M)
CD0 = CD0_table_bspline_expanded(M)
eta = eta_table_bspline_expanded(M)

# Expressions
d2r = ca.pi / 180
r2d = 180 / ca.pi

alpha_hat = alpha
# alpha_hat = alpha * r2d

CD = CD0 + eta * CLalpha * alpha_hat**2
CL = CLalpha * alpha_hat
drag = 0.5 * CD * S * rho * v**2
lift = 0.5 * CL * S * rho * v**2

r = Re + h
mass = w / g0

# Equations of Motion
ocp.add_state(h, v * ca.sin(gam))
ocp.add_state(v, 1/mass * (thrust * ca.cos(alpha) - drag) - mu/r**2 * ca.sin(gam))
ocp.add_state(gam, 1/(mass*v) * (thrust*ca.sin(alpha) + lift) + ca.cos(gam) * (v/r - mu/(v*r**2)))
ocp.add_state(w, -thrust/Isp)

# Inequality Constraints
eps_h = ca.MX.sym('eps_h', 1)

h_min = -10
h_max = 69e3
ocp.add_constant(eps_h, 1e-5)

alpha_max = ca.MX.sym('alpha_max', 1)
eps_alpha = ca.MX.sym('eps_alpha', 1)

ocp.add_constant(alpha_max, 20*d2r)
ocp.add_constant(eps_alpha, 1e-3)

ocp.add_inequality_constraint(
        'path', h,
        lower_limit=h_min, upper_limit=h_max,
        regularizer=giuseppe.problems.automatic_differentiation.regularization.ADiffPenaltyConstraintHandler(
                regulator=eps_h))
ocp.add_inequality_constraint(
        'path', alpha,
        lower_limit=-alpha_max, upper_limit=alpha_max,
        regularizer=giuseppe.problems.automatic_differentiation.regularization.ADiffPenaltyConstraintHandler(
                regulator=eps_alpha))

# Initial Boundary Conditions
t_0 = ca.MX.sym('t_0', 1)
h_0 = ca.MX.sym('h_0', 1)
v_0 = ca.MX.sym('v_0', 1)
gam_0 = ca.MX.sym('gam_0', 1)
w_0 = ca.MX.sym('w_0', 1)

ocp.add_constant(t_0, 0.0)  # s
ocp.add_constant(h_0, 0.0)  # ft
ocp.add_constant(v_0, 424.26)  # ft/s
ocp.add_constant(gam_0, 0.0)  # rad
ocp.add_constant(w_0, 42000.0)  # lb

ocp.add_constraint(location='initial', expr=t - t_0)
ocp.add_constraint(location='initial', expr=h - h_0)
ocp.add_constraint(location='initial', expr=v - v_0)
ocp.add_constraint(location='initial', expr=gam - gam_0)
ocp.add_constraint(location='initial', expr=w - w_0)

# Terminal Boundary Conditions
h_f = ca.MX.sym('h_f', 1)
v_f = ca.MX.sym('v_f', 1)
gam_f = ca.MX.sym('gam_f', 1)

ocp.add_constant(h_f, 65600.0)  # ft
ocp.add_constant(v_f, 968.148)  # ft/s
ocp.add_constant(gam_f, 0.0)  # rad

ocp.add_constraint(location='terminal', expr=h - h_f)
ocp.add_constraint(location='terminal', expr=v - v_f)
ocp.add_constraint(location='terminal', expr=gam - gam_f)

# Objective Function
ocp.set_cost(0, 0, t)

# Compilation
with giuseppe.utils.Timer(prefix='Compilation Time:'):
    adiff_dual = giuseppe.problems.automatic_differentiation.ADiffDual(ocp)
    num_solver = giuseppe.numeric_solvers.SciPySolver(adiff_dual)

if __name__ == '__main__':
    # Guess Generation (overwrites the terminal conditions in order to converge)
    guess = giuseppe.guess_generation.auto_propagate_guess(adiff_dual, control=6*d2r, t_span=1)

    # Continuations (from guess BCs to desired BCs)
    cont = giuseppe.continuation.ContinuationHandler(num_solver, guess)
    cont.add_linear_series(50, {'h_f': 10_000, 'v_f': 1_000, 'gam_f': 35 * np.pi/180})
    cont.add_linear_series(50, {'h_f': 65_600.0, 'v_f': 968.148})
    cont.add_linear_series(50, {'gam_f': 0})
    # cont.add_logarithmic_series(50, {'eps_alpha': 1e-5, 'eps_h': 1e-5})
    sol_set = cont.run_continuation()

    # Save Solution
    sol_set.save('sol_set.data')
