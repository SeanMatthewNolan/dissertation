# -*- coding: utf-8 -*-
# Copyright 2017 Interstellar Technologies Inc. All Rights Reserved.

from __future__ import print_function
import json
import numpy as np
import matplotlib.pyplot as plt
from OpenGoddard.optimize import Problem, Guess, Condition, Dynamics


class Rocket:
    g0 = 1.0  # Gravity at surface [-]

    def __init__(self):
        self.H0 = 1.0  # Initial height
        self.V0 = 0.0  # Initial velocity
        self.M0 = 1.0  # Initial mass
        self.Tc = 3.5  # Use for thrust
        self.Hc = 500  # Use for drag
        self.Vc = 620  # Use for drag
        self.Mc = 0.6  # Fraction of initial mass left at end
        self.c = 0.5 * np.sqrt(self.g0*self.H0)  # Thrust-to-fuel mass
        self.Mf = self.Mc * self.M0               # Final mass
        self.Dc = 0.5 * self.Vc * self.M0 / self.g0  # Drag scaling
        self.T_max = self.Tc * self.g0 * self.M0     # Maximum thrust


def dynamics(prob, obj, section):
    h = prob.states(0, section)
    v = prob.states(1, section)
    m = prob.states(2, section)
    T = prob.controls(0, section)

    Dc = obj.Dc
    c = obj.c
    drag = 1 * Dc * v ** 2 * np.exp(-obj.Hc * (h - obj.H0) / obj.H0)
    g = obj.g0 * (obj.H0 / h)**2

    dx = Dynamics(prob, section)
    dx[0] = v
    dx[1] = (T - drag) / m - g
    dx[2] = - T / c
    return dx()


def equality(prob, obj):
    h = prob.states_all_section(0)
    v = prob.states_all_section(1)
    m = prob.states_all_section(2)
    T = prob.controls_all_section(0)
    tf = prob.time_final(-1)

    ts0 = prob.time_final(0)
    T1 = prob.controls(0, 1)

    result = Condition()

    # event condition
    result.equal(h[0], obj.H0)
    result.equal(v[0], obj.V0)
    result.equal(m[0], obj.M0)
    result.equal(v[-1], 0.0)
    result.equal(m[-1], obj.Mf)

    result.equal(ts0, 0.075)

    return result()


def inequality(prob, obj):
    h = prob.states_all_section(0)
    v = prob.states_all_section(1)
    m = prob.states_all_section(2)
    T = prob.controls_all_section(0)
    ts0 = prob.time_final(0)
    tf = prob.time_final(-1)

    result = Condition()
    # lower bounds
    result.lower_bound(h, obj.H0)
    result.lower_bound(v, 0.0)
    result.lower_bound(m, obj.Mf)
    result.lower_bound(T, 0.0)
    result.lower_bound(tf, 0.1)
    result.lower_bound(ts0, 0.05)
    # upper bounds
    result.upper_bound(m, obj.M0)
    result.upper_bound(T, obj.T_max)

    return result()


def cost(prob, obj):
    h = prob.states_all_section(-1)
    return -h[-1]


def cost_derivative(prob, obj):
    jac = Condition(prob.number_of_variables)
    index_h_end = prob.index_states(0, -1, -1)
    jac.change_value(index_h_end, -1)
    return jac()

# ========================
plt.close("all")
# Program Starting Point
time_init = [0.0, 0.1, 0.15, 0.3]
n = [25, 25, 25]
num_states = [3, 3, 3]
num_controls = [1, 1, 1]
max_iteration = 50

flag_savefig = True
savefig_file = "04_Goddard/05_2knot_"

# ------------------------
# set OpenGoddard class for algorithm determination
prob = Problem(time_init, n, num_states, num_controls, max_iteration)

# ------------------------
# create instance of operating object
obj = Rocket()

prob.set_unit_states_all_section(0, 0.1)

# ========================
# Initial parameter guess

# altitude profile
H_init = Guess.cubic(prob.time_all_section, 1.0, 0.0, 1.010, 0.0)
# Guess.plot(prob.time_all_section, H_init, "Altitude", "time", "Altitude")
# if(flag_savefig):plt.savefig(savefig_file + "guess_alt" + ".png")

# velocity
V_init = Guess.linear(prob.time_all_section, 0.0, 0.0)
# Guess.plot(prob.time_all_section, V_init, "Velocity", "time", "Velocity")

# mass profile
M_init0 = Guess.linear(prob.time[0], 1.0, 0.8)
M_init1 = Guess.linear(prob.time[1], 0.8, 0.6)
M_init2 = Guess.linear(prob.time[2], 0.6, 0.6)
M_init = np.hstack((M_init0, M_init1, M_init2))
# Guess.plot(prob.time_all_section, M_init, "Mass", "time", "Mass")
# if(flag_savefig):plt.savefig(savefig_file + "guess_mass" + ".png")

# thrust profile
T_init0 = Guess.linear(prob.time[0], 3.5, 3.5)
T_init1 = Guess.linear(prob.time[1], 3.0, 3.0)
T_init2 = Guess.linear(prob.time[2], 0.0, 0.0)
T_init = np.hstack((T_init0, T_init1, T_init2))
# Guess.plot(prob.time_all_section, T_init, "Thrust Guess", "time", "Thrust")
# if(flag_savefig):plt.savefig(savefig_file + "guess_thrust" + ".png")

plt.show()

# ========================
# Substitution initial value to parameter vector to be optimized
prob.set_states_all_section(0, H_init)
prob.set_states_all_section(1, V_init)
prob.set_states_all_section(2, M_init)
prob.set_controls_all_section(0, T_init)

# ========================
# Main Process
# Assign problem to SQP solver
prob.dynamics = [dynamics, dynamics, dynamics]
prob.knot_states_smooth = [True, True]
prob.cost = cost
prob.cost_derivative = cost_derivative
prob.equality = equality
prob.inequality = inequality


def display_func():
    h = prob.states_all_section(0)
    print("max altitude: {0:.5f}".format(h[-1]))

prob.solve(obj, display_func, ftol=1e-10)

# ========================
# Post Process
# ------------------------
# Convert parameter vector to variable
h = prob.states_all_section(0)
v = prob.states_all_section(1)
m = prob.states_all_section(2)
T = prob.controls_all_section(0)
time = prob.time_update()

# ------------------------
# Calculate necessary variables
Dc = 0.5 * 620 * 1.0 / 1.0
drag = 1 * Dc * v ** 2 * np.exp(-500 * (h - 1.0) / 1.0)
g = 1.0 * (1.0 / h)**2

data_out = {
    'h': list(h),
    'v': list(v),
    'm': list(m),
    'T': list(T),
    't': list(time),
    'Dc': Dc,
    'drag': list(drag),
    'g': list(g)
}

with open('goddard_2knot.json', 'w') as fp:
    json.dump(data_out, fp)
