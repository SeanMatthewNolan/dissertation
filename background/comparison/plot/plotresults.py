import json
from beluga.utils import load
import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import numpy as np
from math import sqrt

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


def load_goddard(filename):
    with open(filename, 'r') as file:
        raw_data = json.load(file)

    return raw_data


g_0 = 1.0

h_0 = 1.0
v_0 = 0.0
m_0 = 1.0

t_c = 3.5
h_c = 500
v_c = 620
m_c = 0.6

tar_m_f = m_0 * m_c
c = 0.5 * sqrt(g_0 * h_0)
d_c = 0.5 * v_c * m_0 / g_0
thrust_max = t_c * g_0 * m_0

sol_set_beluga = load('goddard_beluga.json')
sol_beluga = sol_set_beluga[-1][-1]

sol_goddard0 = load_goddard('goddard_0knot.json')
sol_goddard1 = load_goddard('goddard_1knot.json')
sol_goddard2 = load_goddard('goddard_2knot.json')

plt.figure()
plt.plot(sol_beluga.t, sol_beluga.y[:, 0])
plt.plot(sol_goddard0['t'], sol_goddard0['h'])
plt.plot(sol_goddard1['t'], sol_goddard1['h'])
plt.plot(sol_goddard2['t'], sol_goddard2['h'])
plt.xlabel('Time [s]')
plt.ylabel('Altitude [nd]')
plt.title('Altitude Profile')
plt.grid(True)

plt.figure()
plt.plot(sol_beluga.t, sol_beluga.y[:, 1])
plt.xlabel('Time [s]')
plt.ylabel('Velocity [nd]')
plt.title('Velocity Profile')
plt.grid(True)

plt.figure()
plt.plot(sol_beluga.t, sol_beluga.y[:, 2])
plt.xlabel('Time [s]')
plt.ylabel('Mass [nd]')
plt.title('Mass Profile')
plt.grid(True)

plt.figure()
plt.plot(sol_beluga.t, sol_beluga.u[:, 0], label='Thrust')
plt.plot(sol_beluga.t, 1 * d_c * sol_beluga.y[:, 1] ** 2 * np.exp(-h_c * (sol_beluga.y[:, 0] - h_0) / h_0), label='Drag')
plt.xlabel('Time [s]')
plt.ylabel('Force [nd]')
plt.title('Forces')
plt.legend()
plt.grid(True)

fig = plt.figure(figsize=(6.5, 5))
ax1 = fig.add_subplot(211)
plt.plot(sol_goddard0['t'], sol_goddard0['T'], label='OpenGoddard 0 knot (Direct)')
plt.plot(sol_goddard1['t'], sol_goddard1['T'], label='OpenGoddard 1 knot (Direct)')
plt.plot(sol_goddard2['t'], sol_goddard2['T'], label='OpenGoddard 2 knot (Direct)')
plt.plot(sol_beluga.t, sol_beluga.u[:, 0], label='beluga (Indirect)')
plt.xlabel('Time [s]')
plt.ylabel('Thrust [nd]')
plt.title('Thrust (Control)')
fig.suptitle('Comparison Between Direct and Indirect Methods: Goddard Rocket Problem')
plt.legend()

ax21 = fig.add_subplot(223)
plt.plot(sol_goddard0['t'], sol_goddard0['h'], label='OpenGoddard 0 knot')
plt.plot(sol_goddard1['t'], sol_goddard1['h'], label='OpenGoddard 1 knot')
plt.plot(sol_goddard2['t'], sol_goddard2['h'], label='OpenGoddard 2 knot')
plt.plot(sol_beluga.t, sol_beluga.y[:, 0], label='beluga')
plt.xlabel('Time [s]')
plt.ylabel('Altitude [nd]')
plt.title('Altitude')

ax22 = fig.add_subplot(224)
plt.plot(sol_goddard0['t'], sol_goddard0['v'], label='OpenGoddard 0 knot')
plt.plot(sol_goddard1['t'], sol_goddard1['v'], label='OpenGoddard 1 knot')
plt.plot(sol_goddard2['t'], sol_goddard2['v'], label='OpenGoddard 2 knot')
plt.plot(sol_beluga.t, sol_beluga.y[:, 1], label='beluga')
plt.xlabel('Time [s]')
plt.ylabel('Velocity [nd]')
plt.title('Velocity')

plt.tight_layout()

plt.show()
