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

fig = plt.figure(figsize=(3.5, 2.75))
ax1 = fig.add_subplot(111)
plt.plot(sol_goddard0['t'], sol_goddard0['T'], label='Direct 0 knot')
plt.plot(sol_goddard1['t'], sol_goddard1['T'], label='Direct 1 knot')
plt.plot(sol_goddard2['t'], sol_goddard2['T'], label='Direct 2 knot')
plt.plot(sol_beluga.t, sol_beluga.u[:, 0], label='Indirect')
plt.xlabel('Time [s]')
plt.ylabel('Thrust [nd]')
# plt.title('Thrust (Control)')
# fig.suptitle('Goddard Rocket Problem:\nComparison Between Direct\nand Indirect Methods')
fig.suptitle('Goddard Rocket Problem')
fig.legend(ncols=2, loc='lower center')

# ax21 = fig.add_subplot(122)
# plt.plot(sol_goddard0['t'], sol_goddard0['h'], label='OpenGoddard 0 knot')
# plt.plot(sol_goddard1['t'], sol_goddard1['h'], label='OpenGoddard 1 knot')
# plt.plot(sol_goddard2['t'], sol_goddard2['h'], label='OpenGoddard 2 knot')
# plt.plot(sol_beluga.t, sol_beluga.y[:, 0], label='beluga')
# plt.xlabel('Time [s]')
# plt.ylabel('Altitude [nd]')
# # plt.title('Altitude')

# ax22 = fig.add_subplot(224)
# plt.plot(sol_goddard0['t'], sol_goddard0['v'], label='OpenGoddard 0 knot')
# plt.plot(sol_goddard1['t'], sol_goddard1['v'], label='OpenGoddard 1 knot')
# plt.plot(sol_goddard2['t'], sol_goddard2['v'], label='OpenGoddard 2 knot')
# plt.plot(sol_beluga.t, sol_beluga.y[:, 1], label='beluga')
# plt.xlabel('Time [s]')
# plt.ylabel('Velocity [nd]')
# plt.title('Velocity')

# plt.tight_layout()

fig.subplots_adjust(
        top=0.847,
        bottom=0.38,
        left=0.146,
        right=0.957,
        hspace=0.2,
        wspace=0.508
)

plt.show()
