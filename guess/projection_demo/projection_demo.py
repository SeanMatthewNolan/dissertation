import numpy as np
import matplotlib.pyplot as plt
from giuseppe.guess_generation.gauss_newton import gauss_newton


def constraint(arr):
    return np.expand_dims(arr[0]**2 + (arr[1] / 2)**2 + arr[2] ** 2 - 1, axis=0)


def goursat_tangle(arr):
    x, y, z = arr
    a, b, c = 0.0, -5.0, 11.8
    return x**4+y**4+z**4+a*(x**2+y**2+z**2)**2+b*(x**2+y**2+z**2)+c


theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
z = np.linspace(-2, 2, 100)
r = z**2 + 1
x = r * np.sin(theta)
y = r * np.cos(theta)

states = np.array((x, y, z))
proj_states = np.empty_like(states)

for idx, states_i in enumerate(states.T):
    proj_states[:, idx] = gauss_newton(constraint, states_i)

fig1 = plt.figure(figsize=(6.5, 5))
ax111 = fig1.add_subplot(111, projection='3d')

ax111.plot(states[0, :], states[1, :], states[2, :])
ax111.plot(proj_states[0, :], proj_states[1, :], proj_states[2, :])

plt.show()
