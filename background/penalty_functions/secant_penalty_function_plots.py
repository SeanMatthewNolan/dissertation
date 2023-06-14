import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams, rc

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

x_min, x_max = -1, 1
x_mid = (x_max + x_min)/2


def utm_penalty(_x):
    return 1/np.cos(np.pi/2 * (_x - x_mid)/(x_max - x_mid))


def utm_penalty_offset(_x):
    return 1/np.cos(np.pi/2 * (_x - x_mid)/(x_max - x_mid)) - 1


# size = (3.5, 3.5)
size = (4, 3.5)

x = np.linspace(x_min, x_max, 100)
y_utm = utm_penalty(x)
y_utm_off = utm_penalty_offset(x)

yticks_val = [1, 5, 10]
# yticks_lab = [r'${0:d}\epsilon_j$'.format(_val) for _val in yticks_val]
yticks_lab = [r'${0:d}\epsilon_i$'.format(_val) for _val in yticks_val]

fig = plt.figure(0, figsize=size)
ax = plt.gca()
ax.set_title('Secant Penalty Functions')
# plt.xticks(np.array([x_min, x_mid, x_max]), (r'$S_j^{min}$', r'$(S_j^{min} + S_j^{max})/2$', r'$S_j^{max}$'))
plt.xticks(np.array([x_min, x_mid, x_max]),
           (r'$C_{i, \text{min}}$', r'$(C_{i, \text{min}} + C_{i, \text{max}})/2$', r'$C_{i, \text{max}}$'))
plt.yticks(yticks_val, yticks_lab)
ax.plot(x, y_utm,  '--', label='Original UTM Penalty Function')
ax.plot(x, y_utm_off,  label='Offset Penalty Function')
ax.set_xlabel(r'$C_i[t, \bm{x}{(t)}, \bm{u}{(t)}]$')
# ax.set_xlabel(r'$S_j[t, x(t), u(t)]$')
ax.set_ylabel(r'Penalty Cost')
# fig.subplots_adjust(bottom=.30, left=.18, right=1-.18)
ax.axis([x_min, x_max, 0, yticks_val[-1]*1.10])
ax.legend()
fig.subplots_adjust(
        top=0.9,
        bottom=0.167,
        left=0.165,
        right=0.916,
        hspace=0.2,
        wspace=0.2
)

# fig.savefig('utm_penalty.png', dpi=300)

plt.show()
