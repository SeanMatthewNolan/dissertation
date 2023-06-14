import matplotlib.pyplot as plt
import numpy as np

import giuseppe as gp
from matplotlib import rc, rcParams, gridspec

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

cmap = plt.colormaps['viridis']

sol_set = gp.data_classes.load_sol_set('crossrange_continuation_scaled.data')

fig_cont = plt.figure(figsize=(3, 1.25))
ax_cont = fig_cont.add_subplot()

xf, yf = [], []
for sol in sol_set:
    phi, theta = sol.x[1, :] * 180 / np.pi, sol.x[2, :] * 180 / np.pi
    ax_cont.plot(phi, theta, color='C0')
    xf.append(phi[-1])
    yf.append(theta[-1])
ax_cont.plot(xf, yf, '--', label='_', color='C1')

skip = 10
offset, mult = 5, 0.95
for xfl, xfr, yfl, yfr in zip(xf[::skip], xf[skip - 1::skip], yf[:-1:skip], yf[skip - 1::skip]):
    offset_theta = np.arctan2(xfr - xfl, yfr - yfl)
    offset_x = np.cos(offset_theta) * offset
    offset_y = -np.sin(offset_theta) * offset
    ax_cont.arrow(xfl + offset_x, yfl + offset_y, mult * (xfr - xfl), mult * (yfr - yfl),
                  width=1, length_includes_head=True, linewidth=0, color='C1')

ax_cont.set_title('Continuation: '
                  r'$J$ = $-\phi_f \cos{\xi} - 3 \theta_f \sin{\xi}$', fontsize=10)
ax_cont.set_aspect('equal')
ax_cont.set_xlabel(r'$\phi$ [deg]')
ax_cont.set_ylabel(r'$\theta$ [deg]')

plt.subplots_adjust(
        top=0.84,
        bottom=0.348,
        left=0.193,
        right=0.926,
        hspace=0.2,
        wspace=0.2
)

plt.show()
