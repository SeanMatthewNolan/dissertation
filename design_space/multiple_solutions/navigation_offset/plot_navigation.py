import matplotlib.pyplot as plt
import giuseppe as gp
import numpy as np
from matplotlib import rc, rcParams, gridspec

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'


sol_set_path = gp.data_classes.load_sol_set('sol_set_nav_from_left.bin')
sol_set_term = gp.data_classes.load_sol_set('sol_set_nav_from_right.bin')

sol_path = sol_set_path[-1]
sol_term = sol_set_term[-1]

xb_path, yb_path = sol_path.k[1], sol_path.k[2]
xb_term, yb_term = sol_term.k[1], sol_term.k[2]

fig_ex = plt.figure(figsize=(6.5, 4))
fig_ex.suptitle('Navigation Cart Problem:\nMultiple Solutions')

gs = gridspec.GridSpec(2, 2)

ax_path = fig_ex.add_subplot(gs[0, 0])
ax_path.plot(sol_path.x[0, :], sol_path.x[1, :], label='Beacon Started on Left')
ax_path.plot(sol_term.x[0, :], sol_term.x[1, :], label='Beacon Started on Right')
ax_path.plot(xb_path, yb_path, '*', label='Beacon', color='C3')
ax_path.set_xlabel(r'Position, $x$ [ft]')
ax_path.set_ylabel(r'Position, $y$ [ft]')

fig_ex.legend(loc='lower center', ncols=4)

ax_vel = fig_ex.add_subplot(gs[0, 1])
ax_vel.plot(sol_path.t, np.rad2deg(sol_path.x[2, :]), label='Beacon Started on Left')
ax_vel.plot(sol_term.t, np.rad2deg(sol_term.x[2, :]), label='Beacon Started on Right')
ax_vel.set_xlabel(r'$t$ [s]')
ax_vel.set_ylabel(r'Heading Angle,' + '\n' + r'$\theta$ [deg]')

ax_theta = fig_ex.add_subplot(gs[1, 0])
ax_theta.plot(sol_path.t, np.rad2deg(sol_path.u[0, :]), label='Beacon Started on Left')
ax_theta.plot(sol_term.t, np.rad2deg(sol_term.u[0, :]), label='Beacon Started on Right')
ax_theta.set_xlabel(r'$t$ [s]')
ax_theta.set_ylabel(r'Turn Rate,' + '\n' + r'$\omega$ [deg / s]')

ax_var_y = fig_ex.add_subplot(gs[1, 1])
ax_var_y.plot(sol_path.t, sol_path.x[3, :] + sol_path.x[6, :], label='Beacon Started on Left')
ax_var_y.plot(sol_term.t, sol_term.x[3, :] + sol_term.x[6, :], label='Beacon Started on Right')
ax_var_y.set_xlabel(r'$t$ [s]')
ax_var_y.set_ylabel('Position Variance,\n' + r'$p_{xx} + p_{yy}$ [$\text{ft}^2$]')

fig_ex.tight_layout()

fig_ex.subplots_adjust(
        top=0.845,
        bottom=0.2,
        left=0.117,
        right=0.977,
        hspace=0.501,
        wspace=0.365
)

plt.show()
