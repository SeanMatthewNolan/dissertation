import matplotlib.pyplot as plt
import giuseppe as gp
import numpy as np
from matplotlib import rc, rcParams, gridspec

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'


cmap = plt.colormaps['viridis']

sol_set_path = gp.data_classes.load_sol_set('sol_set_nav.bin')
sol_set_term = gp.data_classes.load_sol_set('sol_set_nav_final.bin')

sol_path = sol_set_path[-1]
sol_term = sol_set_term[-1]
sol_none = gp.data_classes.load_sol('nav_straight.bin')

xb_path, yb_path = sol_path.k[1], sol_path.k[2]
xb_term, yb_term = sol_term.k[1], sol_term.k[2]
xb_none, yb_none = sol_none.k[1], sol_none.k[2]

fig_ex = plt.figure(figsize=(6.5, 7.5))
fig_ex.suptitle('Navigation Cart Problem')

gs = gridspec.GridSpec(4, 2)

ax_path = fig_ex.add_subplot(gs[0, 0])
ax_path.plot(sol_path.x[0, :], sol_path.x[1, :], label='Path Cost')
ax_path.plot(sol_term.x[0, :], sol_term.x[1, :], label='Terminal Cost')
ax_path.plot(sol_none.x[0, :], sol_none.x[1, :], '--', label='Min Time Straight')
ax_path.plot(xb_path, yb_path, '*', label='Beacon')
ax_path.set_xlabel(r'Position, $x$ [ft]')
ax_path.set_ylabel(r'Position, $y$ [ft]')

fig_ex.legend(loc='lower center', ncols=4)

ax_vel = fig_ex.add_subplot(gs[0, 1])
ax_vel.plot(sol_path.t, np.rad2deg(sol_path.x[2, :]), label='Path Cost')
ax_vel.plot(sol_term.t, np.rad2deg(sol_term.x[2, :]), label='Terminal Cost')
ax_vel.plot(sol_none.t, np.rad2deg(sol_none.x[2, :]), '--', label='Min Time Straight')
ax_vel.set_xlabel(r'$t$ [s]')
ax_vel.set_ylabel(r'Heading Angle,' + '\n' + r'$\theta$ [deg]')

ax_theta = fig_ex.add_subplot(gs[1, 0])
ax_theta.plot(sol_path.t, np.rad2deg(sol_path.u[0, :]), label='Path Cost')
ax_theta.plot(sol_term.t, np.rad2deg(sol_term.u[0, :]), label='Terminal Cost')
ax_theta.plot(sol_none.t, np.rad2deg(sol_none.u[0, :]), '--', label='Min Time Straight')
ax_theta.set_xlabel(r'$t$ [s]')
ax_theta.set_ylabel(r'Turn Rate,' + '\n' + r'$\omega$ [deg / s]')

ax_rho = fig_ex.add_subplot(gs[1, 1])
ax_rho.plot(sol_path.t, np.sqrt((sol_path.x[0, :] - xb_path) ** 2 + (sol_path.x[1, :] - yb_path) ** 2),
            label='Path Cost')
ax_rho.plot(sol_term.t, np.sqrt((sol_term.x[0, :] - xb_term) ** 2 + (sol_term.x[1, :] - yb_term) ** 2),
            label='Terminal Cost')
ax_rho.plot(sol_none.t, np.sqrt((sol_none.x[0, :] - xb_none) ** 2 + (sol_none.x[1, :] - yb_none) ** 2),
            '--', label='Min Time Straight')
ax_rho.set_xlabel(r'$t$ [s]')
ax_rho.set_ylabel(r'Range Measurement,' + '\n' + r'$\rho$ [ft]')

ax_var_x = fig_ex.add_subplot(gs[2, 0])
ax_var_x.plot(sol_path.t, sol_path.x[3, :], label='Path Cost')
ax_var_x.plot(sol_term.t, sol_term.x[3, :], label='Terminal Cost')
ax_var_x.plot(sol_none.t, sol_none.x[3, :], '--', label='Min Time Straight')
ax_var_x.set_xlabel(r'$t$ [s]')
ax_var_x.set_ylabel(r'Variance in $x$,' + '\n' + r'$p_{xx}$ [$\text{ft}^2$]')

ax_var_y = fig_ex.add_subplot(gs[2, 1])
ax_var_y.plot(sol_path.t, sol_path.x[6, :], label='Path Cost')
ax_var_y.plot(sol_term.t, sol_term.x[6, :], label='Terminal Cost')
ax_var_y.plot(sol_none.t, sol_none.x[6, :], '--', label='Min Time Straight')
ax_var_y.set_xlabel(r'$t$ [s]')
ax_var_y.set_ylabel(r'Variance in $y$,' + '\n' + r'$p_{yy}$ [$\text{ft}^2$]')

ax_var_x = fig_ex.add_subplot(gs[3, 0])
ax_var_x.plot(sol_path.t, sol_path.x[8, :] * 180**2 / np.pi**2, label='Path Cost')
ax_var_x.plot(sol_term.t, sol_term.x[8, :] * 180**2 / np.pi**2, label='Terminal Cost')
ax_var_x.plot(sol_none.t, sol_none.x[8, :] * 180**2 / np.pi**2, '--', label='Min Time Straight')
ax_var_x.set_xlabel(r'$t$ [s]')
ax_var_x.set_ylabel(r'Variance in $\theta$,' + '\n' + r'$p_{\theta\theta}$ [$\text{deg}^2$]')

ax_var_y = fig_ex.add_subplot(gs[3, 1])
ax_var_y.plot(sol_path.t, sol_path.x[3, :] + sol_path.x[6, :], label='Path Cost')
ax_var_y.plot(sol_term.t, sol_term.x[3, :] + sol_term.x[6, :], label='Terminal Cost')
ax_var_y.plot(sol_none.t, sol_none.x[3, :] + sol_none.x[6, :], '--', label='Min Time Straight')
ax_var_y.set_xlabel(r'$t$ [s]')
ax_var_y.set_ylabel('Position Variance,\n' + r'$p_{xx} + p_{yy}$ [$\text{ft}^2$]')

fig_ex.tight_layout()

fig_ex.subplots_adjust(
        top=0.94,
        bottom=0.11,
        left=0.117,
        right=0.977,
        hspace=0.471,
        wspace=0.365
)

plt.show()
