import matplotlib.pyplot as plt
import giuseppe as gp
from matplotlib import rcParams, gridspec, rc

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

cmap = plt.colormaps['viridis']

sol_set_alg = gp.data_classes.load_sol_set('sol_set_alg.bin')
sol_set_diff = gp.data_classes.load_sol_set('sol_set_diff.bin')

fig_cont = plt.figure(figsize=(5, 3.5))
ax_cont = fig_cont.add_subplot()

xf, yf = [], []
for sol in sol_set_alg:
    ax_cont.plot(sol.x[0, :], sol.x[1, :], color='C0')
    xf.append(sol.x[0, -1])
    yf.append(sol.x[1, -1])
ax_cont.plot(xf, yf, '--', marker='*', label='', color='C1')

offset, mult = 0.25, 0.95
for xfl, xfr, yfl, yfr in zip(xf[:-1], xf[1:], yf[:-1], yf[1:]):
    ax_cont.arrow(xfl + offset, yfl + offset, mult * (xfr - xfl), mult * (yfr - yfl),
                  width=0.1, length_includes_head=True, linewidth=0, color='C1')

ax_cont.set_title('Continuation Example:\nBrachistochrone Problem')
ax_cont.set_aspect('equal', 'box')
ax_cont.set_xlabel(r'$x$ [ft]')
ax_cont.set_ylabel(r'$y$ [ft]')


sol_alg = sol_set_alg[-1]
sol_diff = sol_set_diff[-1]

fig_ex = plt.figure(figsize=(6.5, 4))
fig_ex.suptitle('Brachistochrone Problem')

gs = gridspec.GridSpec(2, 3)

ax_path = fig_ex.add_subplot(gs[0, 0])
ax_path.plot(sol_alg.x[0, :], sol_alg.x[1, :], '-', label='Algebraic Control Law')
ax_path.plot(sol_diff.x[0, :], sol_diff.x[1, :], '--', label='Differential Control Law')
ax_path.set_xlabel(r'$x$ [ft]')
ax_path.set_ylabel(r'$y$ [ft]')
fig_ex.legend(loc='lower center', ncols=2)

ax_vel = fig_ex.add_subplot(gs[0, 1])
ax_vel.plot(sol_alg.t, sol_alg.x[2, :], '-', label='Algebraic Control Law')
ax_vel.plot(sol_diff.t, sol_diff.x[2, :], '--', label='Differential Control Law')
ax_vel.set_xlabel(r'$t$ [s]')
ax_vel.set_ylabel(r'$v$ [ft/s]')

ax_theta = fig_ex.add_subplot(gs[0, 2])
ax_theta.plot(sol_alg.t, sol_alg.u[0, :] * 180 / 3.141596, '-', label='Algebraic Control Law')
ax_theta.plot(sol_diff.t, sol_diff.u[0, :] * 180 / 3.141596, '--', label='Differential Control Law')
ax_theta.set_xlabel(r'$t$ [s]')
ax_theta.set_ylabel(r'$\theta$ [deg]')

ax_lam_x = fig_ex.add_subplot(gs[1, 0])
ax_lam_x.plot(sol_alg.t, sol_alg.lam[0, :], '-', label='Algebraic Control Law')
ax_lam_x.plot(sol_diff.t, sol_diff.lam[0, :], '--', label='Differential Control Law')
ax_lam_x.set_xlabel(r'$t$ [s]')
ax_lam_x.set_ylabel(r'$\lambda_x$ [1/ft]')
ax_lam_x.set_ylim(-0.1, 0)

ax_lam_y = fig_ex.add_subplot(gs[1, 1])
ax_lam_y.plot(sol_alg.t, sol_alg.lam[1, :], '-', label='Algebraic Control Law')
ax_lam_y.plot(sol_diff.t, sol_diff.lam[1, :], '--', label='Differential Control Law')
ax_lam_y.set_xlabel(r'$t$ [s]')
ax_lam_y.set_ylabel(r'$\lambda_y$ [1/ft]')
ax_lam_y.set_ylim(0, 0.1)

ax_lam_v = fig_ex.add_subplot(gs[1, 2])
ax_lam_v.plot(sol_alg.t, sol_alg.lam[2, :], '-', label='Algebraic Control Law')
ax_lam_v.plot(sol_diff.t, sol_diff.lam[2, :], '--', label='Differential Control Law')
ax_lam_v.set_xlabel(r'$t$ [s]')
ax_lam_v.set_ylabel(r'$\lambda_v$ [s/ft]')

fig_ex.subplots_adjust(
        top=0.9,
        bottom=0.2,
        left=0.15,
        right=0.95,
        hspace=0.5,
        wspace=0.75
)

fig_cont.tight_layout()

plt.show()
