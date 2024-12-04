import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.integrate import ode
import os 

import matplotlib.dates as mdates
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib_defaults import matplotlib_defaults, figsizes, mints_colors
from havok import build_hankel, sHAVOK, make_expM, eval_havok

matplotlib_defaults(font='palatino')


# set up plotting paths
figpath = "./figures/1__havok-lorenz"
if not os.path.exists(figpath):
    os.makedirs(figpath)



# 1: Generate Lorenz data

du = np.zeros(3)
def lorenz(t, u, du, params):
    sigma, rho, beta = params
    x, y, z = u

    du[0] = sigma * (y - x)
    du[1] = (x * (rho - z)) - y
    du[2] = (x * y) - (beta * z)

    return du

params = [10.0, 28.0, 8/3]
u0 = [-8.0, 8.0, 27.0]
dt = 0.001
ts = np.arange(0, 200 + dt, dt)
sol = np.empty((len(u0), len(ts)))

sol[:,0] = u0

mdl = ode(lorenz).set_integrator('dopri5')
mdl.set_initial_value(u0, ts[0])
mdl.set_f_params(du, params)
for i, t in enumerate(ts):
    if i == 0:
        continue
    mdl.integrate(t)
    sol[:,i] = mdl.y


plt.clf()
fig, axs = plt.subplots(3, 1, figsize=figsizes['wide'], sharex=True)

axs[0].plot(ts, sol[0,:], color=mints_colors[0], linewidth=1)
axs[1].plot(ts, sol[1,:], color=mints_colors[1], linewidth=1)
axs[2].plot(ts, sol[2,:], color=mints_colors[2], linewidth=1)

axs[2].set_xlim(ts[0], ts[-1])

axs[0].set_ylabel('x')
axs[1].set_ylabel('y')
axs[2].set_ylabel('z')
axs[2].set_xlabel('t')

axs[0].tick_params(axis='x', which='both', length=0)
axs[1].tick_params(axis='x', which='both', length=0)


# turn on major/minor girds
for ax in axs:
    ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
    ax.grid(which='minor', color='#DDDDDD', linewidth=0.5)
    ax.minorticks_on()

axs[0].set_xscale('linear')
axs[1].set_xscale('linear')
axs[2].set_xscale('linear')

plt.tight_layout()
plt.subplots_adjust(hspace=0.25)

plt.savefig(os.path.join(figpath, "1__lorenz-xyz-timeseries.png"))
plt.close()


# 2: Visualize the attractor
L = ts < 50.0

plt.clf()
fig, ax = plt.subplots(figsize=figsizes['default'], subplot_kw={'projection':'3d'})
ax.plot3D(sol[0,L], sol[1,L], sol[2,L], linewidth=1, color='k', alpha=0.65)
ax.set_axis_off()
ax.view_init(elev=30, azim=-35)

plt.savefig(os.path.join(figpath, "1b_lorenz-attractor.png"))
plt.close()



# 2: fit HAVOK model
n_embedding = 201
r_model = 14
n_control = 1
r = r_model + n_control

# pick out the time series and the times
z = sol[0,:]
ts


# compute HAVOK decomposition
z_x, z_pred, t_x, U, s, V, A, B, fvals = eval_havok(z, ts, n_embedding, r_model, n_control)


plt.clf()
plt.figure()
plt.plot(t_x, z_x)
plt.plot(t_x, z_pred)
plt.xlim(t_x[0], t_x[10000])
plt.show()



# construct Hankel matrix
# visualize the extracted matrices

vmin = -10
vmax = 10
cmap = LinearSegmentedColormap.from_list('list', [mints_colors[2], '#FFFFFF', mints_colors[1]])


plt.clf()
plt.close('all')

fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [r_model, n_control]}, figsize=figsizes['default'])
plt.xscale('linear')

axs[0].set_xticks([])
axs[1].set_xticks([])
axs[0].set_yticks([])
axs[1].set_yticks([])

im_A = axs[0].imshow(A, cmap=cmap, vmin=vmin, vmax=vmax)
im_B = axs[1].imshow(B, cmap=cmap, vmin=vmin, vmax=vmax)

axs[0].set_title('A')
axs[1].set_title('B')

fig.subplots_adjust(right=0.8, wspace=0.05)
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
cbar_ax.set_xticks(range(vmin,vmax, 5))
cb = fig.colorbar(im_A, cax=cbar_ax, extend='both')

plt.savefig(os.path.join(figpath, "2__A-B-operator-heatmap.png"))
plt.close()



