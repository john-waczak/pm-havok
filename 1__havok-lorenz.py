import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.integrate import ode
from scipy import stats
import os 

import matplotlib.dates as mdates
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib_defaults import matplotlib_defaults, figsizes, mints_colors, colored_line
from havok import build_hankel, sHAVOK, make_expM, eval_havok, eval_havok_multi, integrate_havok

matplotlib_defaults(font='palatino')


# set up plotting paths
outpath = "./output/1__havok-lorenz"
if not os.path.exists(outpath):
    os.makedirs(outpath)



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
t = np.arange(0, 200 + dt, dt)
sol = np.empty((len(u0), len(t)))

sol[:,0] = u0

mdl = ode(lorenz).set_integrator('dopri5')
mdl.set_initial_value(u0, t[0])
mdl.set_f_params(du, params)
for i, t_i in enumerate(t):
    if i == 0:
        continue
    mdl.integrate(t_i)
    sol[:,i] = mdl.y


plt.clf()
fig, axs = plt.subplots(3, 1, figsize=figsizes['wide'], sharex=True)

axs[0].plot(t, sol[0,:], color=mints_colors[0], linewidth=1)
axs[1].plot(t, sol[1,:], color=mints_colors[1], linewidth=1)
axs[2].plot(t, sol[2,:], color=mints_colors[2], linewidth=1)

axs[2].set_xlim(t[0], t[-1])

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

plt.savefig(os.path.join(outpath, "1__lorenz-xyz-timeseries.png"))
plt.close()


# 2: Visualize the attractor
L = t < 50.0

plt.clf()
fig, ax = plt.subplots(figsize=figsizes['default'], subplot_kw={'projection':'3d'})
ax.plot3D(sol[0,L], sol[1,L], sol[2,L], linewidth=1, color='k', alpha=0.65)
ax.set_axis_off()
ax.view_init(elev=30, azim=-35)

plt.savefig(os.path.join(outpath, "1b_lorenz-attractor.png"))
plt.close()



# 2: fit HAVOK model
n_embedding = 201
r_model = 14
n_control = 1
r = r_model + n_control

# pick out the time series and the times
z = sol[0,:]


# compute HAVOK decomposition
z_x, z_pred, t_x, U, s, V, A, B, fvals = eval_havok(z, t, n_embedding, r_model, n_control)



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

# plt.show()

plt.savefig(os.path.join(outpath, "2__A-B-operator-heatmap.png"))
plt.close()



# 3 plot the reconstruction
fig, ax = plt.subplots(figsize=figsizes['wide'])

l_orig, = ax.plot(t_x, z_x, label="True")
l_havok, = ax.plot(t_x, z_pred, label="HAVOK")
ax.set_xlim(t_x[0], 50)
ax.set_xscale('linear')

ax.set_xlabel('t')
ax.set_ylabel('x(t)')
ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
ax.grid(which='minor', color='#DDDDDD', linewidth=0.5)
ax.minorticks_on()

# add legend on top of plot
plt.legend(handles=[l_orig, l_havok], bbox_to_anchor=(0, 1.015, 1, 0.15), loc="upper center", frameon=False, borderaxespad=0, ncol=2)

plt.tight_layout()

plt.savefig(os.path.join(outpath, "3__havok-reconstruction.png"))
plt.clf()


# Plot the eigenmodes from the U matrix
fig, ax = plt.subplots(figsize=figsizes['default'])

ls1 = []
ls2 = []
lr = []
for i in range(r):
    if i < 3 :
        l, = ax.plot(np.arange(n_embedding) * dt, U[:,i], color=mints_colors[2], linewidth=1)
        ls1.append(l)
    elif i >= 3 and i < r-1:
        l, = ax.plot(np.arange(n_embedding) * dt, U[:,i], color='tab:gray', alpha=0.5, linewidth=0.5) # gray
        ls2.append(l)
    else:
        l, = ax.plot(np.arange(n_embedding) * dt, U[:,i], color=mints_colors[1], linewidth=1)
        lr.append(l)
        
        
plt.legend([*ls1, ls2[0], lr[0]], [r"$u_1$", r"$u_2$", r"$u_3$", r"$\vdots$", r"$u_r$"], bbox_to_anchor=(1.02, 0.4), loc="center left", frameon=False, borderaxespad=0, fontsize=10, handletextpad=0.2, handlelength=0.8)

ax.set_xlim(0 - 0.1*dt, n_embedding*dt + 0.1*dt)
ax.set_xlabel('t')
ax.set_title('Eigenmodes')
ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
ax.grid(which='minor', color='#DDDDDD', linewidth=0.5)
ax.minorticks_on()

plt.tight_layout()

plt.savefig(os.path.join(outpath, "4__U-eigenmodes.png"))
plt.clf()



# plot with forcing
thresh = 5e-6
forcing_active = np.zeros(fvals.shape[0]) #, dtype=bool)
idx_high = np.where(fvals[:,0] ** 2 >= thresh)[0]
win = 500 # max window size to consider forcing active
forcing_active[idx_high] = 1.0

sum(forcing_active)

for i in range(1, len(idx_high)):
    if idx_high[i] - idx_high[i-1] <= win:
        forcing_active[idx_high[i-1]:idx_high[i]] = 1.0


sum(forcing_active)

on_off = LinearSegmentedColormap.from_list('list', [mints_colors[0], mints_colors[1]], N=2)

lcolor = [mints_colors[1] if fa == 1 else mints_colors[0] for fa in forcing_active ]


plt.close()
fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=figsizes['wide'], sharex=True)

ls = colored_line(t_x, z_x, forcing_active, ax=axs[0], cmap=on_off, linewidth=0.75, alpha=1)
ls = colored_line(t_x, fvals[:,0] ** 2, forcing_active, axs[1], cmap=on_off, linewidth=0.75, alpha=1)

axs[0].set_ylim(-20, 20)
axs[1].set_ylim(-0.00001, 0.0006)
axs[1].set_xlim(t_x[0], 50)

axs[0].grid(which='major', color='#DDDDDD', linewidth=0.8)
axs[0].grid(which='minor', color='#DDDDDD', linewidth=0.5)
axs[0].minorticks_on()

axs[1].grid(which='major', color='#DDDDDD', linewidth=0.8)
axs[1].grid(which='minor', color='#DDDDDD', linewidth=0.5)
axs[1].minorticks_on()

axs[0].tick_params(axis='x', which='both', length=0)

axs[0].set_ylabel(r'$x(t)$')
axs[1].set_ylabel(r'$v_r^2$')
axs[1].set_xlabel('time')

plt.tight_layout()
plt.subplots_adjust(hspace=0.1)

plt.savefig(os.path.join(outpath, '5__timeseries-w-forcing.png'))
plt.clf()



# plot the KDE
kde = stats.gaussian_kde(fvals[:,0])
vs = np.linspace(-0.1, 0.1, 1000)

fig, ax = plt.subplots(figsize=figsizes['default'])
l1, = ax.plot(vs, kde(vs), color=mints_colors[2], label="Estimated PDF")
l2, = ax.plot(vs, stats.norm.pdf(vs, np.mean(fvals[:,0]), np.std(fvals[:,0])), color=mints_colors[1], linestyle='--', label="Gasusian")
ax.set_yscale('log')
ax.set_ylim(1e-1, 1e3)
ax.set_xlim(-0.02, 0.02)
ax.set_xlabel(r'$v_r$')
ax.set_title('Forcing Statistics')

plt.legend(handles=[l1, l2], bbox_to_anchor=(1.04, 0.4), loc="center left", frameon=False, borderaxespad=0, fontsize=10, handletextpad=0.2, handlelength=1.5)

ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
ax.grid(which='minor', color='#DDDDDD', linewidth=0.5)
ax.minorticks_on()

plt.tight_layout()

plt.savefig(os.path.join(outpath, '6__forcing-statistics.png'))
plt.clf()



# test the integrate funtion
# to make sure it works as expected
# z_test, z_pred, t_test = integrate_havok(z[1000:5000], t[1000:5000], n_embedding, r_model, n_control, A, B, U, s)

# plt.figure()
# plt.plot(t_test, z_test)
# plt.plot(t_test, z_pred)
# plt.show()





# Now let's do it for multiple time series

# split the time series into chunks:
N = 5000
skip = 100
zs = [z[i:i+N] for i in range(0, len(z) - N + 1, N + skip)]
ts = [t[i:i+N] for i in range(0, len(t) - N + 1, N + skip)]

# visualize the disjoint time series
fig, ax = plt.subplots(figsize=figsizes['wide'])

ax.set_xscale('linear')
for i in range(len(zs)):
    ax.plot(ts[i], zs[i], linewidth=1, color=mints_colors[2])

ax.set_xlabel("t")    
ax.set_ylabel("x(t)")

ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
ax.grid(which='minor', color='#DDDDDD', linewidth=0.5)
ax.minorticks_on()

ax.set_xlim(ts[0][0], 50)
plt.tight_layout()
ax.set_xscale('linear')

plt.savefig(os.path.join(outpath, "7__timeseries-disjoint.png"))
plt.clf()


# Now let's do the HAVOK bit
zs_x, zs_pred, ts_x, U, s, Vs_out, A, B, fvals = eval_havok_multi(zs, ts, n_embedding, r_model, n_control)


set([len(z) for z in zs_x])

set([len(t) for t in ts_x])



# visualize the disjoint time series
fig, ax = plt.subplots(figsize=figsizes['wide'])

ax.set_xscale('linear')
ls = []
for i in range(len(zs_x)):
    l1, = ax.plot(ts_x[i], zs_x[i], linewidth=1, color=mints_colors[2])
    l2, = ax.plot(ts_x[i], zs_pred[i], linewidth=1, color=mints_colors[1])

    if i == 0 :
        ls.append(l1)
        ls.append(l2)

# add legend on top of plot
plt.legend(ls, ["Original", "HAVOK"], bbox_to_anchor=(0, 1.015, 1, 0.15), loc="upper center", frameon=False, borderaxespad=0, ncol=2)


ax.set_xlabel("t")    
ax.set_ylabel("x(t)")

ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
ax.grid(which='minor', color='#DDDDDD', linewidth=0.5)
ax.minorticks_on()

ax.set_xlim(ts[0][0], 50)
plt.tight_layout()
ax.set_xscale('linear')

plt.savefig(os.path.join(outpath, "8__havok-reconstruction-disjoint.png"))
plt.clf()




