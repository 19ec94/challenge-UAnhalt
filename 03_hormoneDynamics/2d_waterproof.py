import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math

# Domain parameters
Lx = 1.0
Ly = 1.0

Nx = 400  # number of grid cells in x
Ny = 400  # number of grid cells in y

dx = Lx / Nx
dy = Ly / Ny

# Interface parameters
ys = 0.3  # y-coordinate of horizontal interface
xs = 0.5  # x-coordinate splitting Ω2 horizontally

# Velocity in zones for waterproof case:
# Ω1 and Ω3: diagonal flow (1,1)
# Ω2: horizontal flow (1,0)
def velocity(x, y):
    vel_x = np.ones_like(x)
    vel_y = np.ones_like(y)
    # Identify Ω2 zone: x>xs and y<ys
    omega2 = (x >= xs) & (y < ys)
    vel_y[omega2] = 0.0  # vertical velocity zero in Ω2 (waterproof)
    return vel_x, vel_y

# Initial condition: Gaussian bump centered at (cx, cy)
sigma = math.sqrt(0.002) 
cx = 0.3
cy = 0.15

def initial_condition(x, y):
    scale = 1/(2 * np.pi * sigma**2)
    return np.exp(-((x - cx)**2 + (y - cy)**2) / (2 * sigma**2))

# Initialize grid
x = (np.arange(Nx) + 0.5)*dx
y = (np.arange(Ny) + 0.5)*dy
X, Y = np.meshgrid(x, y, indexing='ij')

phi = initial_condition(X, Y)

# CFL condition parameters
vel_x, vel_y = velocity(X, Y)
max_vel_x = np.max(np.abs(vel_x))
max_vel_y = np.max(np.abs(vel_y))
cfl = 0.4
dt = cfl * min(dx/max(max_vel_x,1e-8), dy/max(max_vel_y,1e-8))

# Upwind flux calculation in 2D with waterproof boundary enforcement
def upwind_flux(phi, vel_x, vel_y, dx, dy, dt):
    Nx, Ny = phi.shape
    flux_phi = np.zeros_like(phi)

    # Fluxes in x direction with periodic BC
    phi_x_plus = np.roll(phi, -1, axis=0)
    phi_x_minus = np.roll(phi, 1, axis=0)

    # Fluxes in y direction with zero Dirichlet BC at top boundary y=Ly (waterproof effect)
    phi_y_plus = np.zeros_like(phi)
    phi_y_plus[:, :-1] = phi[:, 1:]
    # For last row, no cells above, Dirichlet zero implicitly

    phi_y_minus = np.roll(phi, 1, axis=1)

    # x-direction flux
    flux_x = np.zeros_like(phi)
    pos_vx = vel_x >= 0
    neg_vx = ~pos_vx

    flux_x[pos_vx] = vel_x[pos_vx] * phi[pos_vx]
    flux_x[neg_vx] = vel_x[neg_vx] * phi_x_plus[neg_vx]

    # y-direction flux
    flux_y = np.zeros_like(phi)
    pos_vy = vel_y >= 0
    neg_vy = ~pos_vy

    flux_y[pos_vy] = vel_y[pos_vy] * phi[pos_vy]
    flux_y[neg_vy] = vel_y[neg_vy] * phi_y_plus[neg_vy]

    # Compute derivatives (flux differences)
    dflux_x = (flux_x - np.roll(flux_x, 1, axis=0)) / dx
    dflux_y = (flux_y - np.roll(flux_y, 1, axis=1)) / dy

    # Update phi according to transport equation
    phi_new = phi - dt * (dflux_x + dflux_y)

    # Enforce impermeable boundary at y=ys for zone Ω2 (x>=xs, y<ys)
    # This means no vertical diffusion of density through y=ys line at Ω2 zone.
    # We zero out the density above ys when x >= xs to mimic Dirichlet zero at upper Ω2 border.
    # It mimics phi being zero just above y=ys and no vertical flux crossing.

    # mask_upper_omega2 = (X >= xs) & (Y >= ys)
    # phi_new[mask_upper_omega2] = 0.0

    return phi_new

# Time integration loop
t = 0.0
t_final = 0.25
nsteps = int(t_final/dt) + 1

# For visualization
fig, ax = plt.subplots(figsize=(6,5))
im = ax.imshow(phi.T, extent=[0,Lx,0,Ly], origin='lower',
               cmap='bwr', vmin=0, vmax=1)
hline = ax.axhline(y=ys, color='black', linewidth=2, label=f'Interface y={ys}')
#vline = ax.axvline(x=xs, ymin=0.0, ymax=0.3, ymax=ys/Ly,  color='black', linewidth=2)
vline = ax.axvline(x=xs, ymin=0.0, ymax=0.3,  color='black', linewidth=2)
title_text = ax.set_title('Waterproof transport, t=0')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend()

def update(frame):
    global phi, t
    phi = upwind_flux(phi, *velocity(X, Y), dx, dy, dt)
    t += dt
    im.set_data(phi.T)
    title_text.set_text(f'Waterproof transport, t={t:.3f}')
    return [im, hline, vline, title_text]

ani = FuncAnimation(fig, update, frames=nsteps, interval=50, blit=True, repeat=False)

plt.show()

