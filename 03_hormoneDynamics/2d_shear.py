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

# Shear interface parameters
ys = 0.3  # y-coordinate of horizontal interface where velocity changes

# Velocity in bottom region: vertical flow (0,1)
# Velocity in top region: oblique flow (1,1)
def velocity(x, y):
    # y is a 2D array
    vel_x = np.zeros_like(y)
    vel_y = np.ones_like(y)
    vel_x[y >= ys] = 1.0
    vel_y[y >= ys] = 1.0
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

# Upwind flux calculation in 2D
def upwind_flux(phi, vel_x, vel_y, dx, dy, dt):
    Nx, Ny = phi.shape
    flux_phi = np.zeros_like(phi)

    # Fluxes in x direction
    # Periodic BC in x (or zero flux at boundaries)
    phi_x_plus = np.roll(phi, -1, axis=0)
    phi_x_minus = np.roll(phi, 1, axis=0)

    # Fluxes in y direction
    phi_y_plus = np.roll(phi, -1, axis=1)
    phi_y_minus = np.roll(phi, 1, axis=1)

    # Compute fluxes using upwind scheme:
    # Flux_x[i,j] = vel_x[i,j]*phi at the upwind location

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

    # Boundary conditions:
    # In x direction: zero inflow BC -> zero flux inflow
    # Set ghost cells to zero inflow
    # Here we use periodic BC for simplicity

    return phi_new

# Time integration loop
t = 0.0
t_final = 0.25
nsteps = int(t_final/dt) + 1

# For visualization
fig, ax = plt.subplots(figsize=(6,5))
im = ax.imshow(phi.T, extent=[0,Lx,0,Ly], origin='lower',
               cmap='bwr', vmin=0, vmax=1)
#ax.axhline(ys, color='red', linestyle='--', label='Interface y=%.2f'%ys)
hline = ax.axhline(y=0.3, color='black', linewidth=2, label='Interface y=%.2f'%ys)
vline = ax.axvline(x=0.5, ymin=0.0, ymax=0.3,  color='black', linewidth=2)
ax.set_title('Shear transport, t=0')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend()

def update(frame):
    global phi, t
    phi = upwind_flux(phi, *velocity(X, Y), dx, dy, dt)
    t += dt
    im.set_data(phi.T)
    ax.set_title(f'Shear transport, t={t:.3f}')
    return [im, hline, vline]

ani = FuncAnimation(fig, update, frames=nsteps, interval=50, blit=True, repeat=False)

plt.show()

