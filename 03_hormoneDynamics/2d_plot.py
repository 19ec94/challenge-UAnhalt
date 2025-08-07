import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parameters
num_frames = 499# For example, number of CSV files: data_0.csv ... data_9.csv

# Read the first CSV to get x, y coordinates and shape info
initial_data = pd.read_csv('output_step_1.csv')
x_unique = np.sort(initial_data['x'].unique())
y_unique = np.sort(initial_data['y'].unique())
Nx = x_unique.size
Ny = y_unique.size

# Plot setup
fig, ax = plt.subplots(figsize=(6,5))
# Initialize with the first frame data
phi = initial_data['phi'].values.reshape((Nx, Ny))
im = ax.imshow(phi.T, extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
               origin='lower', cmap='bwr', vmin=phi.min(), vmax=phi.max())

title_text = ax.set_title('Waterproof transport, t=0')
ax.set_xlabel('x')
ax.set_ylabel('y')

def update(frame):
    filename = f'output_step_{frame+1}.csv'
    data = pd.read_csv(filename)
    phi = data['phi'].values.reshape((Nx, Ny))
    im.set_data(phi.T)
    title_text.set_text(f'Waterproof transport, t={frame}')
    return [im, title_text]

ani = FuncAnimation(fig, update, frames=num_frames, interval=200, blit=True, repeat=False)

plt.show()

