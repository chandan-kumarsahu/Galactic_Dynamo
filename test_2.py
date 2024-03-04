import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# Constants and parameters
eta_T = 2e-4  # magnetic diffusivity
t_max = 500  # total simulation time
z_max = 1.0  # thickness of the disc
r_max = 1.0  # radius of the disc
dt = 0.05  # time step
dz = 0.01  # spatial step in z direction
dr = 0.01  # spatial step in r direction

# Discretization
Nt = int(t_max / dt) + 1    # number of time steps
Nz = int(z_max / dz) + 1    # number of spatial steps in z direction
Nr = int(r_max / dr) + 1    # number of spatial steps in r direction
t_values = np.linspace(0, t_max, Nt)    # array of time values
z_values = np.linspace(0, z_max, Nz)    # array of z values
r_values = np.linspace(0, r_max, Nz)    # array of r values

# Initialize temperature array with the specified quadratic initial condition
B_z = np.zeros((Nt, Nz))
B_r = np.zeros((Nt, Nz))
B_phi = np.zeros((Nt, Nz))

# Initial conditions for the magnetic field in z, r and phi directions
B_z[0, :] = -300*(z_values-0.5)*np.exp(-2*(z_values-0.5)**2)
B_r[0, :] = 100*np.exp(-5*(z_values)**2)
B_phi[0, :] = -1000*(z_values-0.5)*np.exp(-20*(z_values-0.5)**2)
Omega = 50*np.exp(-5*(r_values)**2)

# Boundary conditions in the z direction
B_z[:, 0] = B_z[:, 1]  # Insulated top end of the disc
B_z[:, -1] = B_z[:, -2]  # Insulated bottom end of the disc

# Boundary conditions in the r direction
B_r[:, 0] = B_r[:, 1]  # Insulated centre of the disc
B_r[:, -1] = B_r[:, -2]  # Insulated at the edge of the disc

# Boundary conditions in the phi direction
B_phi[:, 0] = B_phi[:, 1]  # Insulated left end of the disc
B_phi[:, -1] = B_phi[:, -2]  # Insulated right end of the disc

# Finite difference method in z direction
for n in (range(0, Nt - 1)):
    for i in range(1, Nz - 1):
        B_z[n + 1, i] = B_z[n, i] + eta_T * dt / dz**2 * (B_z[n, i + 1] - 2 * B_z[n, i] + B_z[n, i - 1])

    # Maintain the temperature at the insulated left end
    B_z[n + 1, 0] = B_z[n, 1]
    B_z[n + 1, -1] = B_z[n, -2]

# Finite difference method in r direction
for n in (range(0, Nt - 1)):
    for i in range(1, Nz - 1):
        B_r[n + 1, i] = B_r[n, i] + eta_T * dt / dz**2 * (B_r[n, i + 1] - 2 * B_r[n, i] + B_r[n, i - 1])

    # Maintain the temperature at the insulated left end
    B_r[n + 1, 0] = B_r[n, 1]
    B_r[n + 1, -1] = B_r[n, -2]

# Finite difference method in phi direction
for n in (range(0, Nt - 1)):
    for i in range(1, Nz - 1):
        B_phi[n + 1, i] = B_phi[n, i] + dt / dr * r_values[i]*(Omega[i]-Omega[i-1])*B_r[n, i] + eta_T * dt / dz**2 * (B_phi[n, i + 1] - 2 * B_phi[n, i] + B_phi[n, i - 1])

    # Maintain the temperature at the insulated left end
    B_phi[n + 1, 0] = B_phi[n, 1]
    B_phi[n + 1, -1] = B_phi[n, -2]

# Create a 3D surface plot
Z_values, T_values = np.meshgrid(z_values, t_values)
R_values, T_values = np.meshgrid(z_values, t_values)
Phi_values, T_values = np.meshgrid(z_values, t_values)


# make three 3D subplots for z, r and phi
fig = plt.figure(figsize=(20, 7))
ax1 = fig.add_subplot(131, projection='3d')
ax1.plot_surface(T_values, Z_values, B_z, cmap='viridis')
ax1.set_xlim(0, t_max)
ax1.set_ylim(0, z_max)
ax1.set_zlim(np.min(B_z)-10, np.max(B_z)+10)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Distance (m)')
ax1.set_zlabel('Magnetic Field Strength')
ax1.set_title('Diffusion of Magnetic field in z direction')

ax2 = fig.add_subplot(132, projection='3d')
ax2.plot_surface(T_values, R_values, B_r, cmap='viridis')
ax2.set_xlim(0, t_max)
ax2.set_ylim(0, z_max)
ax2.set_zlim(-np.max(B_r)-10, np.max(B_r)+10)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Distance (m)')
ax2.set_zlabel('Magnetic Field Strength')
ax2.set_title('Diffusion of Magnetic field in r direction')

ax3 = fig.add_subplot(133, projection='3d')
ax3.plot_surface(T_values, Phi_values, B_phi, cmap='viridis')
ax3.set_xlim(0, t_max)
ax3.set_ylim(0, z_max)
ax3.set_zlim(np.min(B_phi)-10, np.max(B_phi)+10)
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Distance (m)')
ax3.set_zlabel('Magnetic Field Strength')
ax3.set_title('Diffusion of Magnetic field in phi direction')

plt.show()
