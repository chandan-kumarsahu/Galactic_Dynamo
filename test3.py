import matplotlib.pyplot as plt
import numpy as np


def source_term(x, t):
    """
    Source term function.

    Parameters:
    - x: Spatial coordinate
    - t: Time coordinate

    Returns:
    - Source term value at (x, t)
    """
    # You can modify this function to suit your specific source term requirements
    return 10 * np.sin(np.pi * x) * np.exp(-0.1 * t)

def get_diff_matrix(N, sigma):
    diag_values = 2 + 2 * sigma
    off_diag_values = -sigma

    A = np.diag(diag_values * np.ones(N)) + np.diag(off_diag_values * np.ones(N - 1), 1) + np.diag(off_diag_values * np.ones(N - 1), -1)
    B = np.diag((2 - 2 * sigma) * np.ones(N)) + np.diag(sigma * np.ones(N - 1), 1) + np.diag(sigma * np.ones(N - 1), -1)

    # Boundary conditions
    A[0, 0] = 1 + sigma
    B[0, 0] = 1 - sigma
    A[-1, -1] = 1 + sigma
    B[-1, -1] = 1 - sigma

    return A, B

def crank_nicolson_diffusion_with_source(x_max, t_max, dx, dt, Diff, init_cond, source_term):

    alpha = Diff * dt / (dx**2)

    # Spatial grid
    x = [i*dx for i in range(int(x_max/dx)+1)]
    t = [j*dt for j in range(int(t_max/dt)+1)]

    # Initialize temperature array
    Temp = np.zeros((len(x), len(t)))

    # Initial condition
    for i in range(len(x)):
        Temp[i][0] = init_cond(x[i])

    # Get the matrices for solving the matrix using crank-nicolson method
    A, B = get_diff_matrix(len(x), alpha)

    A = np.array(A)
    B = np.array(B)

    for j in range(1, len(t)):
        source_vector = np.array([source_term(xi, t[j]) for xi in x])
        Temp[:, j] = np.linalg.solve(A, np.dot(B, Temp[:, j - 1]) + dt * source_vector)

    return Temp, np.array(x), np.array(t)


def init_cond(x):
    return -500*(x-0.5)*np.exp(-20*(x-0.5)**2)

# Constants and parameters
eta_T = 1e-4    # thermal diffusivity
t_max = 500     # total simulation time
z_max = 1.0     # thickness of the disc
dt = 0.05       # time step
dz = 0.01       # spatial step in z direction

solution_with_source, spatial_grid, time_grid = crank_nicolson_diffusion_with_source(z_max, t_max, dz, dt, eta_T, init_cond, source_term)

plt.figure(figsize=(12, 6))
# Plot the diffusion equation solution with source term
plt.subplot(1, 2, 1)
plt.imshow(solution_with_source, aspect='auto', extent=[0, t_max, 0, z_max], origin='lower', cmap='Spectral_r')
plt.colorbar(label='Temperature')
plt.xlabel('Time')
plt.ylabel('Spatial coordinate')
plt.title('Diffusion Equation with Source Term')

# plot without source term
solution_without_source, spatial_grid, time_grid = crank_nicolson_diffusion_with_source(z_max, t_max, dz, dt, eta_T, init_cond, lambda x, t: 0)

plt.subplot(1, 2, 2)
plt.imshow(solution_without_source, aspect='auto', extent=[0, t_max, 0, z_max], origin='lower', cmap='Spectral_r')
plt.colorbar(label='Temperature')
plt.xlabel('Time')
plt.ylabel('Spatial coordinate')
plt.title('Diffusion Equation without Source Term')


plt.show()
