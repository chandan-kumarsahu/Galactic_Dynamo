import matplotlib.pyplot as plt
import numpy as np

"""
Plot the solution in both 1D and Heatmap format.

Parameters:
- time_grid: Time grid
- spatial_grid: Spatial grid
- solution: Solution of the diffusion equation

Returns:
- Makes the plots
"""

def plot_diff(time_grid, spatial_grid, solution_r, solution_phi):

    # Create 2D plots
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 2, 1)
    for i in (range(0, len(time_grid), int(len(time_grid)/5))):
        plt.plot(spatial_grid, solution_r[:, i], label=f'time = {time_grid[i]:.1f}')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'Magnetic Field Strength ($B_r$)')
    plt.title('Diffusion of Magnetic field in radial direction')
    # plt.ylim(np.min(solution_r), np.max(solution_r))
    plt.grid()
    plt.legend()

    # Create imshow plot
    plt.subplot(2, 2, 2)
    plt.contourf(*np.meshgrid(spatial_grid, time_grid), solution_r.T, 40, cmap='Spectral_r')
    plt.colorbar(label=r'Magnetic Field Strength ($B_r$)')
    plt.title(r'Diffusion of Magnetic field in radial direction')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'Time (Myr)')
    plt.grid()


    # Create 2D plots
    plt.subplot(2, 2, 3)
    for i in (range(0, len(time_grid), int(len(time_grid)/5))):
        plt.plot(spatial_grid, solution_phi[:, i], label=f'time = {time_grid[i]:.1f}')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'Magnetic Field Strength ($B_\phi$)')
    plt.title(r'Diffusion of Magnetic field in azimuthal direction')
    # plt.ylim(np.min(solution_phi), np.max(solution_phi))
    plt.grid()
    plt.legend()

    # Create imshow plot
    plt.subplot(2, 2, 4)
    plt.contourf(*np.meshgrid(spatial_grid, time_grid), solution_phi.T, 40, cmap='Spectral_r')
    plt.colorbar(label=r'Magnetic Field Strength ($B_\phi$)')
    plt.title('Diffusion of Magnetic field in azimuthal direction')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel('Time (Myr)')
    plt.grid()

    plt.tight_layout(pad=3)



def plot_pitch(time_grid, spatial_grid, B, pitch):

    # Create 2D plots
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 2, 1)
    for i in (range(0, len(time_grid), int(len(time_grid)/5))):
        plt.plot(spatial_grid, B[:, i], label=f'time = {time_grid[i]:.1f}')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'Magnetic Field Strength ($B_r$)')
    plt.title('Diffusion of Magnetic field in radial direction')
    # plt.ylim(np.min(solution_r), np.max(solution_r))
    plt.grid()
    plt.legend()

    # Create imshow plot
    plt.subplot(2, 2, 2)
    plt.contourf(*np.meshgrid(spatial_grid, time_grid), B.T, 40, cmap='Spectral_r')
    plt.colorbar(label=r'Magnetic Field Strength ($B_r$)')
    plt.title(r'Diffusion of Magnetic field in radial direction')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'Time (Myr)')
    plt.grid()


    # Create 2D plots
    plt.subplot(2, 2, 3)
    for i in (range(0, len(time_grid), int(len(time_grid)/5))):
        plt.plot(spatial_grid[1:-1], pitch[1:-1, i], label=f'time = {time_grid[i]:.1f}')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'Magnetic Field Strength ($B_\phi$)')
    plt.title(r'Diffusion of Magnetic field in azimuthal direction')
    # plt.ylim(np.min(solution_phi), np.max(solution_phi))
    plt.grid()
    plt.legend()

    # Create imshow plot
    plt.subplot(2, 2, 4)
    plt.contourf(*np.meshgrid(spatial_grid[1:-1], time_grid), pitch.T[:, 1:-1], 40, cmap='Spectral_r')
    plt.colorbar(label=r'Magnetic Field Strength ($B_\phi$)')
    plt.title('Diffusion of Magnetic field in azimuthal direction')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel('Time (Myr)')
    plt.grid()

    plt.tight_layout(pad=3)
