import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def create_animation(B_array, z_values, t_values, filename='animation.gif', B_label='B(z)', z_label='z'):
    fig, ax = plt.subplots()

    def update(frame):
        ax.clear()
        ax.plot(z_values, B_array[:, frame])
        ax.set_title(f"Time = {t_values[frame]:.2f}")
        ax.set_xlabel(z_label)
        ax.set_ylabel(B_label)
        ax.set_xlim(z_values.min(), z_values.max())
        ax.set_ylim(B_array.min(), B_array.max())

    ani = animation.FuncAnimation(fig, update, frames=len(t_values), interval=100)

    ani.save(filename, writer='pillow')
    plt.close(fig)


def plot_init_cond(z, init_cond_Br, init_cond_Bphi, title1, title2, global_title):

    plt.figure(figsize=(11, 3.5))
    plt.subplot(121)
    plt.plot(z, init_cond_Br(z))
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'$B_r$')
    plt.title(title1)

    plt.subplot(122)
    plt.plot(z, init_cond_Bphi(z))
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'$B_{\phi}$')
    plt.title(title2)

    plt.suptitle(global_title)
    plt.tight_layout(pad=1)


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
    plt.contourf(*np.meshgrid(spatial_grid, time_grid), solution_r.T, 50, cmap='Spectral_r')
    plt.colorbar(label=r'Magnetic Field Strength ($B_r$)')
    plt.title(r'Diffusion of Magnetic field in radial direction')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'Time (Myr)')


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
    plt.contourf(*np.meshgrid(spatial_grid, time_grid), solution_phi.T, 50, cmap='Spectral_r')
    plt.colorbar(label=r'Magnetic Field Strength ($B_\phi$)')
    plt.title('Diffusion of Magnetic field in azimuthal direction')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel('Time (Myr)')
    plt.tight_layout(pad=3)



def plot_pitch(time_grid, spatial_grid, B, pitch):

    # Create 2D plots
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 2, 1)
    for i in (range(0, len(time_grid), int(len(time_grid)/5))):
        plt.plot(spatial_grid, B[:, i], label=f'time = {time_grid[i]:.1f}')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'$B_{total}$')
    plt.title('Diffusion of total magnetic field')
    plt.grid()
    plt.legend()

    # Create imshow plot
    plt.subplot(2, 2, 2)
    plt.contourf(*np.meshgrid(spatial_grid, time_grid), B.T, 40, cmap='Spectral_r')
    plt.colorbar(label=r'($B_{total}$)')
    plt.title(r'Diffusion of total magnetic field')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'Time (Myr)')


    # Create 2D plots
    plt.subplot(2, 2, 3)
    for i in (range(0, len(time_grid), int(len(time_grid)/5))):
        plt.plot(spatial_grid[1:-1], pitch[1:-1, i], label=f'time = {time_grid[i]:.1f}')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel(r'Pitch angle $p_B$ (in degrees)')
    plt.title(r'Variation of pitch angle with time')
    plt.grid()
    plt.legend()

    # Create imshow plot
    plt.subplot(2, 2, 4)
    plt.contourf(*np.meshgrid(spatial_grid[1:-1], time_grid), pitch.T[:, 1:-1], 40, cmap='Spectral_r')
    plt.colorbar(label=r'Pitch angle $p_B$ (in degrees)')
    plt.title('Variation of pitch angle with time')
    plt.xlabel(r'$z$ (normalized to 100 pc)')
    plt.ylabel('Time (Myr)')

    plt.tight_layout(pad=3)



def plot_decay(time_grid, B_mid, m, c):
    # Plot the log of magnetic field strength at midplane and the slope of the logplot
    plt.figure(figsize=(6, 4))
    plt.plot(time_grid, B_mid, 'b-')
    # plot another line with the slope and intercept m and c
    plt.plot(time_grid[-50:], m*time_grid[-50:] + c, 'r:', linewidth=3, label=r'Slope ($\gamma$) = {:.3e}'.format(m))
    plt.xlabel('Time (Myr)')
    plt.ylabel('log$(B_{total})$ at midplane')
    plt.title(r'Magnetic field strength at midplane')
    # plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.tight_layout()
