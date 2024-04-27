import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from my_code import *
from plotting import *


def init_cond_Br(x):
    return -(1-x**2) #(1-x**2)*np.cos(2*np.pi*x)

def init_cond_Bphi(x):
    return (1-x**2)

def source_term(x, t):
    return 0

z = np.linspace(-1, 1, 101)


# Constants and parameters
h = 7                               # 100pc
eta_T = 1.3#3.48e-2/h**2                     # magnetic diffusivity in (100pc)^2/Myr
Omega = 100*MYR*KM/(1000*PC)        # angular velocity, converted from km/s/kpc to 1/Myr
q = 0.98                            # shear parameter
t_max = 100                         # total simulation time
z_min = -1.0                        # minimum thickness of the disc
z_max = 1.0                         # thickness of the disc
dt = t_max/200                      # time step
dz = 0.01                           # spatial step in z direction
alpha_0 = 4/h#10*1e3*MYR/(100*PC)      # alpha effect, converted from km/s to 100pc/Myr

print('eta', eta_T, '\t',  'alpha',alpha_0, '\t', 'omega', Omega)
print('Rw', -q*Omega*h**2/eta_T)
print('Ra', alpha_0*h/eta_T)
print('Dynamo number ', -alpha_0*q*Omega*h**3/eta_T**2)


# Spatial grid
z = np.linspace(z_min, z_max, int((z_max - z_min) / dz) + 1)
t = np.linspace(0, t_max, int(t_max / dt) + 1)

# Coefficients for the matrix A and B
rho = eta_T*dt/(2*dz**2)
sigma = dt/(2*dz)
alpha = alpha_0*np.sin(np.pi*z/2)
dalpha_dt = np.gradient(alpha, z)

A = mod_matrix(z, 1+2*rho, 1, q*Omega*dt/2, 1+2*rho, -rho, sigma, 0, -rho, -rho, 0, 0, -rho, alpha, dalpha_dt)
B = mod_matrix(z, 1-2*rho, -1, -q*Omega*dt/2, 1-2*rho, rho, -sigma, 0, rho, rho, 0, 0, rho, alpha, dalpha_dt)

# Solve the diffusion equation in radial direction
solution = crank_nicolson_mod(len(z), len(t), init_cond_Br(z), init_cond_Bphi(z), A, B)

B_r = solution[:len(z), :]
B_phi = solution[len(z):, :]


# Create a figure and axis
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharex=True)

# Initialize the line objects for B_r and B_phi
line_br, = ax1.plot([], [], color='blue', label='$B_r$ (in $\mu B$)')
line_bphi, = ax2.plot([], [], color='red', label='$B_{\phi}$ (in $\mu B$)')



# Set the super title
fig.suptitle(r'Galactic magnetic field evolution of $B_r$ and $B_{\phi}$ for Dynamo number $D = $'+str(np.round(-alpha_0*q*Omega*h**3/eta_T**2, 4)))

# Set the axis limits
ax1.set_xlim(z_min, z_max)
ax1.set_ylim(np.min(B_r), np.max(B_r))
ax1.set_xlabel('$z$')
ax1.set_ylabel('Magnetic Field Strength $B_r$ (in $\mu G$)')
ax1.set_title('Animation of $B_r$')
ax1.grid()

ax2.set_xlim(z_min, z_max)
ax2.set_ylim(np.min(B_phi), np.max(B_phi))
ax2.set_xlabel('$z$')
ax2.set_ylabel('Magnetic Field Strength $B_{\phi}$ (in $\mu G$)')
ax2.set_title('Animation of $B_{\phi}$')
ax2.grid()

# Create the update function for the animation
def update(frame):
    # Clear the previous lines
    line_br.set_data(z, B_r[:, frame])
    line_bphi.set_data(z, B_phi[:, frame])
    
    # Set the legend with changing time and LaTeX form
    ax1.legend([r'$B_r$'+f'\nTime: ${t[frame]:.2f}$'])
    ax2.legend([r'$B_{\phi}$'+f'\nTime: ${t[frame]:.2f}$'])
    
    return line_br, line_bphi

# Create the animation
animation = FuncAnimation(fig, update, frames=len(t), interval=50, blit=True)

# Set tight layout
plt.tight_layout()

# Display the animation
animation.save('ani_1.gif', writer='pillow')

plt.show()
