import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from my_code import *
from plotting import *


def init_cond_S(x):
    return -(1-x**2) #(1-x**2)*np.cos(2*np.pi*x)

def init_cond_T(x):
    return (1-x**2)

def source_term(x, t):
    return 0

z = np.linspace(-1, 1, 101)


# Constants and parameters
h = 5                               # 100pc
eta_T = 0.7#3.48e-2/h**2                     # magnetic diffusivity in (100pc)^2/Myr
Omega = 40*MYR*KM/(1000*PC)        # angular velocity, converted from km/s/kpc to 1/Myr
q = 0.98                            # shear parameter
t_max = 10                          # total simulation time
z_min = -1.0                        # minimum thickness of the disc
z_max = 1.0                         # thickness of the disc
dt = t_max/200                      # time step
dz = 0.01                           # spatial step in z direction
alpha_0 = 1/h#10*1e3*MYR/(100*PC)      # alpha effect, converted from km/s to 100pc/Myr

# print('eta', eta_T, '\t',  'alpha',alpha_0, '\t', 'omega', Omega)
# print('Rw', -q*Omega*h**2/eta_T)
# print('Ra', alpha_0*h/eta_T)
# print('Dynamo number ', -alpha_0*q*Omega*h**3/eta_T**2)

# Spatial grid
z = np.linspace(z_min, z_max, int((z_max - z_min) / dz) + 1)
t = np.linspace(0, t_max, int(t_max / dt) + 1)

# Coefficients for the matrix A and B
nu = eta_T*dt/(2*dz**2)
mu = dt/(2*dz)
alpha = alpha_0*np.sin(np.pi*z/2)
dalpha_dt = np.gradient(alpha, z)

A = mod_matrix(len(z), 1+2*nu, dz*mu, q*Omega*mu, 1+2*nu, -nu, 0, -q*Omega*mu, -nu, -nu, 0, 0, -nu, alpha, dalpha_dt, dz, dt)
B = mod_matrix(len(z), 1-2*nu, -dz*mu, -q*Omega*mu, 1-2*nu, nu, 0, q*Omega*mu, nu, nu, 0, 0, nu, alpha, dalpha_dt, dz, dt)

# Solve the diffusion equation in radial direction
solution = crank_nicolson_mod(len(z), len(t), init_cond_S(z), init_cond_T(z), A, B)

S = solution[:len(z), :]
T = solution[len(z):, :]


# Create a figure and axis
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharex=True)

# Initialize the line objects for B_r and B_phi
line_br, = ax1.plot([], [], color='blue', label='$\psi$ (in $\mu G$/100pc)')
line_bphi, = ax2.plot([], [], color='red', label='$T$ (in $\mu G$/100pc)')



# Set the super title
fig.suptitle(r'Galactic magnetic field evolution of $\psi$ and $T$ for Dynamo number $D = $'+str(np.round(-alpha_0*q*Omega*h**3/eta_T**2, 4)))

# Set the axis limits
ax1.set_xlim(z_min, z_max)
ax1.set_ylim(np.min(S), np.max(S))
ax1.set_xlabel('$z$')
ax1.set_ylabel('Strength of $\psi$ (in $\mu G$/100pc)')
ax1.set_title('Animation of $\psi$')
ax1.grid()

ax2.set_xlim(z_min, z_max)
ax2.set_ylim(np.min(T), np.max(T))
ax2.set_xlabel('$z$')
ax2.set_ylabel('Magnetic Field Strength $T$ (in $\mu G$/100pc)')
ax2.set_title('Animation of $T$')
ax2.grid()

# Create the update function for the animation
def update(frame):
    # Clear the previous lines
    line_br.set_data(z, S[:, frame])
    line_bphi.set_data(z, T[:, frame])
    
    # Set the legend with changing time and LaTeX form
    ax1.legend([r'$\psi$'+f'\nTime: ${t[frame]:.2f}$'])
    ax2.legend([r'$T$'+f'\nTime: ${t[frame]:.2f}$'])
    
    return line_br, line_bphi

# Create the animation
animation = FuncAnimation(fig, update, frames=len(t), interval=50, blit=True)

# Set tight layout
plt.tight_layout()

# Display the animation
animation.save('ani_3.gif', writer='pillow')

plt.show()
