# Importing required libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

##########################################################################################################

# CONSTANTS
KM = 1e3
PC = 3.086e16
MYR = 1e6*365*24*3600

##########################################################################################################

"""
Get the matrices A and B for solving the diffusion equation using Crank-Nicolson method.
This function is used for vacuum boundary conditions.

Parameters:
- N: Number of spatial grid points
- sigma: alpha*dt/dx^2

Returns:
- A: Matrix A
- B: Matrix B
"""

def diff_matrix_vacuum_boundary(N, sigma):
    # Initialize matrices A and B with zeros
    A = [[0] * N for _ in range(N)]
    B = [[0] * N for _ in range(N)]

    # Interior points
    for i in range(0, N):
        A[i][i] = 2 + 2 * sigma  # Diagonal element of A
        B[i][i] = 2 - 2 * sigma  # Diagonal element of B

        # Connect to the left neighbor (if not on the left edge)
        if i > 0:
            A[i][i - 1] = -sigma
            B[i][i - 1] = sigma

        # Connect to the right neighbor (if not on the right edge)
        if i < N - 1:
            A[i][i + 1] = -sigma
            B[i][i + 1] = sigma

    return A, B




"""
Get the matrices A and B for solving the diffusion equation using Crank-Nicolson method.
This function is used for isolated boundary conditions.

Parameters:
- N: Number of spatial grid points
- sigma: alpha*dt/dx^2

Returns:
- A: Matrix A
- B: Matrix B
"""
def diff_matrix_isolated_boundary(N, sigma):
    # Initialize matrices A and B with zeros
    A = [[0] * N for _ in range(N)]
    B = [[0] * N for _ in range(N)]

    # Fill diagonal and off-diagonal values for matrices A and B
    for i in range(N):
        A[i][i] = 2 + 2 * sigma  # Diagonal element of A
        B[i][i] = 2 - 2 * sigma  # Diagonal element of B

        # Connect to the left neighbor (if not on the left edge)
        if i > 0:
            A[i][i - 1] = -sigma
            B[i][i - 1] = sigma

        # Connect to the right neighbor (if not on the right edge)
        if i < N - 1:
            A[i][i + 1] = -sigma
            B[i][i + 1] = sigma

    # Boundary conditions
    A[0][0] = 2 + sigma
    B[0][0] = 2 - sigma
    A[-1][-1] = 2 + sigma
    B[-1][-1] = 2 - sigma

    return A, B




"""
Solve 1D diffusion equation using Crank-Nicolson method.

Parameters:
- x_max: Extent of the spatial domain
- t_max: Total simulation time
- dx: Spatial step size
- dt: Time step size
- Diff: Thermal diffusivity
- init_cond: Initial condition function
- source_term: Source term function
- boundary: Boundary condition function

Returns:
- u: Temperature distribution over space and time
- x: Spatial grid
- t: Time grid
"""
def crank_nicolson_diffusion(x_min, x_max, t_max, dx, dt, Diff, init_cond, source_term, boundary):

    alpha = Diff * dt / (dx**2)

    # Spatial grid
    x = np.linspace(x_min, x_max, int((x_max - x_min) / dx) + 1)
    t = np.linspace(0, t_max, int(t_max / dt) + 1)

    # Initialize temperature array
    Temp = np.zeros((len(x), len(t)))

    # Initial condition
    for i in range(len(x)):
        Temp[i][0] = init_cond(x[i])

    # Get the matrices for solving the matrix using crank-nicolson method
    A, B = boundary(len(x), alpha)

    A = np.array(A)
    B = np.array(B)

    for j in range(1, len(t)):
        source_vector = np.array([source_term(xi, t[j]) for xi in x])
        Temp[:, j] = np.linalg.solve(A, np.dot(B, Temp[:, j - 1]) + dt * source_vector)

    return Temp, np.array(x), np.array(t)




"""
Calculation of the pitch angle

Parameters:
- Br: Radial component of the magnetic field
- Bphi: Azimuthal component of the magnetic field

Returns:
- B: Total magnetic field
- p: Pitch angle
"""
def get_B_and_pitch(Br, Bphi):
    B = np.sqrt(Br**2 + Bphi**2)
    p = np.zeros(Br.shape)
    for i in range(Br.shape[0]):
        for j in range(Br.shape[1]):
            if Bphi[i, j]!=0:
                p[i, j] = 180/np.pi*np.arctan(Br[i, j]/Bphi[i, j])
            elif Br[i, j]>0:
                p[i, j] = 90
            elif Br[i, j]<0:
                p[i, j] = -90
            else:
                p[i, j] = 0
    return B, p




"""
Get the matrices A and B for solving the diffusion equation using Crank-Nicolson method.
This function is used for vacuum boundary conditions.

Parameters:
- N: Number of spatial grid points
- a1, b1, ... etc: Coefficients of the matrix
- alpha: Alpha effect
- dalpha_dt: Derivative of alpha with respect to time

Returns:
- A: Matrix A
"""
def mod_matrix(N, a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, alpha, dalpha_dt):
    A = np.zeros((2*N, 2*N))
    for i in range(N):
        A[i, i] = a1
        A[i, i+N] = (dalpha_dt[i]/2 - b2*alpha[i])*a2
        A[i+N, i] = a3
        A[i+N, i+N] = a4
    for i in range(N-1):
        A[i, i+1] = b1
        A[i, i+N+1] = b2
        A[i+N, i+1] = b3
        A[i+N, i+N+1] = b4
        A[i+1, i] = c1
        A[i+1, i+N] = c2
        A[i+N+1, i] = c3
        A[i+N+1, i+N] = c4
    return A




"""
Solve 1D diffusion equation using Crank-Nicolson method with modified matrices.

Parameters:
- N_x: Number of spatial grid points
- N_t: Number of time grid points
- init_cond_Br: Initial condition for Br
- init_cond_Bphi: Initial condition for Bphi
- A: Coefficient matrix A
- B: Coefficient matrix B

Returns:
- U: Magnetic field distribution over space and time
"""
def crank_nicolson_mod(N_x, N_t, init_cond_Br, init_cond_Bphi, A, B):

    # Initialize temperature array
    U = np.zeros((2*N_x, N_t))

    # Initial condition
    for i in range(N_x):
        U[i, 0] = init_cond_Br[i]
        U[N_x+i, 0] = init_cond_Bphi[i]

    for j in range(1, N_t):
        U[:, j] = np.dot(np.linalg.inv(A), np.dot(B, U[:, j - 1]))

    return U




"""
Function to find the local maxima of a curve.

Parameters:
- x: x values
- y: y values

Returns:
- x_maxima: x values of the local maxima
- y_maxima: y values of the local maxima
"""
def find_local_maxima(x, y):
    x = x[100:]
    y = y[100:]
    x_maxima = []
    y_maxima = []
    for i in range(1, len(y) - 1):
        if y[i] > y[i - 1] and y[i] > y[i + 1]:
            x_maxima.append(x[i])
            y_maxima.append(y[i])
    return np.array(x_maxima), np.array(y_maxima)




"""
Function to find the decay rate of a curve.

Parameters:
- x: x values
- y: y values

Returns:
- slope: Decay rate of the curve
"""
def get_decay_rate(x, y):
    # x, y = find_local_maxima(x, y)
    y = np.log(y)
    slope, intercept = np.polyfit(x, y, 1)
    return slope




"""
Function to find roots of a function using bisection method.

Parameters:
- f: Function for which roots are to be found
- a: Lower bound of the interval
- b: Upper bound of the interval
- eps: Desired accuracy

Returns:
- c: Root of the function
"""
def bisection(f, a, b, eps):
    counter = 1
    COUNT = []
    VAL = []
    if f(a)*f(b) == 0.0:
        if f(a)==0.0:
            return a
        else:
            return b

    c = (a+b)/2
    while np.abs(f(c)) > eps: # checking if the accuracy is achieved

        c = (a+b)/2
        if (f(a)*f(c)) <= 0.0: # Check if the root is properly bracketted
            b = c
        else:
            a = c
        if counter > 100:
            print('Maximum iterations reached.')
            break
        counter += 1
        COUNT.append(counter)
        VAL.append(c)
        print(np.round(c,6), round(f(c),6))

    return c, COUNT, VAL
