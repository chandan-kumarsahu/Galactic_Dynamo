# Importing required libraries
import matplotlib.pyplot as plt
import numpy as np

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


# Pitch angle
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




def matrix_A(N, a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4):
    # return a 2N x 2N matrix with the given terms in each block of the matrix
    A = np.zeros((2*N, 2*N))
    for i in range(N):
        A[i, i] = a1
        A[i, i+N] = a2
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

def matrix_B(N, a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4):
    B = np.zeros((2*N, 2*N))
    for i in range(N):
        B[i, i] = a1
        B[i, i+N] = a2
        B[i+N, i] = a3
        B[i+N, i+N] = a4
    for i in range(N-1):
        B[i, i+1] = b1
        B[i, i+N+1] = b2
        B[i+N, i+1] = b3
        B[i+N, i+N+1] = b4
        B[i+1, i] = c1
        B[i+1, i+N] = c2
        B[i+N+1, i] = c3
        B[i+N+1, i+N] = c4
    return B



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

