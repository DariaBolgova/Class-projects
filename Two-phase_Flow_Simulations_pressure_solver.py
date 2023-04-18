import numpy as np
import matplotlib.pyplot as plt
from itertools import repeat

#Desired number of nodes across the domain input
print('Please enter the number of nodes across the domain:')
Ny = int(input())
Nx = Ny * 2

#Domain parameters and flow conditions
Ly = 0.02 #Domain height, m
Lx = 0.04 #Domain lenght, m
hy = Ly/Ny #Cell size in y direction, m
hx = Lx/Nx #Cell size in x direction, m

uav = 0.0125 #average velocity, m/s
umax = 0.01875 #maximum velocity, m/s
p = -0.375 #pressure gradient, Pa/m
vis = 0.000001 #kinematic viscosity, m^2/s
rho = 1000 #density, kg/m^3
t = 0.1 #time step, s
iter = 0 #iteration number counter for velocity solution
iterp = 0 #iteration number counter for pressure solution
L = 2 #velocity solution error on n time step
LL = 1 #velocity solution error on n+1 time step
eps = 1 #pressure convergence criteria

#Setting the initial velocity
print('Please enter the initial velocity value, m/s:')
uin = float(input())

#Pressure assumption
print('Please enter the initial pressure value assumption, Pa:')
pin = float(input())

#initial x-component velocity field
U = np.concatenate((np.full((1, Nx+3), -uin), np.full((Ny, Nx+3), uin), np.full((1, Nx+3), -uin)), axis=0)

#arrays creation
A = np.zeros([Ny, Nx+1]) #advection term
D = np.zeros([Ny, Nx+1]) #diffusion term
US = np.zeros([Ny+2, Nx+3]) #temporary x-component velosity (U*)
P = np.full((Ny, Nx+2), pin) #pressure field at n time step
PP = np.zeros([Ny, Nx+2]) #pressure field at n+1 time step
ep = np.full((Ny, Nx+2), 1) #pressure solution error on n+1 time step
UU = np.zeros([Ny+2, Nx+3]) #velosity field at n+1 time step
U_an = np.zeros([Ny, 1]) #analytical solutuion
y = np.zeros([Ny, 1]) #y-axis coordinates
Pl = np.zeros([Ny, 1]) #numerical solution plot
dU = np.zeros([Ny, 1]) #difference between numerical and analytical solution

#analytical solutuion
for j in range(0, Ny):
    y[j,0] = (hy/2 + j*hy)
    U_an[j, 0] = - 187.5 * (y[j,0] - Ly/2) ** 2 + 0.01875

#numerical solution
while L > 0.0001:
    iter = iter + 1
    # advection term at n time step
    for j in range(0, Ny):
        for i in range(0, Nx + 1):
            A[j, i] = (1 / hx) * (((U[j + 1, i] + U[j + 1, i + 1]) / 2) ** 2 - ((U[j + 1, i + 1] + U[j + 1, i + 2]) / 2) ** 2)

    # diffusion term at n time step
    for j in range(0, Ny):
        for i in range(0, Nx + 1):
            D[j, i] = (1 / hx ** 2) * (U[j + 1, i] + U[j + 1, i + 2] - 2 * U[j + 1, i + 1]) + (1 / hy ** 2) * (U[j, i + 1] + U[j + 2, i + 1] - 2 * U[j + 1, i + 1])

    # temporary (U-star) x-component velosity
    for j in range(1, Ny + 1):
        for i in range(1, Nx + 2):
            US[j, i] = U[j, i] + t * (-A[j - 1, i - 1] + vis * D[j - 1, i - 1] - p/rho)
            US[0, i] = -US[1, i]
            US[Ny + 1, i] = -US[Ny, i]
    for j in range(0, Ny + 2):
        US[j, 0] = US[j, 1]
        US[j, Nx + 2] = US[j, Nx + 1]

    #pressure
    while eps > 0.0001:
        iterp = iterp + 1

        for j in range(1, Ny-1):
            for i in range(1, Nx+1):
                RHS = (rho/(t * hx)) * (US[j+1, i+2] - US[j+1, i+1])
                PP[j, i] = -(1/4) * (RHS * hx**2 - (P[j, i+1] + P[j, i-1] + P[j+1, i] + P[j-1, i]))
                PP[0, i] = -(1/3) * (RHS * hx**2 - (P[j, i+1] + P[j, i-1] + P[j+1, i]))
                PP[Ny-1, i] = -(1/3) * (RHS * hx**2 - (P[j, i+1] + P[j, i-1] + P[j-1, i]))
        for j in range(0, Ny):
            PP[j, 0] = PP[j, 1]
            PP[j, Nx+1] = PP[j, Nx]

        for j in range(0, Ny):
            for i in range(0, Nx+2):
                ep[j, i] = PP[j, i] - P[j, i]
        eps = ep.max()
        for j in range(0, Ny):
            for i in range(0, Nx+2):
                P[j, i] = PP[j, i]

    # x-component velosity field at n+1 time step
    for j in range(1, Ny + 1):
        for i in range(1, Nx + 2):
            UU[j, i] = US[j, i] - (t/(rho * hx)) * (PP[j-1, i-1] - PP[j-1, i-2])
            UU[0, i] = -UU[1, i]
            UU[Ny + 1, i] = -UU[Ny, i]
    for j in range(0, Ny + 2):
        UU[j, 0] = UU[j, 1]
        UU[j, Nx + 2] = UU[j, Nx + 1]

    #velocity field reassignment before calculation the new time step
    for j in range(0, Ny + 2):
        for i in range(0, Nx + 3):
            U[j, i] = UU[j, i]

    for j in range(0, Ny):
        Pl[j, 0] = U[j + 1, round(Nx / 2)]
        dU[j, 0] = Pl[j, 0] - U_an[j, 0]
    LL = ((1 / Ny) * (sum(dU) ** 2)) ** 0.5
    L = LL

print(PP)
print ('Velosity distribution in nodes in y direction is ', Pl)
print('Number of inerations is ', iter)
print('The solution error is ', L)
plt.plot(U_an, y)
plt.plot(Pl, y)
plt.title("Velocity profile")
plt.ylabel('u(y)')
plt.xlabel('y')
plt.show()

