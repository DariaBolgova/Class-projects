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
iter = 0 #iteration number counter
e = 0
L = 2
LL = 1

#Setting the initial velocity
print('Please enter the initial velocity value, m/s:')
uin = float(input())

#initial x-component velocity field
U = np.concatenate((np.full((1, Nx+4), -uin), np.full((Ny, Nx+4), uin), np.full((1, Nx+4), -uin)), axis=0)

#arrays creation
A = np.zeros([Ny, Nx+2]) #advection term
D = np.zeros([Ny, Nx+2]) #diffusion term
US = np.zeros([Ny+2, Nx+4]) #temporary x-component velosity (U*)
UU = np.zeros([Ny+2, Nx+4]) #velosity field at n+1 time step
U_an = np.zeros([Ny, 1]) #analytical solutuion
y = np.zeros([Ny, 1]) #y-axis coordinates
Pl = np.zeros([Ny, 1]) #numerical solution plot
dU = np.zeros([Ny, 1]) #difference between numerical and analytical solution

#analytical solutuion
for j in range(0, Ny):
    y[j,0] = (hy/2 + j*hy)
    U_an[j, 0] = - 187.5 * (y[j,0] - Ly/2) ** 2 + 0.01875

#numerical solution
while L > 0.000001:
    iter = iter + 1
    # advection term at n time step
    for j in range(0, Ny):
        for i in range(0, Nx + 2):
            A[j, i] = (1 / hx) * (((U[j + 1, i] + U[j + 1, i + 1]) / 2) ** 2 - ((U[j + 1, i + 1] + U[j + 1, i + 2]) / 2) ** 2)

    # diffusion term at n time step
    for j in range(0, Ny):
        for i in range(0, Nx + 2):
            D[j, i] = (1 / hx ** 2) * (U[j + 1, i] + U[j + 1, i + 2] - 2 * U[j + 1, i + 1]) + (1 / hy ** 2) * (U[j, i + 1] + U[j + 2, i + 1] - 2 * U[j + 1, i + 1])

    # temporary (U-star) x-component velosity
    for j in range(1, Ny + 1):
        for i in range(1, Nx + 3):
            US[j, i] = U[j, i] + t * (-A[j - 1, i - 1] + vis * D[j - 1, i - 1])
            US[0, i] = -US[1, i]
            US[Ny + 1, i] = -US[Ny, i]
    for j in range(0, Ny + 2):
        US[j, 0] = US[j, 1]
        US[j, Nx + 3] = US[j, Nx + 2]

    # x-component velosity field at n+1 time step
    for j in range(1, Ny + 1):
        for i in range(1, Nx + 3):
            UU[j, i] = US[j, i] - (t * p/(rho))
            UU[0, i] = -UU[1, i]
            UU[Ny + 1, i] = -UU[Ny, i]
    for j in range(0, Ny + 2):
        UU[j, 0] = UU[j, 1]
        UU[j, Nx + 3] = UU[j, Nx + 2]

    #velocity field reassignment before calculation the new time step
    for j in range(0, Ny + 2):
        for i in range(0, Nx + 4):
            U[j, i] = UU[j, i]

    for j in range(0, Ny):
        Pl[j, 0] = U[j + 1, round(Nx / 2)]
        dU[j, 0] = Pl[j, 0] - U_an[j, 0]
    LL = ((1 / Ny) * (sum(dU) ** 2)) ** 0.5
    e = LL / L
    L = LL


print ('Velosity distribution in nodes in y direction is ', Pl)
print('Number of inerations is ', iter)
print(e)
print(L)
plt.plot(U_an, y)
plt.plot(Pl, y)
plt.title("Velocity profile")
plt.ylabel('u(y)')
plt.xlabel('y')
plt.show()

