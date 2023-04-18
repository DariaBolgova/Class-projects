import numpy as np
import matplotlib.pyplot as plt
import math as m

#Domain parameters and flow conditions
L = 0.015 #domain lenght, m
N = int(80) #number of elements
h = L/N #element lenght, m
rhol = 1000 #liquid density, kg/m^3
rhog = 1 #gas density, kg/m^3
M = 3 #interface thickness
m_d_an = rhol * L/2 #analytical mass of the droplet

#arrays creation
Phi = np.zeros(N) #level-set distance field
F = np.zeros(N) #smoothed Heaviside function
Rho = np.zeros(N) #density field
H = np.zeros(N) #cell numbers (for plotting)

M_d_1 = np.zeros(N) #numerical mass of the liquid object, option 1
M_d_2 = np.zeros(N) #numerical mass of the liquid object, option 2
M_d_3 = np.zeros(N) #numerical mass of the liquid object, option 3
DW = np.zeros(N) #density weight

for i in range(0, N):
    H[i] = i+1

#1D level-set distance field
for i in range(0, int(N/2)):
    Phi[i] = (i - (N/4 - 1/2)) * h
for i in range(int(N/2), int(N)):
    Phi[i] = (-i + (N*3/4 - 1/2)) * h

#1D density field
for i in range(0, N):
    if Phi[i] < (-M*h):
        F[i] = 0
    if abs(Phi[i]) <= (M*h):
        F[i] = 0.5 * (1 + Phi[i]/(M*h) + (1/m.pi) * m.sin(m.pi * Phi[i]/(M*h)))
    if Phi[i] > (M*h):
        F[i] = 1

#density field
for i in range(0, N):
    Rho[i] = rhol * F[i] + rhog*(1 - F[i])

#numerical mass calculation
#option 1
for i in range(0, N):
    if Phi[i] > 0:
        M_d_1[i] = (rhol * h)
m_d_1 = sum(M_d_1)
RE1 = abs((m_d_an - m_d_1)*100/m_d_an)

#option 2
for i in range(0, N):
    if Phi[i] > 0:
        M_d_2[i] = (Rho[i] * h)
m_d_2 = sum(M_d_2)
RE2 = abs((m_d_an - m_d_2)*100/m_d_an)

#option 3
for i in range(0, N):
    DW[i] = (Rho[i] - rhog)/(rhol - rhog)
    M_d_3[i] = (rhol * DW[i] * h)
m_d_3 = sum(M_d_3)
RE3 = abs((m_d_an - m_d_3)*100/m_d_an)

plt.plot(H, Rho)
plt.title("Density field, M = 3")
plt.ylabel('Rho, kg/m^3')
plt.xlabel('N')
plt.show()

print("The analytical mass of the liquid object is m = ", round(m_d_an, 3))
print("The numerical mass of the liquid object is:")
print("Option 1: m = ", round(m_d_1, 3),", RE = ", round(RE1, 2), "%")
print("Option 2: m = ", round(m_d_2, 3),", RE = ", round(RE2, 2), "%")
print("Option 3: m = ", round(m_d_3, 3),", RE = ", round(RE3, 2), "%")