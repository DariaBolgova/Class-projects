# import of modules
import matplotlib.pyplot as plt
import numpy as np

# Numerical solution without reflector

# initial data
Za = 1.3072 # neutron absorbing macroscopic cross section in cm^(-1)
vZf = 1.3727 # neutron fission macroscopic cross section in cm^(-1) * number of neutrons released per fission
Zf = 0.5649 # neutron fission macroscopic cross section in cm^(-1)
D = 0.1947 # diffusion coefficient in cm
Ef = 200 # amount of energy released per fission in MeV
Q = 6.2414 * 10**14 # power in MeV/(s*cm^2)
al = 2.5 # half of the slab thickness in cm
del_al = 0.01 # grid step in cm
e = 10**(-6) # permissible error

# coefficients
a = (2* D/(del_al**2) + Za)
ai = (D/(del_al**2) + Za)
af = (D/(del_al**2) + 1/(2.1312*del_al) + Za)
b = -D/del_al**2

# number of points or the order of the M Matrix
m = int(al/del_al)

# M Matrix creation
M = np.zeros((m, m))

# creating an array filled with a certain value:
for i in range(0, 1):
    M[i, i] = float(ai)
    M[i, i + 1] = float(b)
for i in range(1, 2):
    M[i, i - 1] = float(b)
for i in range(len(M)-1, len(M)):
    M[i, i] = float(af)
    M[i, i - 1] = float(b)
for i in range(len(M)-2, len(M)-1):
    M[i, i + 1] = float(b)
for i in range(1, len(M)-1):
    M[i, i] = float(a)
    M[i, i - 1] = float(b)
    M[i, i + 1] = float(b)

# inverse M Matrix
M_inv = np.linalg.inv(M)

# initial flux and eigenvalue guess
f0 = 10**7
k0 = 1

# initial F matrix creation
F0 = np.zeros((m, 1))
for i in range(0, m):
    F0[i, 0] = float(f0)

# matrix of coefficients creation
VZf = np.zeros((m, m))
for i in range(0, m):
    VZf[i, i] = vZf

ek = 1
es = 1

while ek > e or np.all(es > e):
    # flux normalization
    N = (Q / 2) / (Ef * Zf * del_al * np.sum(F0))
    Fn = F0 * N

    # initial source matrix creation
    S0 = np.dot(VZf, Fn)

    # F matrix creation
    F = np.dot(np.dot(M_inv * 1/k0, VZf), Fn)

    # estimation of the source value
    S = np.dot(VZf, F)

    # estimation of k(n+1)
    k = k0 * np.sum(S) / np.sum(S0)

    # error estimation
    ek = np.abs((k - k0) / k)
    for i in range(0, m):
        es = np.max(np.abs((S[i, 0] - S0[i, 0]) / S[i, 0]))

    F0 = F
    k0 = k

print('Numerical solution, k = ', k)

#Analytical solution
a_til = al + 2.1312*D # extrapolated boundary in cm
Bm = np.pi / (2 * a_til)
A = Q * Bm / (2 * Ef * Zf * np.sin(Bm * al))

F2 = np.zeros((m, 1))
for i in range(0, len(F2)):
    F2[i, 0] = A * np.cos((np.pi * (del_al/2 + del_al * i)) / (2* a_til))

print(Bm)

k2 = vZf / (D * np.pi**2 / (4 * a_til**2) + Za)

print('Analytical solution, k = ', k2)

# Error
R = np.zeros((m, 1))

for i in range(0, len(R)):
    R[i, 0] = np.fabs(1 - (F[i] / F2[i]))

# Plotting
x = np.linspace(del_al / 2, al - del_al / 2, m)
plt.ylabel('Flux, 1/(cm^2 * s)')
plt.xlabel('x, cm')
plt.grid(True)
plt.plot(x, F)
plt.plot(x, F2)
plt.show()
x = np.linspace(del_al / 2, al - del_al / 2, m)
plt.ylabel('Relative error')
plt.xlabel('x, cm')
plt.grid(True)
plt.plot(x, R)
plt.show()