# import of modules for plotting
import matplotlib.pyplot as plt
import numpy as np

#initial data
f0 = 10**12 # flux density at x = 0, cm^(-2)*s^(-1)
Za = 1.25 # neutron absorbing macroscopic cross section, cm^(-1)
al = 2.5 # slab thickness, cm
del_al = 0.01 # grid step, cm
D = 1/(3*Za) # diffusion coefficient, cm

#Matrix coefficients
a = (2* D/(del_al**2) + Za)
ai = (3*D/(del_al**2) + Za)
af = (D/(del_al**2/2 + 2.1312*del_al*D) + D/(del_al**2) + Za)
b = -D/del_al**2
s = 2*D*f0/del_al**2

# number of points (order of the M Matrix)
m = int(al/del_al)

# M Matrix
M = np.zeros((m,m))

# M Matrix filling
for i in range(0,1):
    M[i, i] = float(ai)
    M[i, i + 1] = float(b)
for i in range(1,2):
    M[i, i - 1] = float(b)
for i in range(len(M)-1,len(M)):
    M[i, i] = float(af)
    M[i, i - 1] = float(b)
for i in range(len(M)-2,len(M)-1):
    M[i, i + 1] = float(b)
for i in range(1,len(M)-1):
    M[i, i] = float(a)
    M[i, i - 1] = float(b)
    M[i, i + 1] = float(b)
np.savetxt('result.txt', M, fmt='%.2e')

# inverse M Matrix
M_inv = np.linalg.inv(M)

# S matrix
S = np.zeros((m,1))

# S Matrix filling
for i in range(0,1):
    S[i, 0] = float(s)

# F matrix (solution)
F = np.dot(M_inv, S)

# Analytical solution
e = al + 2.1312*D
L = np.sqrt(D/Za)

F2 = np.zeros((m,1))
for i in range(0,len(F2)):
    F2[i, 0] = f0 * (np.exp((([i]) / (L * 100)) + 0.005/L) - np.exp(e / L) * np.sinh((([i]) / (L * 100)) + 0.005/L) / np.sinh(e / L))

# Error
R = np.zeros((m,1))
for i in range(0,len(R)):
    R[i, 0] = np.fabs(1 - (F[i] / F2[i]))

# Plotting
x = np.linspace(del_al/2, al - del_al/2, 250)

fig, ax = plt.subplots()
ax.plot(x, F, color="blue", label="f(x)_num")
ax.plot(x, F2, color="green", label="f(x)_an")

fig, ar = plt.subplots()
ar.plot(x, R, color="red", label="R")
plt.show()

