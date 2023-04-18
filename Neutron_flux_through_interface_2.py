# import of modules
import matplotlib.pyplot as plt
import numpy as np

# Without reflector
# Numerical solution

# initial data
# fuel parameters
Za_f = 1.3072 # neutron absorbing macroscopic cross section in cm^(-1)
vZf = 1.3727 # neutron fission macroscopic cross section in cm^(-1) * number of neutrons released per fission
Zf = 0.5649 # neutron fission macroscopic cross section in cm^(-1)
D_f = 0.1947 # diffusion coefficient in cm
Ef = 200 # amount of energy released per fission in MeV
af = 2.5 # half of the slab thickness in cm

# reflector parameters
Za_r = 0.0221 # neutron absorbing macroscopic cross section in cm^(-1)
D_r = 0.3278 # diffusion coefficient in cm
L = np.sqrt(D_r/Za_r) # diffusion length in cm
ar = 5 # reflector thickness in cm

al = af + ar # total thickness, cm
Q = 6.2414 * 10**14 # power in MeV/(s*cm^2)
del_al = 0.01 # grid step in cm
e = 10**(-6) # permissible error
D_int = (D_f + D_r)/2 # diffusion coefficient in the interface in cm

# coefficients
a_i = (D_f/(del_al**2) + Za_f)
a_end = (D_r/(del_al**2) + 1/(2.1312*del_al) + Za_r)
a_f = (2*D_f/(del_al**2) + Za_f)
a_r = (2*D_r/(del_al**2) + Za_r)
a_int_f = ((D_f + D_int)/(del_al**2) + Za_f)
a_int_r = ((D_r + D_int)/(del_al**2) + Za_r)
b_f = -D_f/del_al**2
b_r = -D_r/del_al**2
b_int = -D_int/del_al**2

# number of points or the order of the M Matrix
m = int(al/del_al)
m_f = int(af/del_al)
m_r = int(ar/del_al)

# M Matrix creation
M = np.zeros((m, m))

# M Matrix filling
M[0, 0] = float(a_i)
M[0, 1] = float(b_f)
M[1, 0] = float(b_f)
for i in range(1, m_f - 1):
    M[i, i] = float(a_f)
    M[i, i - 1] = float(b_f)
    M[i, i + 1] = float(b_f)
M[m_f - 1, m_f - 2] = float(b_f)
M[m_f - 1, m_f - 1] = float(a_int_f)
M[m_f - 1, m_f] = float(b_int)
M[m_f, m_f - 1] = float(b_int)
M[m_f, m_f] = float(a_int_r)
M[m_f, m_f + 1] = float(b_r)
M[m_f + 1, m_f] = float(b_r)
for i in range(m_f + 1, len(M)-1):
    M[i, i] = float(a_r)
    M[i, i - 1] = float(b_r)
    M[i, i + 1] = float(b_r)
M[len(M) - 1, len(M) - 2] = float(b_r)
M[len(M) - 1, len(M) - 1] = float(a_end)

np.savetxt('M1.txt', M, fmt='%.5e')

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
for i in range(0, m_f):
    VZf[i, i] = vZf

ek = 1
es = 1

while ek > e or np.all(es > e):
    # flux normalization
    N = (Q / 2) / (Ef * Zf * del_al * np.sum(F0[0:m_f]))
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
    for i in range(0, m_f):
        es = np.max(np.abs((S[i, 0] - S0[i, 0]) / S[i, 0]))
    F0 = F
    k0 = k

print('Numerical solution, k = ', k)

#Analytical solution
b_til = ar + 2.1312*D_r # extrapolated boundary in cm
Bm = 0.372
A = Q * Bm / (2 * Ef * Zf * np.sin(Bm * af))
B = A * np.cos(Bm * af) / np.sinh(b_til / L)
F2 = np.zeros((m_f, 1))
F3 = np.zeros((m_r, 1))
for i in range(0, len(F2)):
    F2[i, 0] = A * np.cos((Bm * del_al * i))
for i in range(0, len(F3)):
    F3[i, 0] = B * np.sinh((b_til - del_al * i)/L)
F_com = np.append(F2, F3)
print(B)

# Error
R = np.zeros((m, 1))
for i in range(0, len(R)):
    R[i, 0] = np.fabs(1 - (F[i] / F_com[i]))

# Plotting
x = np.linspace(del_al / 2, al - del_al / 2, m)
plt.ylabel('Flux, 1/(cm^2 * s)')
plt.xlabel('x, cm')
plt.grid(True)
plt.plot(x, F)
plt.plot(x, F_com)
plt.show()

x = np.linspace(del_al / 2, al - del_al / 2, m)
plt.ylabel('Relative error')
plt.xlabel('x, cm')
plt.grid(True)
plt.plot(x, R)
plt.show()