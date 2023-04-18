# import of modules for plotting
import matplotlib.pyplot as plt
import numpy as np

# initial data
D_fuel_fast = 0.7936
D_mod_fast = 0.9493
D_fuel_th = 0.1025
D_mod_th = 0.2782
D_int_fast = (D_fuel_fast + D_mod_fast) / 2
D_int_th = (D_fuel_th + D_mod_th) / 2

ZR_fuel_fast = 0.0232
ZR_mod_fast = 0.0451
ZR_fuel_th = 2.8301
ZR_mod_th = 0.0669

vZf_fuel_fast = 0.0368
vZf_fuel_th = 2.8243

Zf_fast = 0.0139
Zf_th = 1.1623

Zs_12_fuel = 0.0016
Zs_12_mod = 0.045

r_fuel = 1
r_mod = 1
r_tot = r_fuel + r_mod
dr = 0.005

e = 10**(-6)

Q = 3.121 * 10**14
Ef = 200

# number of points or the order of the M Matrix
m = int(r_tot/dr)
m_f = int(r_fuel/dr)
m_r = int(r_mod/dr)

# M Matrix creation
M = np.zeros((2 * m, 2 * m))

# M Matrix filling
for i in range(1, m_f-1):
    M[i, i] = float((2 * D_fuel_fast)/dr**2 + ZR_fuel_fast)
    M[i, i - 1] = float(-D_fuel_fast/dr**2 + D_fuel_fast/(2 * dr * (dr * i + dr/2)))
    M[i, i + 1] = float(-D_fuel_fast/dr**2 - D_fuel_fast/(2 * dr * (dr * i + dr/2)))

M[0, 0] = float(D_fuel_fast/dr**2 + D_fuel_fast/(2 * dr * (dr/2)) + ZR_fuel_fast)
M[0, 1] = float(-D_fuel_fast/dr**2 - D_fuel_fast/(2 * dr * (dr/2)))
M[1, 0] = float(-D_fuel_fast/dr**2 + D_fuel_fast/(2 * dr * (2*dr - dr/2)))

M[m_f - 1, m_f - 2] = float(-D_fuel_fast/dr**2 + D_fuel_fast/(2 * dr * (dr * (m_f - 1) + dr/2)))
M[m_f - 1, m_f - 1] = float(D_int_fast/dr**2 + D_int_fast/(2 * dr * (dr * (m_f - 1) + dr/2)) + D_fuel_fast/dr**2 - D_fuel_fast/(2 * dr * (dr * (m_f - 1) + dr/2)) + ZR_fuel_fast)
M[m_f - 1, m_f] = float(-D_int_fast/dr**2 - D_int_fast/(2 * dr * (dr * (m_f - 1) + dr/2)))
M[m_f, m_f - 1] = float(-D_int_fast/dr**2 + D_int_fast/(2 * dr * (dr * m_f + dr/2)))
M[m_f, m_f] = float(D_mod_fast/dr**2 + D_mod_fast/(2 * dr * (dr * m_f + dr/2)) + D_int_fast/dr**2 - D_int_fast/(2 * dr * (dr * m_f + dr/2)) + ZR_mod_fast)
M[m_f, m_f + 1] = float(-D_mod_fast/dr**2 - D_mod_fast/(2 * dr * (dr * m_f + dr/2)))
M[m_f + 1, m_f] = float(-D_mod_fast/dr**2 + D_mod_fast/(2 * dr * (dr * (m_f + 1) + dr/2)))

for i in range(m_f + 1, m - 1):
    M[i, i] = float(2 * D_mod_fast/dr**2 + ZR_mod_fast)
    M[i, i - 1] = float(-D_mod_fast/dr**2 + D_mod_fast/(2 * dr * (dr * i +dr/2)))
    M[i, i + 1] = float(-D_mod_fast/dr**2 - D_mod_fast/(2 * dr * (dr * i + dr/2)))
M[m - 1, m - 2] = float(-D_mod_fast/dr**2 + D_mod_fast/(2 * dr * (dr * (m - 1) +dr/2)))
M[m - 1, m - 1] = float(D_mod_fast/dr**2 - D_mod_fast/(2 * dr * (dr * (m - 1) +dr/2)) + ZR_mod_fast)

for i in range(m + 1, m + m_f - 1):
    M[i, i] = float((2 * D_fuel_th)/dr**2 + ZR_fuel_th)
    M[i, i - 1] = float(-D_fuel_th/dr**2 + D_fuel_th/(2 * dr * (dr * (i - m) + dr/2)))
    M[i, i + 1] = float(-D_fuel_th/dr**2 - D_fuel_th/(2 * dr * (dr * (i - m) + dr/2)))
M[m, m] = float(D_fuel_th/dr**2 + D_fuel_th/(2 * dr * (dr/2)) + ZR_fuel_th)
M[m, m + 1] = float(-D_fuel_th/dr**2 - D_fuel_th/(2 * dr * (dr/2)))
M[m + 1, m] = float(-D_fuel_th/dr**2 + D_fuel_th/(2 * dr * (2*dr - dr/2)))

M[m + m_f - 1, m + m_f - 2] = float(-D_fuel_th/dr**2 + D_fuel_th/(2 * dr * (dr * (m_f - 1) + dr/2)))
M[m + m_f - 1, m + m_f - 1] = float(D_int_th/dr**2 + D_int_th/(2 * dr * (dr * (m_f - 1) + dr/2)) + D_fuel_th/dr**2 - D_fuel_th/(2 * dr * (dr * (m_f - 1) + dr/2)) + ZR_fuel_th)
M[m + m_f - 1, m + m_f] = float(-D_int_th/dr**2 - D_int_th/(2 * dr * (dr * (m_f - 1) + dr/2)))
M[m + m_f, m + m_f - 1] = float(-D_int_th/dr**2 + D_int_th/(2 * dr * (dr * m_f + dr/2)))
M[m + m_f, m + m_f] = float(D_mod_th/dr**2 + D_mod_th/(2 * dr * (dr * m_f + dr/2)) + D_int_th/dr**2 - D_int_th/(2 * dr * (dr * m_f + dr/2)) + ZR_mod_th)
M[m + m_f, m + m_f + 1] = float(-D_mod_th/dr**2 - D_mod_th/(2 * dr * (dr * m_f + dr/2)))
M[m + m_f + 1, m + m_f] = float(-D_mod_th/dr**2 + D_mod_th/(2 * dr * (dr * (m_f + 1) + dr/2)))

for i in range(m + m_f + 1, m + m - 1):
    M[i, i] = float(2 * D_mod_th/dr**2 + ZR_mod_th)
    M[i, i - 1] = float(-D_mod_th/dr**2 + D_mod_th/(2 * dr * (dr * (i - m) +dr/2)))
    M[i, i + 1] = float(-D_mod_th/dr**2 - D_mod_th/(2 * dr * (dr * (i - m) + dr/2)))

M[2*m - 1, 2*m - 2] = float(-D_mod_th/dr**2 + D_mod_th/(2 * dr * (dr * (m - 1) +dr/2)))
M[2*m - 1, 2*m - 1] = float(D_mod_th/dr**2 - D_mod_th/(2 * dr * (dr * (m - 1) +dr/2)) + ZR_mod_th)

for i in range(m, m + m_f):
    M[i, i - m] = float(-Zs_12_fuel)
for i in range(m + m_f, 2 * m):
    M[i, i - m] = float(-Zs_12_mod)

# inverse M Matrix
M_inv = np.linalg.inv(M)
# X Matrix creation
X = np.zeros((2 * m, 2 * m))

for i in range(0, m_f):
    X[i, i] = float(vZf_fuel_fast)
    X[i, i + m] = float(vZf_fuel_th)

# initial flux and eigenvalue guess
f0 = 10**7
k0 = 1

# initial F matrix creation
F0 = np.zeros((2 * m, 1))
for i in range(0, 2 * m):
    F0[i, 0] = float(f0)
ek = 1
es = 1

# radius matrix creation
R = np.zeros((m, 1))
for i in range(0, m):
    R[i, 0] = float(dr * i +dr/2)

while ek > e or np.all(es > e):
    # flux normalization
    N = Q / (Ef * 2 * np.pi * dr * (Zf_fast * np.sum(F0[0: m_f] * R[0: m_f]) + Zf_th * np.sum(F0[m: m + m_f] * R[0 : m_f])))
    Fn = F0 * N
    # initial source matrix creation
    S0 = np.dot(X, Fn)
    # F matrix creation
    F = np.dot(M_inv * (1 / k0), S0)
    # estimation of the source value
    S = np.dot(X, F)
    # estimation of k(n+1)
    k = k0 * np.sum(S) / np.sum(S0)
    # error estimation
    ek = np.abs((k - k0) / k)
    es = np.abs((np.max(S) - np.max(S0)) / np.max(S))
    F0 = F
    k0 = k

print('k = ', k)
np.savetxt('F.txt', F, fmt='%.5e')

Flux_fast = np.zeros((m, 1))
for i in range(0, m):
    Flux_fast[i, 0] = float(F[i])

Flux_th = np.zeros((m, 1))
for i in range(0, m):
    Flux_th[i, 0] = float(F[i + m])

# Plotting
x = np.linspace(dr / 2, r_tot - dr / 2, m)
plt.ylabel('Fast Neutron Flux, 1/(cm^2 * s)')
plt.xlabel('r, cm')
plt.grid(True)
plt.plot(x, Flux_fast)
plt.title('Fast neutron flux')
plt.show()

x = np.linspace(dr / 2, r_tot - dr / 2, m)
plt.ylabel('Thermal Neutron Flux, 1/(cm^2 * s)')
plt.xlabel('r, cm')
plt.grid(True)
plt.plot(x, Flux_th)
plt.title('Thermal neutron flux')
plt.show()

x = np.linspace(dr / 2, r_tot - dr / 2, 2* m)
plt.ylabel('Source, 1/(cm^3 * s)')
plt.xlabel('r, cm')
plt.grid(True)
plt.plot(x, (1/k) *S)
plt.title('Source')
plt.show()