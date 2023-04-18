import numpy as np
import math
import matplotlib.pyplot as plt

Ntry = 1000000 #number of trials
Nbin = 50 #number of bins
S = 1 #macroscopic cross-section value
dx = 1/10 #delta x value
tally = np.zeros(Nbin) #creation of bins for sorting
f = np.zeros(Nbin) #creation of frequency bins
for i in range (0,Ntry-1): #counting limit
    R = np.random.uniform(0,1) #random choice of xye
    x = (-1)*np.log(R)/S #corresponding x determination
    if x <= 5: #x value limit
        bin = 1 + math.floor(x / dx)  #selection of proper bin
        tally[bin-1] = tally[bin-1] + 1  #record +1 in the corresponding bin
for i in range (0,Nbin-1): #in range of bins number
    f[i] = 1/Ntry * tally[i]/dx #frequency values for i bin

fig, ax = plt.subplots()
g = np.linspace(dx/2, (Nbin-1)*dx+dx/2, 50)
y = np.exp(-g) #decreasing the argument by dx times to alogn the scales
ax.plot(g, y)
ax.plot(g, f)
plt.title('Random sample')
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.show()

