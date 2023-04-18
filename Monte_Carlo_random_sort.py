import numpy as np
import math
import matplotlib.pyplot as plt

Ntry = 10000000 #amount of random numbers
Nbin = 10 #number of cells
tally = np.zeros(Nbin) #empty cells for sorting random numbers
for i in range (0,Ntry-1): #counting limit
    R = np.random.uniform(0,1) #random number creation
    bin = 1 + math.floor(Nbin*R) #selection of cell for sorting
    tally[bin-1] = tally[bin-1] + 1 #record +1 in the desired cell
    if (i % 1000000) == 0: #result output every 1M for tally
        print(tally)
freq = tally / Ntry #frequency calculation
print(freq) #result output for frequency

plt.figure(1)
plt.plot(tally)
plt.grid(True)
plt.show()

plt.figure(2)
plt.plot(freq)
plt.grid(True)
plt.show()
