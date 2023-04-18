import numpy as np
import math
import matplotlib.pyplot as plt

Ntry = 1000000 #number of particles
S = 1 #macroscopic cross-section value in cm^(-1)
db = 2 #distance to the boundary in cm
Ntr = 0 #counter of particles transmitted through the slab
for i in range (0,Ntry-1): #counting limit
    R = np.random.uniform(0,1) #random choice of xye
    M = np.random.uniform(0,1) #random choice of mu
    x = (-1)*np.log(R)/S #corresponding x determination
    if x > db/M: #x exceeds the distance to the boundary
        Ntr = Ntr + 1  # counter
P = Ntr*100/Ntry
print('The percent transmission through the slab is', P)

