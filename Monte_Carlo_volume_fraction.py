import numpy as np

Ntry = 10000000 #number of throws
Na = 0 #number of hits
for i in range (0,Ntry-1): #counting limit
    x = np.random.uniform(-2,2) #random x
    y = np.random.uniform(-2,2)  # random y
    z = np.random.uniform(-2,2)  # random z
    S1 = ((x+0.25)**2 + y**2 + z**2) #equation of sphere 1
    S2 = ((x-0.25)**2 + y**2 + z**2) #equation of sphere 2
    if S1 <= 1 or S2 <= 1: #hit condition
        Na = Na +1 #counter
F = Na/Ntry #volume fraction
print('The fraction of the volume occupied by the spheres is', F)
