import numpy as np
import time

Ntry = 10000000 #number of throws
Na = 0 #number of hits
for i in range (0,Ntry-1): #counting limit
    x = np.random.uniform(0,1) #random x
    y = np.random.uniform(0,1)  # random y
    S = float(x**2 + y**2)
    if S <= 1: #hit condition
        Na = Na +1 #counter
Pi = 4*Na/Ntry #Pi calculation
end = time.time() #turning off the clock
T = end - start #time calculation
print(Pi)
print(T)