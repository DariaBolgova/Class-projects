import numpy as np

Ntry = 1000000 #number of particles
ST = 0.1 #total macroscopic cross-section value in cm^(-1)
SS = 0.075 #scattering macroscopic cross-section value in cm^(-1)
SA = 0.025 #absorbtion macroscopic cross-section value in cm^(-1)
PA = SA/ST #probability of absorption
PS = SS/ST #probability of scattering
Xmax = 5 #right border of the slab, cm
Xmin = - 5 #left border of the slab, cm
NL = 0 #leakage count
NA = 0 #adsorption count
for i in range (0,Ntry): #counting limit
    XC = 0  #start x coordinate
    M = 2*np.random.uniform(0, 1) - 1 #random choice of mu
    for j in range (0, Ntry):
        #determination of the distance to the border depending on the angle
        if M > 0:
            db = (Xmax - XC)/M
        else:
            db = (Xmin - XC) / M

        x = (-1) * np.log(np.random.uniform(0, 1)) / ST  # determining the length of the flight path to collision
        if x>db: #if the particle is outside the slab
            NL = NL + 1 #count +1 to leaks
            break
        else: #if the particle is inside the slab
            XC = XC + x*M #defining a new starting coordinate

            if np.random.uniform(0, 1) < PA: #determination of the reaction type, if it is absorbtion
                NA = NA +1 #count +1 to absorbtions
                break
            else: #if it is scattering
                M = 2*np.random.uniform(0, 1) - 1 #defining a new angle


Leak = NL/Ntry #leakage fraction
Absorb = NA/Ntry #absorbtion fraction
print('The absorption fraction is', Absorb,'. The leaks fraction is', Leak, '.')

