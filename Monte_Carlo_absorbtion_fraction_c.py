import numpy as np

Ntry = 1000000 #number of particles
ST = 0.1 #total macroscopic cross-section value in cm^(-1)
SS = 0.075 #scattering macroscopic cross-section value in cm^(-1)
SA = 0.025 #absorbtion macroscopic cross-section value in cm^(-1)
PA = SA/ST #probability of absorption
PS = SS/ST #probability of scattering
Xmax = 5 #left border of the slab, cm
Xmin = - 5 #right border of the slab, cm
NL1 = 0 #problem tallies for scores
NA1 = 0 #roblem tallies for scores^2
NL2 = 0 #problem tallies for scores
NA2 = 0 #roblem tallies for scores^2
wgt = 1 # particle weight
for i in range (0,Ntry): #counting limit
    XC = 0  #start x coordinate
    M = 2*np.random.uniform(0, 1) - 1 #random choice of mu
    NL = 0  # leakage count
    NA = 0  # adsorption count

    for j in range (0, Ntry):
        #determination of the distance to the border depending on the angle
        if M > 0:
            db = (Xmax - XC)/M
        else:
            db = (Xmin - XC) / M

        x = (-1) * np.log(np.random.uniform(0, 1)) / ST  # determining the length of the flight path to collision
        if x>db: #if the particle is outside the slab
            NL = NL + wgt #count +wgt to leaks
            break
        else: #if the particle is inside the slab
            XC = XC + x*M #defining a new starting coordinate

            if np.random.uniform(0, 1) < PA: #determination of the reaction type, if it is absorbtion
                NA = NA + wgt #count +wht to absorbtions
                break
            else: #if it is scattering
                M = 2*np.random.uniform(0, 1) - 1 #defining a new angle

    NA1 = NA1 + NA  # number of absorbtion
    NA2 = NA2 + NA ** 2
    NL1 = NL1 + NL  # number of leaks
    NL2 = NL2 + NL ** 2

Leak_ave = NL1 / Ntry  # fraction of leaks
Leak_std = ((NL2 / (Ntry - 1) - Leak_ave ** 2) / Ntry) ** (1 / 2)  # standard deviation of leaks
Abs_ave = NA1 / Ntry  # fraction of absorbtion
Abs_std = ((NA2 / (Ntry - 1) - Abs_ave ** 2) / Ntry) ** (1 / 2)  # standard deviation of absorbtion

print('The absorption fraction is', Abs_ave,'+/-', Abs_std, '. The leaks fraction is', Leak_ave, '+/-', Leak_std)


