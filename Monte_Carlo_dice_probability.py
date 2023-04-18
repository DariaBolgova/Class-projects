import numpy as np
import random
import matplotlib.pyplot as plt

Ntry = 1000000000 #the maximum number of trials
Nbin = 11 #number of possible vasues
tally = np.zeros(Nbin) #empty cells for sorting dropped values
T = 0 #dice roll counter
D = [1, 2, 3, 4, 5, 6] #possible values on a dice
n = 3 #required precision of the probability distribution, that is, the number of decimal places
for i in range (1,Ntry): #counting limit
    R1 = random.choice(D)  #dice roll 1
    R2 = random.choice(D)  #dice roll 2
    R = R1 + R2 #the sum of the values rolled on each dice
    bin = R - 1 #selection of cell for sorting
    tally[bin-1] = tally[bin-1] + 1 #record +1 in the desired cell
    T = T + 1 #record +1 to the dice rolls
    P = tally*100/ T  #probability calculation
    if (i % 100000) == 0: #result output every 100.000 rolls
        print(P)
    if abs(P[0]-P[10])<10**(-n) and abs(P[1]-P[9])<10**(-n) and abs(P[2]-P[8])<10**(-n) and abs(P[3]-P[7])<10**(-n) and abs(P[4]-P[6])<10**(-n) and tally[5]>10000:
        break #condition for stopping dice rolls, provided that the probability distribution is symmetric up to n decimal places (but not less than 10,000 rolls)

print(P) #the probability distribution
print(T) #number of rolls

fig, ax = plt.subplots()
g = np.linspace(2, 12, 11)
ax.plot(g, P)
plt.title('Probability of getting different values on two dice')
plt.xlabel("The sum of the values on two dice")
plt.ylabel("The probability, %")
plt.grid(True)
plt.show()