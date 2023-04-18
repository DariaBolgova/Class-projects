
import random

print('Monty Hall presents you with three doors.')
print('Behind one of the doors is a valuable prize, while the other two doors conceal goats.')
print('You need to choose a door, but before it is opened, Monty, who knows what is behind each door, opens one of the other two doors to reveal a goat.')
print('You are then given the choice to stick with your initial choice or switch to the other unopened door.')
print('The question is: should you stick with the initial choice or switch to the other door to maximize your chances of winning the prize?')
Ng = 10000000 #number of games
S = input("What strategy will you use? Enter 1 if always switch and enter 2 if never switch: ") #you should choose a strategy
Nw = 0 #counter of wins before start
for i in range (0,Ng-1): #the program will play Ng games in turn
    doors = [1, 2, 3]  # array of doors
    R = [1, 2, 3]  # copy of the array of doors that is used for Monty's choice
    V = [1, 2, 3]  # copy of the array of doors that is used for the final choice
    P = random.choice(doors)  # random selection of the prize door
    C = random.choice(doors)  # random selection of the your door
    if P == C:  # if a prize door was chosen
        R.remove(P)  # then the prize door is removed from the door list
        M = random.choice(R)  # Monty chooses from the two remaining doors
    else:  # if a prize door was not chosen
        R.remove(C)  # then the prize door is removed from the door list
        R.remove(P)  # and your door is also removed from the list
        M = random.choice(R)  # Monty chooses the remaining door
    if S == '1': #if the strategy is always switch
        V.remove(C)  #then the initial choice is removed from the door list
        V.remove(M)  #and the Monty's choice is removed from the door list as well
        A = random.choice(V)  #so the final choice is the remaining door
    else: #if the strategy is never switch
        A = C #then the final choice is the same as the initial choice
    if A == P: #if the final choice is the prize door
        Nw = Nw + 1 #then +1 is added to the wins

W = Nw*100/Ng #calculation of the percentage of wins
print('You won', Nw, 'games or', W, '% from', Ng, 'games')







