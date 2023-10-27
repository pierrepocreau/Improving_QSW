from game import Game
from hierarchy import Hierarchy
import cvxpy as cp
import numpy as np

# Each player has 3 operators: 0 corresponds to the identity operator,
# the others correspond to two different basis of measurement.
operatorsP1 = [0, 1, 2]
operatorsP2 = [0, 3, 4]
operatorsP3 = [0, 5, 6]
operatorsP4 = [0, 7, 8]
operatorsP5 = [0, 9, 10]

nbPlayers = 3 # Change to 5 for the C5 games
P3 = [operatorsP1, operatorsP2, operatorsP3]
P5 = [operatorsP1, operatorsP2, operatorsP3, operatorsP4, operatorsP5]

v0 = cp.Parameter()
v1 = cp.Parameter()

game = Game(nbPlayers, v0, v1, sym=False) # Set sym=True for NC_01(C5)
prob = Hierarchy(game, P3, level=1) # Change to P5 for C5 games
prob.setNashEqConstraints()

# Compute the NPA hierarchy bound
qswList = []
for k in np.arange(0, 2.01, 0.01):
    v0.value = k
    v1.value = 2-k
    qsw = prob.optimize(verbose=True, warmStart=False, solver="MOSEK")
    qswList.append(qsw)

print(qswList)