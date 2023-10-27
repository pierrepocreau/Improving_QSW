Scripts of the Matlab program:
- **runSeeSaw.m**: run the seesaw algorithm to compute a lower bound on the quantum social welfare, while veryfing either some Nash or quantum equilibrium constraints.
- **loadNPAData.m**: load the data from the NPA hierarchy contained in the data folder and computed with the Python program of the NPA folder.

Functions of the Matlab program:
- **isQuantumEquilibrium**: verify that a given solution (state and measurement) is a quantum equilbrium for a game.
- **isQuantumCorrelatedEquilibrium.m**: verify that a given solution (state and measurement) is a Nash equilbrium for a game.

The folders data and utils should be added to the Matlab path for the code to function.