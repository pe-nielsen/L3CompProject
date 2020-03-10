simResultDir = r'simulationResults/high_densities'

# densities = [3, 4, 5, 5.25, 5.5, 5.75, 6]
# redTemps = [1, 11]
densities = [3, 3.5, 4, 4.5, 5, 5.25, 5.5, 5.75, 6, 6.25]
redTemps = [5]

numPools = 10

numParts = 10**3
dim = 3

# Lennard-Jones Parameters:
sigma = 3.4e-10
well_depth = 1.65e-21
Boltzmann_const = 1.38064852e-23  # m^2.kg.s^-2.K-1

numEquilMoves = int(1e7)
numEvolveMoves = int(2e7)

numEquilIter = 1000
initDisplacement = 0.1
configSampleRate = 1000
