simResultDir = r'resultsNCC/shortRun_badConvergence'
# simResultDir = r'simulationResults/debugRun2'

# densities = [3, 4, 5, 5.25, 5.5, 5.75, 6]
# redTemps = [1, 11]
densities = [1, 2, 3, 3.5, 4, 4.5, 5, 5.2, 5.4, 5.6, 5.8, 6, 6.2]
redTemps = [1, 10]
# densities=[4.5, 5, 5.5, 6]
# redTemps=[5]

numPools = 26

numParts = 10**3
dim = 3

# Lennard-Jones Parameters:
sigma = 3.4e-10
# sigma = 3.345e-10  # +/- 0.04e-10, Chem. Phys. 111, 9352 (1999); https://doi.org/10.1063/1.479848
well_depth = 1.65e-21
Boltzmann_const = 1.38064852e-23  # m^2.kg.s^-2.K-1

numEquilMoves = int(5e6)
numEvolveMoves = int(1e7)


numEquilIter = 2000
initDisplacement = 0.1
configSampleRate = 1000
