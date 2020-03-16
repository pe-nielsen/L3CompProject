simResultDir = r'resultsNCC/NIST_params_3'

densities = [1e-3, 3e-3, 5e-3, 7e-3, 9e-3, 7.76e-1, 7.80e-1, 8.2e-1, 8.60e-1, 9e-1]
redTemps = [0.85, 0.9]

numPools = 20

numParts = 10**3
dim = 3

# Lennard-Jones Parameters:
sigma = 3.4e-10
# sigma = 3.345e-10  # +/- 0.04e-10, Chem. Phys. 111, 9352 (1999); https://doi.org/10.1063/1.479848
well_depth = 1.65e-21
Boltzmann_const = 1.38064852e-23  # m^2.kg.s^-2.K-1

numEquilMoves = int(1e7)
numEvolveMoves = int(1.5e7)

numEquilIter = 1000
initDisplacement = 0.1
configSampleRate = 2000
