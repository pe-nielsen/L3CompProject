#simResultDir = r'resultsNCC/dataCollection_lowTemp2'

simResultDir = r'resultsNCC/BWF_results'
#densities = [1e-3, 3e-3, 5e-3, 7e-3, 9e-3, 7.76e-1, 7.80e-1, 8.2e-1, 8.60e-1, 9e-1]
#redTemps = [0.85, 0.9]
#numPools = 20

densities = [0.6]
redTemps = [0.75, 1.5]
numPools = 2

#dataCollectionFullRange 0-1.5
#densities = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
#redTemps = [1, 2, 3, 5, 10]
#numPools = 30

#####
#densities = [0.1, 0.2, 0.3 ,0.4 ,0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
#redTemps = [1]
#redTemps = [2]
#redTemps = [3]
#redTemps = [5]
#redTemps = [10]

#densities = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.55, 0.575, 0.6, 0.625, 0.65, 0.7, 0.8, 0.9, 1.0]
#densities = [0.25, 0.35, 0.45, 1.1, 1.2, 1.3, 1.4, 1.5]
#redTemps = [0.75, 1.00, 1.25, 1.50]
#redTemps = [1.25, 1.50]
#numPools = 32



#####

#dataCollectionSmallDens 0-0.35
#densities = [0.05, 0.15, 0.25, 0.35, 0.45]
#redTemps = [1, 5, 10]
#numPools = 15

#dataCollectionMidDens 0.35-0.75
#densities = [0.55, 0.65, 0.675, 0.725, 0.75]
#redTemps = [1, 5, 10]
#numPools = 15

#dataCollectionHighDens 0.75-1
#densities = [0.76, 0.77, 0.78, 0.79, 0.825, 0.85, 0.875, 0.925, 0.95, 0.975]
#redTemps = [1, 5, 10]
#numPools = 15


# numPools = 20

numParts = 10**3
dim = 3

# Lennard-Jones Parameters:
sigma = 3.4e-10
# sigma = 3.345e-10  # +/- 0.04e-10, Chem. Phys. 111, 9352 (1999); https://doi.org/10.1063/1.479848
well_depth = 1.65e-21
Boltzmann_const = 1.38064852e-23  # m^2.kg.s^-2.K-1

numEquilMoves = int(2e7)
numEvolveMoves = int(2e7)

numEquilIter = 1000
initDisplacement = 0.1
configSampleRate = 2000
