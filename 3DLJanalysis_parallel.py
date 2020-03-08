from os import listdir
from os.path import isfile, join
from multiprocessing import Pool

import matplotlib.pyplot as plt

from resultClasses import simResult, compResult
from commonFunctions import getDistanceSquared, redLJ, diff_redLJ

from time import time
from datetime import timedelta
from pathlib import Path
import pickle

import numpy as np
np.random.seed(0)

plt.style.use('seaborn-deep')
plt.style.use(r'C:\Users\splb68\comp_proj\scripts\PaperDoubleFig.mplstyle')


def computeQuantities(sR):
    binWidthFactor = 0.1
    dr = binWidthFactor * sR.partRad
    rMax = 10 * sR.partRad

    edges = np.arange(0.0, rMax + 1.1*dr, dr)
    numIncrements = len(edges) - 1
    numSamples = sR.configs.shape[2]
    nCount = np.zeros(shape=(sR.numParts, numIncrements, numSamples))  # counted particles
    radii = np.zeros(numIncrements)
    numberDensity = sR.numParts / (sR.containerLength**3)

    potnEnergies = np.zeros(shape=(sR.numParts, sR.configs.shape[2]))  # potnEnergies for each [particle, sample]

    for p in range(sR.numParts):
        if p % 250 == 0:
            print(f'p: {p/1000}')
        distsSq = getDistanceSquared(sR.configs[p, :, :], sR.configs, sR.containerLength, dim=3)  # for []
        potnEnergies[p, :] = np.nansum(redLJ(distsSq))

        distsSq[p, :] = (1000*rMax)**2

        for s in range(distsSq.shape[1]):  # for each sample compute as normal
            distances = np.sqrt(distsSq[:, s])
            hist, bins = np.histogram(distances, bins=edges, density=False)
            nCount[p, :, s] = hist

    avPCF = np.zeros(shape=(numIncrements, numSamples))
    for i in range(numIncrements):
        radii[i] = 0.5*(edges[i]+edges[i+1])
        nIdeal = numberDensity * 4*np.pi * (edges[i] - dr/2)**2 * dr
        avPCF[i, :] = np.mean(nCount[:, i, :], axis=0)/nIdeal

    totAvPCF = np.mean(avPCF[:, :], axis=1)
    totAvPCF_stdev = np.std(avPCF[:, :], axis=1)

    cR = compResult(sR.numParts, sR.sigma, sR.wellDepth, sR.redTemp, sR.numEvolveMoves, sR.displacement, sR.mode, sR.density,
                 sR.partRad, sR.containerLength, totAvPCF, totAvPCF_stdev, radii)

    # calculating the energy and pressure
    LJpotns = redLJ(cR.PCFradii ** 2)  # the LJ potential at the PCF radii
    diff_LJpotns = diff_redLJ(cR.PCFradii)  # du/dr at the PCF radii

    intEn_is = LJpotns * cR.PCF * (cR.PCFradii ** 2) * dr  # array of internal energy sum terms at PCF radii
    pressure_is = diff_LJpotns * cR.PCF * (cR.PCFradii ** 3) * dr  # array of pressure terms at PCF radii

    intEn_sum = intEn_is.sum()  # sum of all the internal energy terms
    pressure_sum = pressure_is.sum()  # sum of all the pressure terms

    intEn_perPart = 2 * np.pi * numberDensity * intEn_sum
    pressure_minRhokT = (2 / 3) * np.pi * numberDensity ** 2 * pressure_sum

    cR.intEn_perPart = intEn_perPart
    cR.pressure_minRhokT = pressure_minRhokT
    print(f'\nfor density = {sR.density},\n'
          f'internal energy per particle: {intEn_perPart},\n'
          f'P - rho.kT = {pressure_minRhokT}\n')

    return cR


def main(sR):
    """calculate and save quantities for a single simResult"""
    cR = computeQuantities(sR)

    #saving the compResult:

    cR_sub_dir = f'compResult\\redTemp{cR.redTemp}'
    cR_parent_dir = join(parent_dir, cR_sub_dir)
    Path(cR_parent_dir).mkdir(parents=True, exist_ok=True)
    cr_data_path = cR_parent_dir + '\\' + f'den{int(cR.density * 100)}'

    outfile = open(cr_data_path, 'wb')
    toDump = cR
    pickle.dump(toDump, outfile)
    outfile.close()

    print(f'''\n\n
            Results computed saved for:\n
            density: {sR.density},\n
            reduced Temperature: {sR.redTemp}
            At time ---{timedelta(seconds=(time() - startTime))}--- from start
            \n\n''')


if __name__ == '__main__':
    startTime = time()
    print(f'Start Time: ---{timedelta(seconds=startTime)}---')

    # finding all simResult files
    simResults = []  # holds list of all simResults

    parent_dir = r'C:\Users\splb68\comp_proj\simulationResults\smallRun_isotherms'
    redTemp_dirs = [f for f in listdir(parent_dir) if not isfile(join(parent_dir, f)) and 'redTemp' in f]

    for i in range(len(redTemp_dirs)):
        rT_dir = redTemp_dirs[i]
        dps = join(parent_dir, rT_dir)
        data_files = [f for f in listdir(dps) if
                      isfile(join(dps, f)) and ('den' in f)]
        for j in range(len(data_files)):
            df = data_files[j]
            data_path = join(dps, df)
            infile = open(data_path, 'rb')
            sR = pickle.load(infile)
            infile.close()
            simResults.append(sR)

    print(f'simResults: {simResults}')

    numPools = 10  # number of cores used for calculations
    with Pool(numPools) as p:
        p.map(main, simResults)

    endTime = time() - startTime
    print(f'Run Time:   ---{timedelta(seconds=endTime)}---')
    exit()
