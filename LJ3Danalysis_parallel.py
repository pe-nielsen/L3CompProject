from os import listdir
from os.path import isfile, join
from multiprocessing import Pool

import matplotlib.pyplot as plt

from resultClasses import simResult, compResult
from commonFunctions import getDistanceSquared, redLJ, diff_redLJ
from initial_parameters import simResultDir

from time import time
from datetime import timedelta
from pathlib import Path
import pickle
import bz2

import numpy as np
import numexpr as ne
# np.random.seed(0)

plt.style.use('seaborn-deep')
plt.style.use(r'PaperDoubleFig.mplstyle')


def computeQuantities(sR):
    binWidthFactor = 0.02
    # binWidthFactor = 0.1
    dr = binWidthFactor * sR.partRad
    rMax = 6 * sR.partRad

    edges = np.arange(0.0, rMax + 1.1*dr, dr)
    numIncrements = len(edges) - 1
    numSamples = sR.configs.shape[2]
    # nCount = np.zeros(shape=(sR.numParts, numIncrements, numSamples), dtype=np.float32)  # counted particles
    radii = np.zeros(numIncrements)
    numberDensity = sR.numParts / (sR.containerLength**3)

    avnCount = np.zeros(shape=(numIncrements, numSamples))
    avPCF = np.zeros(shape=(numIncrements, numSamples))

    # potnEnergies = np.zeros(shape=(sR.numParts, sR.configs.shape[2]))  # potnEnergies for each [particle, sample]

    for p in range(sR.numParts):
        if p % 250 == 0:
            print(f'p: {p/1000}')
        distsSq = getDistanceSquared(sR.configs[p, :, :], sR.configs, sR.containerLength, dim=3)  # for []
        # potnEnergies[p, :] = np.nansum(redLJ(distsSq))

        distsSq[p, :] = (1000*rMax)**2

        for s in range(distsSq.shape[1]):  # for each sample compute as normal
            distSq_i = distsSq[:, s]
            # distances = np.sqrt(distsSq[:, s])
            distances = ne.evaluate("sqrt(distSq_i)")
            hist, bins = np.histogram(distances, bins=edges, density=False)
            # nCount[p, :, s] = hist

            if p == 0:
                avnCount[:, s] = hist
            else:
                avnCount[:, s] = (avnCount[:, s] + hist)/2



    for i in range(numIncrements):
        radii[i] = 0.5*(edges[i]+edges[i+1])
        nIdeal = numberDensity * 4*np.pi * (edges[i] - dr/2)**2 * dr
        # avPCF[i, :] = np.mean(nCount[:, i, :], axis=0)/nIdeal
        avPCF[i, :] = avnCount[i, :]/nIdeal


    # for tracking convergence:
    indexAtSigma = int(np.ceil(1/binWidthFactor))
    gAtSigma = avPCF[indexAtSigma, :]

    # find a way to save gAtSigma with all of the relevant information
    # plot gAtSigma against first index ~ and relate first index to the number of MC steps taken





    totAvPCF = np.mean(avPCF[:, :], axis=1)
    totAvPCF_stdev = np.std(avPCF[:, :], axis=1)

    cR = compResult(sR.numParts, sR.sigma, sR.wellDepth, sR.redTemp, sR.numEvolveMoves, sR.displacement, sR.mode, sR.density,
                 sR.partRad, sR.containerLength, totAvPCF, totAvPCF_stdev, radii)

    # calculating the energy and pressure
    LJpotns = redLJ(cR.PCFradii ** 2)  # the LJ potential at the PCF radii
    diff_LJpotns = diff_redLJ(cR.PCFradii)  # du/dr at the PCF radii

    PCF = cR.PCF
    PCFradii = cR.PCFradii
    PCFERRORin = cR.PCFstdev

    intEn_is = ne.evaluate("LJpotns * PCF * (PCFradii ** 2) * dr")  # array of internal energy sum terms at PCF radii
    pressure_is = ne.evaluate("diff_LJpotns * PCF * (PCFradii ** 3) * dr")  # array of pressure terms at PCF radii

    intEn_sum = intEn_is.sum()  # sum of all the internal energy terms
    pressure_sum = pressure_is.sum()  # sum of all the pressure terms

    pi = np.pi
    intEn_perPart = ne.evaluate("2 * pi * numberDensity * intEn_sum")
    pressure_minRhokT = ne.evaluate("-(2/3) * pi * (numberDensity**2) * pressure_sum")

    # error from standard deviation
    stdev_intEn = 2*np.pi*numberDensity*np.std(intEn_is)
    stdev_pressure = (2/3) * pi * (numberDensity**2) * np.std(pressure_is)


    # intEn_isERRORinSq = ne.evaluate("(LJpotns * PCFERRORin * (PCFradii ** 2) * dr)**2")
    # pressure_isERRORinSq = ne.evaluate("(diff_LJpotns * PCFERRORin * (PCFradii ** 3) * dr)**2")

    # intEn_sumERRORin = intEn_isERRORinSq.sum()
    # pressure_sumERRORin = pressure_isERRORinSq.sum()

    # intEn_perPartERRORin = ne.evaluate("2 * pi * numberDensity * intEn_sumERRORin")
    # pressure_minRhokTERRORin = ne.evaluate("-(2/3) * pi * (numberDensity**2) * pressure_sum")

    cR.intEn_perPart = intEn_perPart
    cR.pressure_minRhokT = pressure_minRhokT
    cR.intEnstdev = stdev_intEn
    cR.pressurestdev = stdev_pressure
    print(f'\nfor density = {sR.density},\n'
          f'internal energy per particle: {intEn_perPart},\n'
          f'P - rho.kT = {pressure_minRhokT}\n')

    return cR, gAtSigma


def runComp(sR):
    """calculate and save quantities for a single simResult"""
    cR, gAtSigma = computeQuantities(sR)

    #saving the compResult:
    # parent_dir = r'simulationResults/smallRunIsotherms3'
    parent_dir = simResultDir
    cR_sub_dir = f'compResultLargedr/redTemp{cR.redTemp}'
    cR_parent_dir = join(parent_dir, cR_sub_dir)
    Path(cR_parent_dir).mkdir(parents=True, exist_ok=True)
    cR_data_path = join(cR_parent_dir, f'den{int(cR.density * 100)}')

    outfile = bz2.BZ2File(cR_data_path, 'w')
    toDump = cR
    pickle.dump(toDump, outfile)
    outfile.close()

    print(f'''\n\n
            Results computed saved for:\n
            density: {sR.density},\n
            reduced Temperature: {sR.redTemp}
            At time ---{timedelta(seconds=(time()))}---
            \n\n''')


if __name__ == '__main__':
    startTime = time()
    print(f'Start Time: ---{timedelta(seconds=startTime)}---')

    # finding all simResult files
    simResults = []  # holds list of all simResults

    parent_dir = simResultDir
    redTemp_dirs = [f for f in listdir(parent_dir) if not isfile(join(parent_dir, f)) and 'redTemp' in f]

    for i in range(len(redTemp_dirs)):
        rT_dir = redTemp_dirs[i]
        dps = join(parent_dir, rT_dir)
        data_files = [f for f in listdir(dps) if
                      isfile(join(dps, f)) and ('den' in f)]
        for j in range(len(data_files)):
            df = data_files[j]
            data_path = join(dps, df)
            infile = bz2.BZ2File(data_path, 'r')
            sR = pickle.load(infile)
            # infile = open(data_path, 'rb')
            # sR = pickle.load(infile)
            infile.close()
            simResults.append(sR)

    print(f'simResults: {simResults}')

    numPools = 6  # number of cores used for calculations
    with Pool(numPools) as p:
        p.map(runComp, simResults)

    endTime = time() - startTime
    print(f'Run Time:   ---{timedelta(seconds=endTime)}---')
    exit()
