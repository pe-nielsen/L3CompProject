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
import bz2


import numpy as np
import numexpr as ne
np.random.seed(0)

plt.style.use('seaborn-deep')
plt.style.use(r'PaperDoubleFig.mplstyle')

def computeQuantities(sRbWF):
    sR = sRbWF[0]
    binWidthFactor = sRbWF[1]
    # binWidthFactor = 0.1
    dr = binWidthFactor * sR.partRad
    rMax = 6 * sR.partRad

    edges = np.arange(0.0, rMax + 1.1*dr, dr)
    numIncrements = len(edges) - 1
    numSamples = sR.configs.shape[2]
    # nCount = np.zeros(shape=(sR.numParts, numIncrements, numSamples))  # counted particles
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

            if p == 0:
                avnCount[:, s] = hist
            else:
                avnCount[:, s] = (avnCount[:, s] + hist)/2

    # avPCF = np.zeros(shape=(numIncrements, numSamples))
    for i in range(numIncrements):
        radii[i] = 0.5 * (edges[i] + edges[i + 1])
        nIdeal = numberDensity * 4 * np.pi * (edges[i] - dr / 2) ** 2 * dr
        # avPCF[i, :] = np.mean(nCount[:, i, :], axis=0)/nIdeal
        avPCF[i, :] = avnCount[i, :] / nIdeal

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

    cR.intEn_perPart = intEn_perPart
    cR.pressure_minRhokT = pressure_minRhokT
    cR.intEnstdev = stdev_intEn
    cR.pressurestdev = stdev_pressure
    print(f'\nfor density = {sR.density},\n'
          f'internal energy per particle: {intEn_perPart},\n'
          f'P - rho.kT = {pressure_minRhokT}\n')

    return cR


def main(sRbWF):
    """calculate and save quantities for a single simResult"""
    sR = sRbWF[0]
    binWidthFactor = sRbWF[1]
    cR = computeQuantities(sRbWF)

    #saving the compResult:

    parent_dir = r'resultsNCC/BWF_results'
    cR_sub_dir = f'compResultLargedr/redTemp{cR.redTemp}'
    cR_parent_dir = join(parent_dir, cR_sub_dir)
    Path(cR_parent_dir).mkdir(parents=True, exist_ok=True)
    cR_data_path = join(cR_parent_dir, f'den{int(cR.density)}bWF{int(binWidthFactor)}')


    #
    # cR_sub_dir = f'compResults'
    # parent_dir = r'report/figures/binWidth'
    # cR_parent_dir = join(parent_dir, cR_sub_dir)
    # Path(cR_parent_dir).mkdir(parents=True, exist_ok=True)
    # # cr_data_path = cR_parent_dir + '\\' + f'den{int(cR.density * 100)}bWF{int(binWidthFactor*10000)}'
    # cr_data_path = join(cR_parent_dir, f'den{int(cR.density)}bWF{int(binWidthFactor)}')
    # print(f'den{int(cR.density)}bWF{int(binWidthFactor)}')

    outfile = bz2.BZ2File(cR_data_path, 'w')
    toDump = (binWidthFactor, cR)
    pickle.dump(toDump, outfile)
    outfile.close()

    print(f"""\n\n
            Results computed saved for:
            density: {sR.density},
            bin width factor: {binWidthFactor},
            reduced Temperature: {sR.redTemp}
            At time ---{timedelta(seconds=(time()))}---
            \n\n""")


if __name__ == '__main__':
    startTime = time()
    print(f'Start Time: ---{timedelta(seconds=startTime)}---')

    # finding all simResult files
    simResults = []  # holds list of all simResults

    parent_dir = r'resultsNCC/BWF_results'
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

    # df = 'rT1.5den0.6'
    # data_path = join(parent_dir, df)
    # infile = bz2.BZ2File(data_path, 'r')
    # sR = pickle.load(infile)
    # infile.close()

    print(f'simResults: {simResults}')

    # binWidthFactors = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15]
    # binWidthFactors = [0.005, 0.01, 0.015, 0.02, 0.025, 0.04, 0.045, 0.05]
    binWidthFactors = [0.001, 0.002, 0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    print(f'binwidthfactors: {binWidthFactors}')
    sRbWFs = []
    for sR in simResults:
        for bWF in binWidthFactors:
            sRbWF = [sR, bWF]
            sRbWFs.append(sRbWF)
            print(sR)
            print(bWF)

    numProccesses = 28  # number of cores used for calculations
    with Pool(processes=numProccesses) as p:
        p.map(main, sRbWFs)

    endTime = time() - startTime
    print(f'Run Time:   ---{timedelta(seconds=endTime)}---')
    exit()
