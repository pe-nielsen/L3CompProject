from resultClasses import *
from commonFunctions import *
from initial_parameters import *
from LJ3D_parallel import Container, genMoves_rndV, initialiseSystem, equilibrateSystem, evolveSystem, runSim
from LJ3Danalysis_parallel import computeQuantities, runComp
from plotCompResults import importCompResults, plotCompResults

from os import listdir
from os.path import isfile, join
from pathlib import Path
from multiprocessing import Pool

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from time import time
from datetime import timedelta
import pickle
import bz2

import numpy as np
import numexpr as ne

# np.random.seed(0)
plt.style.use('seaborn-deep')
plt.style.use(r'PaperDoubleFig.mplstyle')


def model(tempDens):
    redTemp = tempDens[0]
    density = tempDens[1]
    displacement = initDisplacement
    mode = 'LJ'

    print(f"""
    current density: {density}
    current redTemp: {redTemp}
    """)

    container = initialiseSystem(numParts, density, dim, sigma, redTemp, mode)
    configEquil = equilibrateSystem(container, numEquilMoves, numEquilIter, displacement)

    # saving the equilibrated configuration
    compResultDir = join(r'compResult', f'redTemp{redTemp}')
    fileDir = join(simResultDir, compResultDir)
    Path(fileDir).mkdir(parents=True, exist_ok=True)

    filePathEquilConfig = join(fileDir, f'equilConfig_T{redTemp}_d{density}')
    outfile = bz2.BZ2File(filePathEquilConfig, 'w')
    toDump = configEquil
    pickle.dump(toDump, outfile)
    outfile.close()

    configEquilDEPRICATED, configs = evolveSystem(container, numEvolveMoves, configSampleRate)

    sR = simResult(container.numParts, container.sigma, well_depth, container.redTemp,
                   numEvolveMoves, container.displacement, container.mode, container.density,
                   container.partRad, container.length, configs)


    # ####SAVE SINGLE SR FOR CALCULATING BIN WIDTHS####
    # filepath1 = join(fileDir, f'simResult')
    # Path(filepath1).mkdir(parents=True, exist_ok=True)
    sR_dir = join(simResultDir, f'redTemp{sR.redTemp}')
    Path(sR_dir).mkdir(parents=True, exist_ok=True)
    filepath = join(sR_dir, f'den{sR.density}')
    outfile = bz2.BZ2File(filepath, 'w')
    toDump = sR
    pickle.dump(toDump, outfile)
    outfile.close
    # #### REMOVE THIS ASAP ####




    print(f"""
    Evolution completed for:
    density: {sR.density},
    redTemp: {sR.redTemp},
    at time: ---{timedelta(seconds=(time()))}---''',
    """)

    cR = computeQuantities(sR)

    # saving computation results
    filePath = join(fileDir, f'den{cR.density}')
    outfile = bz2.BZ2File(filePath, 'w')
    toDump = cR
    pickle.dump(toDump, outfile)
    outfile.close()



    print(f"""
        Computation completed and saved for:
        density: {sR.density},
        redTemp: {sR.redTemp},
        at time: ---{timedelta(seconds=(time()))}---''',
        """)


if __name__ == '__main__':
    startTime = time()
    print(f'Started at: ---{timedelta(seconds=startTime)}---')

    tempDenss = []
    for rT in redTemps:
        for d in densities:
            tempDens = [rT, d]
            tempDenss.append(tempDens)

    with Pool(numPools) as p:
        p.map(model, tempDenss)

    cRDir = join(simResultDir, r'compResult')
    cRs = importCompResults(cRDir)
    plotCompResults(cRs, cRDir)
    print(f"""
    Results plotted.
    File saved at: {cRDir}
    """)

    endTime = time() - startTime
    print(f'Time to run: ---{timedelta(seconds=endTime)}---')

    quit()
