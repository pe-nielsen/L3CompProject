from os import listdir
from os.path import isfile, join
from resultClasses import compResult
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm
from time import time
from datetime import timedelta
import pickle
import bz2
from pathlib import Path
import numpy as np

# plt.style.use('seaborn-deep')
# plt.style.use(r'PaperDoubleFig.mplstyle')


def plotCompResults(compResults, cR_parent_dir):
    plt.style.use('seaborn-deep')
    plt.style.use(r'PaperDoubleFig.mplstyle')
    # first: unpack cRs into lists of results for each reduced temp
    # lists of intEn_perPart, pressure_minRhokT

    print(f'compRes: {compResults}')

    # create list of redTemps:
    redTemps = []
    for cR in compResults:
        # print(f'cR.density: {cR.redTemp} ')
        if cR.redTemp not in redTemps:
            redTemps.append(cR.redTemp)

    fig, (axU, axP) = plt.subplots(nrows=2, ncols=1, figsize=[6.4, 8])
    axU.set(
        xlabel=r'volume fraction',
        ylabel=r'configuration energy per particle, u'
    )
    axP.set(
        xlabel=r'volume fraction',
        ylabel=r'P-$\rho$kT'
    )
    axU.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
    axP.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)

    cm1 = mcol.LinearSegmentedColormap.from_list("BlueToRed", ["b", "r"])
    cnorm = mcol.Normalize(vmin=0, vmax=len(redTemps))  # normalising the colourmap
    cpick = cm.ScalarMappable(norm=cnorm, cmap=cm1)  # object which maps redTemp to colour
    # cpick.set_array([])

    for i in range(len(redTemps)):
        rT = redTemps[i]
        intEns_perPart = []
        pressures_minRhokT = []
        densities = []
        for j in range(len(compResults)):
            cR = compResults[j]
            if cR.redTemp == rT:
                intEns_perPart.append(cR.intEn_perPart)
                pressures_minRhokT.append(cR.pressure_minRhokT)
                densities.append(cR.density)

        axU.plot(densities, intEns_perPart,
                 marker='x', ms=4, lw=1, color=cpick.to_rgba(rT))
        axP.plot(densities, pressures_minRhokT,
                 marker='x', ms=4, lw=1, color=cpick.to_rgba(rT))
    fig.savefig(join(cR_parent_dir, r'result.png'), format='png', dpi=1200)
    # fig.show()

def importCompResults(cR_parent_dir):
    compResults = []
    redTemp_dirs = [f for f in listdir(cR_parent_dir) if not isfile(join(cR_parent_dir, f)) and 'redTemp' in f]

    for i in range(len(redTemp_dirs)):
        rT_dir = redTemp_dirs[i]
        rT_path = join(cR_parent_dir, rT_dir)
        data_files = [f for f in listdir(rT_path) if isfile(join(rT_path, f)) and ('den' in f)]
        for j in range(len(data_files)):
            df = data_files[j]
            data_path = join(rT_path, df)
            infile = bz2.BZ2File(data_path, 'r')
            cR = pickle.load(infile)
            infile.close()
            compResults.append(cR)
    print('Computation results successfully imported')
    return compResults


if __name__ == '__main__':
    startTime = time()

    cR_parent_dir = r'simulationResults/smallRunIsotherms3/compResultLargedr'
    compResults = importCompResults(cR_parent_dir)
    plotCompResults(compResults, cR_parent_dir)

    endTime = time() - startTime
    print(f'---{timedelta(seconds=endTime)}---')
    quit()




