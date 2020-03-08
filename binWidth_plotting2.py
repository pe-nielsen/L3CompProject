from os import listdir
from os.path import isfile, join
from resultClasses import compResult
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm
from time import time
from datetime import timedelta
import pickle
from pathlib import Path
import numpy as np

plt.style.use('seaborn-deep')
plt.style.use(r'C:\Users\splb68\comp_proj\scripts\PaperDoubleFig.mplstyle')

def plotCompResults(binWidthFactors, compResults):
    fig, (axU, axP) = plt.subplots(nrows=2, ncols=1, figsize=[6.4, 8])
    axU.set(
        xlabel=r'bin width factor',
        ylabel=r'configuration energy per particle, u'
    )
    axP.set(
        xlabel=r'bin width factor',
        ylabel=r'P-$\rho$kT'
    )
    axU.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
    axP.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)

    for i in range(len(compResults)):
        bWF = binWidthFactors[i]
        cR = compResults[i]

        axU.plot(bWF, cR.intEn_perPart,
                 marker='x', ms=6, c='k')
        axP.plot(bWF, cR.pressure_minRhokT,
                 marker='x', ms=6, c='k')

    fig.savefig(r'C:\Users\splb68\comp_proj\simulationResults\binWidthTesting\compResultWithBWF3\afo_binWidth.png', format='png', dpi=1200)
    fig.show()

def importCompResults(cR_parent_dir):
    binWidthFactors = []
    compResults = []
    data_files = [f for f in listdir(cR_parent_dir) if isfile(join(cR_parent_dir, f)) and ('den' in f)]

    for j in range(len(data_files)):
        df = data_files[j]
        data_path = join(cR_parent_dir, df)
        infile = open(data_path, 'rb')
        bWF, cR = pickle.load(infile)
        infile.close()
        compResults.append(cR)
        binWidthFactors.append(bWF)

    return binWidthFactors, compResults


if __name__ == '__main__':
    startTime = time()
    cR_parent_dir = r'C:\Users\splb68\comp_proj\simulationResults\binWidthTesting\compResultWithBWF3'

    binWidthFactors, compResults = importCompResults(cR_parent_dir)
    plotCompResults(binWidthFactors, compResults)

    endTime = time() - startTime
    print(f'---{timedelta(seconds=endTime)}---')
    quit()