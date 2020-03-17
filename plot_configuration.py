from os.path import isfile, join

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from resultClasses import simResult, compResult
from commonFunctions import getDistanceSquared, redLJ, diff_redLJ

from time import time
from datetime import timedelta
from pathlib import Path
import pickle
import bz2


import numpy as np
np.random.seed(0)
plt.style.use('seaborn-deep')
plt.style.use(r'PaperDoubleFig.mplstyle')


def plotConfig(sR):
    config = sR.configs[:, :, 100]  # guessing right indices: [particle num, xyz, sample num]
    containerLength = sR.containerLength
    partRad = sR.partRad
    plotTitle = f'density: {sR.density}, redTemp: {sR.redTemp:.2E}'

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set(
        xlabel='x',
        ylabel='y',
        zlabel='z',
        title=plotTitle,
        xlim=(0, containerLength),
        ylim=(0, containerLength),
        zlim=(0, containerLength),
    )
    xs = config[:, 0]  # picking out the x, y, z values for the locations for plotting
    ys = config[:, 1]
    zs = config[:, 2]

    size = (2*partRad)**4
    ax.scatter(xs, ys, zs,
               s=size, color='grey')
    fig.show()


def plotPCF(cR):
    fig, (axPCF) = plt.subplots(nrows=1, ncols=1, figsize=[6.4, 4.8])
    plotTitle = f'density: {cR.density}, redTemp: {cR.redTemp:.2E}'

    axPCF.set(
        xlabel=r'r/$\sigma{}$',
        ylabel=r'g(r/$\sigma{}$)',
        title=plotTitle
    )
    axPCF.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)

    axPCF.plot(cR.PCFradii, cR.PCF,
               marker='', ms=4, c='k', ls='-', lw=1.5)

    fig.show()


def importResult(filepath):
    # infile = open(filepath, 'rb')
    # result = pickle.load(infile)
    # infile.close()

    infile = bz2.BZ2File(filepath, 'r')
    result = pickle.load(infile)
    infile.close()

    return result


if __name__ == '__main__':
    startTime = time()

    cR_filepath = r'resultsNCC/dataCollection/compResult/redTemp10/den0.1'
    # sR_filepath = r'C:\Users\Peter Nielsen\Documents\l3\comp_proj\simulationResults\testRun\simResult_den75'

    cR = importResult(cR_filepath)
    # sR = importResult(sR_filepath)


    # plotConfig(sR)  # this iff simResult
    plotPCF(cR)  # this iff compResult

