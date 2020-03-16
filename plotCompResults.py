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
        xlabel=r'$\rho$*',
        ylabel=r'configuration energy per particle, u',
        # yscale='log'
    )
    axP.set(
        xlabel=r'$\rho$*',
        ylabel=r'P-$\rho$kT',
        # yscale='log'
    )
    axU.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
    axP.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)

    cm1 = mcol.LinearSegmentedColormap.from_list("BlueToRed", ["b", "r"])
    cnorm = mcol.Normalize(vmin=min(redTemps), vmax=max(redTemps))  # normalising the colourmap
    cpick = cm.ScalarMappable(norm=cnorm, cmap=cm1)  # object which maps redTemp to colour
    # cpick.set_array([])

    for i in range(len(redTemps)):
        rT = redTemps[i]
        intEns_perPart = []
        pressures_minRhokT = []
        densities = []
        numberDensities = []
        errorsIn_u = []
        errorsIn_pressure = []
        for j in range(len(compResults)):
            cR = compResults[j]
            if cR.redTemp == rT:
                intEns_perPart.append(cR.intEn_perPart)
                pressures_minRhokT.append(cR.pressure_minRhokT)
                densities.append(cR.density)
                numberDensity = cR.numParts/(cR.containerLength**3)

                numberDensities.append(numberDensity)

                # sum_ErrorGiSq = np.sum(cR.PCFstdev**2)
                errorsIn_u.append(cR.intEnstdev)
                errorsIn_pressure.append(cR.pressurestdev)
                print(f'at redTemp {cR.redTemp}, den {cR.density}:\n'
                      f'    U*: {cR.intEn_perPart}\n'
                      f'    +/-: {cR.intEnstdev}\n'
                      f'    p*: {cR.pressure_minRhokT}\n'
                      f'    +/-: {cR.pressurestdev}\n')

        axU.plot(densities, intEns_perPart,
                 ls='', marker='x', ms=4, lw=1, color=cpick.to_rgba(rT))
        # axU.errorbar(densities, intEns_perPart, yerr=errorsIn_u,
        #              marker='', ls='', ecolor=cpick.to_rgba(rT))
        axP.plot(densities, pressures_minRhokT,
                 ls='', marker='x', ms=4, lw=1, color=cpick.to_rgba(rT))
        # axP.errorbar(densities, pressures_minRhokT, yerr=errorsIn_pressure,
        #              marker='', ls='', ecolor=cpick.to_rgba(rT))


    # ###NIST RESULTS###
    NIST_rho = [1e-3, 3e-3, 5e-3, 7e-3, 9e-3, 7.76e-1, 7.80e-1, 8.2e-1, 8.60e-1, 9e-1]
    NIST_T085_u = [-1.0317E-2, -3.1019E-2, -5.1901E-2, -7.2834E-2, -9.3973E-02, -5.5121E+00, -5.5386E+00, -5.7947E+00, -6.0305E+00, -6.2391E+00]
    NIST_T085_p = [8.4402E-04, 2.4965E-03, 4.1003E-03, 5.6565E-03, 7.1641E-03, 6.7714E-03, 4.7924E-02, 5.5355E-01, 1.2660E+00, 2.2314E+00]
    NIST_T09_u = [-9.9165E-03, -2.9787E-02, -4.9771E-02, -6.9805E-02, -8.9936E-02, -5.4689E+00, -5.4956E+00, -5.7456E+00, -5.9753E+00, -6.1773E+00]
    NIST_T09_p = [8.9429E-04, 2.6485E-03, 4.3569E-03, 6.0193E-03, 7.6363E-03, 2.4056E-01, 2.7851E-01, 8.2386E-01, 1.5781E+00, 2.5848E+00]

    axU.plot(NIST_rho, NIST_T085_u,
             ls='', marker='o', ms=4, lw=1, color='dodgerblue')
    axU.plot(NIST_rho, NIST_T09_u,
             ls='', marker='o', ms=4, lw=1, color='crimson')
    axP.plot(NIST_rho, NIST_T085_p,
             ls='', marker='o', ms=4, lw=1, color='dodgerblue')
    axP.plot(NIST_rho, NIST_T09_p,
             ls='', marker='o', ms=4, lw=1, color='crimson')





    fig.savefig(join(cR_parent_dir, r'result.png'), format='png', dpi=600)
#    fig.show()

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

    cR_parent_dir = r'resultsNCC/NIST_params_2/compResult'
    compResults = importCompResults(cR_parent_dir)
    plotCompResults(compResults, cR_parent_dir)

    endTime = time() - startTime
    print(f'---{timedelta(seconds=endTime)}---')
    quit()




