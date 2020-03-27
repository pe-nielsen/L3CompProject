from os import listdir
from os.path import isfile, join
from resultClasses import compResult
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import colorsys
import matplotlib.cm as cm
from matplotlib import gridspec
from time import time
from datetime import timedelta
import pickle
import bz2
from pathlib import Path
import numpy as np
from numpy import convolve
plt.style.use('seaborn-deep')
plt.style.use(r'PaperDoubleFig.mplstyle')


def ReePressure(density, redTemp, sigma, depth):
    numDens = density/(sigma**3)
    beta = 1/(redTemp*depth)

    x = density/(redTemp**0.25)
    Bs = [3.629, 7.2641, 10.4924, 11.459, 0, 0, 0, 0, 0, 2.1769]
    Cs = [5.3692, 6.5797, 6.1745, -4.2685, 1.6841]
    Ds = [-3.4921, 18.6980, -35.5049, 31.8151, -11.1953]

    term1 = 1 + (Bs[0]*x**1) + (Bs[1]*x**2) + (Bs[2]*x**3) + (Bs[3]*x**4) + (Bs[9]*x**10)
    term2 = (1/redTemp)**(1/2) * np.sum([(i+1)*Cs[i]*x**(i+1) for i in range(len(Cs))])
    term3 = (1/redTemp) * np.sum([Ds[i]*x**(i+1) for i in range(len(Cs))])

    betaPorho = term1 + term2 + term3
    return betaPorho


def plotCompResults(compResults, cR_parent_dir):
    plt.style.use('seaborn-deep')
    plt.style.use(r'PaperDoubleFig.mplstyle')

    print(f'compRes: {compResults}')

    # create list of redTemps:
    redTemps = []
    for cR in compResults:
        # print(f'cR.density: {cR.redTemp} ')
        if cR.redTemp not in redTemps:
            redTemps.append(cR.redTemp)

    # fig, (axU, axP) = plt.subplots(nrows=2, ncols=1, figsize=[6.4, 8])
    figU, axU = plt.subplots(nrows=1, ncols=1, figsize=[6.4, 4.8])
    figP, axP = plt.subplots(nrows=1, ncols=1, figsize=[6.4, 4.8])

    left, bottom, width, height = [0.275, 0.55, 0.5, 0.35]
    axPinset = figP.add_axes([left, bottom, width, height])

    axU.set(
        xlabel=r'$\rho{}*$',
        # ylabel=r'configuration energy per particle, u',
        ylabel=r'$u*$',

        # yscale='log',
        # xscale='log',
    )
    axP.set(
        xlabel=r'$\rho{}*$',
        # ylabel=r'P-$\rho$kT',
        ylabel=r'$p*$',
        # yscale='log',
        # xscale='log',
    )
    # axPinset.set(
    #     xlabel=r'$\rho{}*$',
    #     # ylabel=r'P-$\rho$kT',
    #     ylabel=r'$p*$',
    #     ylim=[-1.25, 0.5],
    #     xlim=[0.5, 0.7],
    # )
    axPinset.set(
        xlabel=r'$\rho{}*$',
        # ylabel=r'P-$\rho$kT',
        ylabel=r'$p*$',
        ylim=[-1.25, 0.3],
        xlim=[0, 0.9],
    )

    axU.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
    axP.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
    axPinset.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)


    cm1 = mcol.LinearSegmentedColormap.from_list("BlueToRed", ["b", "r"])
    # cnorm = mcol.Normalize(vmin=min(redTemps), vmax=max(redTemps))  # normalising the colourmap
    cnorm = mcol.LogNorm(vmin=min(redTemps), vmax=max(redTemps))  # normalising the colourmap
    cpick = cm.ScalarMappable(norm=cnorm, cmap=cm1)  # object which maps redTemp to colour
    # cpick.set_array([])

    redTemps.sort()
    for i in range(len(redTemps)):
        rT = redTemps[i]

        ReePressures = []

        intEns_perPart = []
        pressures_minRhokT = []
        densities = []
        numberDensities = []
        errorsIn_u = []
        errorsIn_pressure = []
        for j in range(len(compResults)):
            cR = compResults[j]
            if cR.redTemp == rT:
                cutOffR = cR.containerLength/2
                u_LRC = (8/9)*np.pi*cR.density*((1/cutOffR)**9 - 3*(1/cutOffR)**3)
                p_LRC = (32/9)*np.pi*(cR.density**2)*((1/cutOffR)**9 - (3/2)*(1/cutOffR)**3)
                print(f'uLRC: {u_LRC},\n'
                      f'pLRC: {p_LRC}\n')

                intEns_perPart.append(cR.intEn_perPart + u_LRC)
                pressures_minRhokT.append(cR.pressure_minRhokT + p_LRC)
                densities.append(cR.density)
                numberDensity = cR.numParts/(cR.containerLength**3)

                numberDensities.append(numberDensity)

                # sum_ErrorGiSq = np.sum(cR.PCFstdev**2)
                errorsIn_u.append(cR.intEnstdev)
                errorsIn_pressure.append(cR.pressurestdev)

                ReePressures.append(ReePressure(cR.density, cR.redTemp, cR.sigma, cR.wellDepth))

                # print(f'at redTemp {cR.redTemp}, den {cR.density}:\n'
                #       f'    U*: {cR.intEn_perPart:.5E}\n'
                #       f'    +/-: {cR.intEnstdev:.5E}\n'
                #       f'    p*: {cR.pressure_minRhokT:.5E}\n'
                #       f'    +/-: {cR.pressurestdev:.5E}\n'
                #       )

        axU.plot(densities, intEns_perPart,
                 ls=':', marker='.', ms=4, lw=1, color=cpick.to_rgba(rT), label=f'T*={rT}')
        axP.plot(densities, pressures_minRhokT,
                 ls=':', marker='.', ms=4, lw=1, color=cpick.to_rgba(rT), label=f'T*={rT}')

        axPinset.plot(densities, pressures_minRhokT,
                      ls=':', marker='.', ms=4, lw=1.5, color=cpick.to_rgba(rT), label=f'T*={rT}')
        if rT == 0.75 or rT == 1.5:
            axPinset.errorbar(densities, pressures_minRhokT, yerr=errorsIn_pressure,
                         marker='', ls='', lw=1, ecolor=cpick.to_rgba(rT), capsize=2)

        # axP.plot(densities, np.array(ReePressures),
        #          ls='-', lw=1.5, c=cpick.to_rgba(rT))

        axU.errorbar(densities, intEns_perPart, yerr=errorsIn_u,
                     marker='', ls='', lw=1, ecolor=cpick.to_rgba(rT), capsize=2)
        axP.errorbar(densities, pressures_minRhokT, yerr=errorsIn_pressure,
                     marker='', ls='', lw=1, ecolor=cpick.to_rgba(rT), capsize=2)

    # Johnson_density = [0.05, 0.6, 0.7, 0.8, 0.9, 0.95]
    # Johnson_p = [0.0368, -2.272, 0.015, 1.011, 3.28, 5.131]
    # Johnson_u = [-0.483, -4.224, -4.887, -5.535, -6.055, -6.2321]
    #
    # axU.plot(Johnson_density, Johnson_u,
    #          ls=':', marker='x', ms=4, lw=1.5, color='k', label=f'T*=1.0')
    # axP.plot(Johnson_density, Johnson_p,
    #          ls=':', marker='x', ms=4, lw=1.5, color='k', label=f'T*=1.0')


    # axU.legend(loc='upper left',
    #             frameon=True, framealpha=1, edgecolor='black', fancybox=False)
    # axP.legend(loc='upper left',
    #             frameon=True, framealpha=1, edgecolor='black', fancybox=False)

    figU.savefig(join(cR_parent_dir, r'excessEnergy.png'), format='png', dpi=1200)
    figP.savefig(join(cR_parent_dir, r'excessPressure.png'), format='png', dpi=1200)

    # figU.show()
    figP.show()


def plotJohnsoncomparison(compResults, saveLocation):
    plt.style.use('seaborn-deep')
    plt.style.use(r'PaperDoubleFig.mplstyle')

    print(f'compRes: {compResults}')

    # create list of redTemps:
    redTemps = []
    for cR in compResults:
        if cR.redTemp not in redTemps:
            redTemps.append(cR.redTemp)

    cm1 = mcol.LinearSegmentedColormap.from_list("BlueToRed", ["b", "r"])
    cnorm = mcol.LogNorm(vmin=min(redTemps), vmax=max(redTemps))  # normalising the colourmap
    cpick = cm.ScalarMappable(norm=cnorm, cmap=cm1)  # object which maps redTemp to colour

    figsizes = [6.4, 7]

    fig = plt.figure(figsize=figsizes)
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1,1])
    gs.update(hspace=0.0)
    axU = fig.add_subplot(gs[0])
    axP = fig.add_subplot(gs[1], sharex=axU)

    axU.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True,
                    labelbottom=False)
    axP.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)

    axU.set(
        ylabel=r'$u*$'
    )
    axP.set(
        xlabel=r'$\rho{}*$',
        ylabel=r'$p*$',
    )

    redTemps.sort()
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
                cutOffR = cR.containerLength / 2
                u_LRC = (8 / 9) * np.pi * cR.density * ((1 / cutOffR) ** 9 - 3 * (1 / cutOffR) ** 3)
                p_LRC = (32 / 9) * np.pi * (cR.density ** 2) * ((1 / cutOffR) ** 9 - (3 / 2) * (1 / cutOffR) ** 3)
                print(f'uLRC: {u_LRC},\n'
                      f'pLRC: {p_LRC}\n')

                intEns_perPart.append(cR.intEn_perPart + u_LRC)
                pressures_minRhokT.append(cR.pressure_minRhokT + p_LRC)
                densities.append(cR.density)
                numberDensity = cR.numParts / (cR.containerLength ** 3)
                numberDensities.append(numberDensity)
                errorsIn_u.append(cR.intEnstdev)
                errorsIn_pressure.append(cR.pressurestdev)

        limIndex = -5
        axU.plot(densities[:limIndex], intEns_perPart[:limIndex],
                 ls='', marker='o', ms=4, mfc='none', color=cpick.to_rgba(rT))
        axP.plot(densities[:limIndex], pressures_minRhokT[:limIndex],
                 ls='', marker='o', ms=4, mfc='none', color=cpick.to_rgba(rT))

        axU.errorbar(densities[:limIndex], intEns_perPart[:limIndex], yerr=errorsIn_u[:limIndex],
                     marker='', ls='', lw=1, ecolor=cpick.to_rgba(rT), capsize=2)
        axP.errorbar(densities[:limIndex], pressures_minRhokT[:limIndex], yerr=errorsIn_pressure[:limIndex],
                     marker='', ls='', lw=1, ecolor=cpick.to_rgba(rT), capsize=2)

    Johnson_density100 = [0.05, 0.6, 0.7, 0.8, 0.9, 0.95]
    Johnson_p100 = [0.0368, -2.272, 0.015, 1.011, 3.28, 5.131]
    Johnson_u100 = [-0.483, -4.224, -4.887, -5.535, -6.055, -6.2321]

    Johnson_density075 = [0.005, 0.005, 0.01, 0.01, 0.7, 0.8, 0.9]
    Johnson_p075 = [0.00359, 0.00358, 0.00684, 0.00685, -0.834, -0.256, 1.503]
    Johnson_u075 = [-0.0554, -0.0613, -0.115, -0.113, -5.070, -5.765, -6.365]

    axU.plot(Johnson_density100, Johnson_u100,
             ls='', marker='x', ms=10, c=cpick.to_rgba(1.0))
    axP.plot(Johnson_density100, Johnson_p100,
             ls='', marker='x', ms=10, c=cpick.to_rgba(1.0))
    axU.plot(Johnson_density075, Johnson_u075,
             ls='', marker='x', ms=10, c=cpick.to_rgba(0.75))
    axP.plot(Johnson_density075, Johnson_p075,
             ls='', marker='x', ms=10, c=cpick.to_rgba(0.7))

    fig.savefig(join(saveLocation, r'Johnson.png'), format='png', dpi=1200)
    fig.show()


def plotNISTcomparison(compResults, saveLocation):
    plt.style.use('seaborn-deep')
    plt.style.use(r'PaperDoubleFig.mplstyle')

    print(f'compRes: {compResults}')

    # create list of redTemps:
    redTemps = []
    for cR in compResults:
        if cR.redTemp not in redTemps:
            redTemps.append(cR.redTemp)

    cm1 = mcol.LinearSegmentedColormap.from_list("BlueToRed", ["b", "r"])
    cnorm = mcol.LogNorm(vmin=min(redTemps), vmax=max(redTemps))  # normalising the colourmap
    cpick = cm.ScalarMappable(norm=cnorm, cmap=cm1)  # object which maps redTemp to colour

    figsizes = [6.4, 4.8]
    height_ratios = [3, 2]

    figU = plt.figure(figsize=figsizes)
    gsU = figU.add_gridspec(2, 2, height_ratios=height_ratios)
    gsU.update(hspace=0.0)
    axU1 = figU.add_subplot(gsU[0])
    axU2 = figU.add_subplot(gsU[1])
    axRU1 = figU.add_subplot(gsU[2], sharex=axU1)
    axRU2 = figU.add_subplot(gsU[3], sharex=axU2)


    figP = plt.figure(figsize=figsizes)
    gsP = figP.add_gridspec(2, 2, height_ratios=height_ratios)
    gsP.update(hspace=0.0)
    axP1 = figP.add_subplot(gsP[0])
    axP2 = figP.add_subplot(gsP[1])
    axRP1 = figP.add_subplot(gsP[2], sharex=axP1)
    axRP2 = figP.add_subplot(gsP[3], sharex=axP2)


    axU1.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True,
                    labelbottom=False)
    axU2.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True,
                     labelbottom=False)
    axP1.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True,
                     labelbottom=False)
    axP2.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True,
                     labelbottom=False)

    axRU1.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
    axRU2.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
    axRP1.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
    axRP2.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)


    axU1.set(
        ylabel=r'$u*\times{}1000$',
    )
    axU2.set(
        ylabel=r'$u*$',
    )
    axP1.set(
        # xlabel=r'$\rho{}*\times{}1000$',
        ylabel=r'$p*\times{}1000$',
    )
    axP2.set(
        # xlabel=r'$\rho{}*$',
        ylabel=r'$p*$',
    )

    axRU1.set(
        xlabel=r'$\rho{}*\times{}1000$',
        ylabel=r'$\Delta{}u*\times{}1000$',
    )
    axRU2.set(
        xlabel=r'$\rho{}*$',
        ylabel=r'$\Delta{}u*$',
    )
    axRP1.set(
        xlabel=r'$\rho{}*\times{}1000$',
        ylabel=r'$\Delta{}p*\times{}1000$',
    )
    axRP2.set(
        xlabel=r'$\rho{}*$',
        ylabel=r'$\Delta{}p*$',
    )

    # ###NIST RESULTS###
    NIST_rho = [1e-3, 3e-3, 5e-3, 7e-3, 9e-3, 7.76e-1, 7.80e-1, 8.2e-1, 8.60e-1, 9e-1]
    NIST_T085_u = [-1.0317E-2, -3.1019E-2, -5.1901E-2, -7.2834E-2, -9.3973E-02, -5.5121E+00, -5.5386E+00, -5.7947E+00,
                   -6.0305E+00, -6.2391E+00]
    NIST_T085_p = [8.4402E-04, 2.4965E-03, 4.1003E-03, 5.6565E-03, 7.1641E-03, 6.7714E-03, 4.7924E-02, 5.5355E-01,
                   1.2660E+00, 2.2314E+00]
    NIST_T09_u = [-9.9165E-03, -2.9787E-02, -4.9771E-02, -6.9805E-02, -8.9936E-02, -5.4689E+00, -5.4956E+00,
                  -5.7456E+00, -5.9753E+00, -6.1773E+00]
    NIST_T09_p = [8.9429E-04, 2.6485E-03, 4.3569E-03, 6.0193E-03, 7.6363E-03, 2.4056E-01, 2.7851E-01, 8.2386E-01,
                  1.5781E+00, 2.5848E+00]
    #
    axU1.plot(1000 * np.array(NIST_rho[:5]), 1000 * np.array(NIST_T085_u[:5]),
              ls='', marker='x', ms=4, lw=1, color=cpick.to_rgba(0.85))
    axU1.plot(1000 * np.array(NIST_rho[:5]), 1000 * np.array(NIST_T09_u[:5]),
              ls='', marker='x', ms=4, lw=1, color=cpick.to_rgba(0.90))
    axP1.plot(1000 * np.array(NIST_rho[:5]), 1000 * np.array(NIST_T085_p[:5]),
              ls='', marker='x', ms=4, lw=1, color=cpick.to_rgba(0.85))
    axP1.plot(1000 * np.array(NIST_rho[:5]), 1000 * np.array(NIST_T09_p[:5]),
              ls='', marker='x', ms=4, lw=1, color=cpick.to_rgba(0.90))

    axU2.plot(NIST_rho[5:], NIST_T085_u[5:],
              ls='', marker='x', ms=4, lw=1, color=cpick.to_rgba(0.85))
    axU2.plot(NIST_rho[5:], NIST_T09_u[5:],
              ls='', marker='x', ms=4, lw=1, color=cpick.to_rgba(0.90))
    axP2.plot(NIST_rho[5:], NIST_T085_p[5:],
              ls='', marker='x', ms=4, lw=1, color=cpick.to_rgba(0.85))
    axP2.plot(NIST_rho[5:], NIST_T09_p[5:],
              ls='', marker='x', ms=4, lw=1, color=cpick.to_rgba(0.90))

    # plot line across 0 on residuals
    axRP1.axhline(y=0, xmin=0, xmax=10, lw=1, ls='--', c='k', alpha=1)
    axRP2.axhline(y=0, xmin=0, xmax=1, lw=1, ls='--', c='k', alpha=1)
    axRU1.axhline(y=0, xmin=0, xmax=10, lw=1, ls='--', c='k', alpha=1)
    axRU2.axhline(y=0, xmin=0, xmax=1, lw=1, ls='--', c='k', alpha=1)



    # added long range corrections!

    redTemps.sort()
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
                cutOffR = cR.containerLength/2
                u_LRC = (8/9)*np.pi*cR.density*((1/cutOffR)**9 - 3*(1/cutOffR)**3)
                p_LRC = (32/9)*np.pi*(cR.density**2)*((1/cutOffR)**9 - (3/2)*(1/cutOffR)**3)
                print(f'uLRC: {u_LRC},\n'
                      f'pLRC: {p_LRC}\n')

                intEns_perPart.append(cR.intEn_perPart + u_LRC)
                pressures_minRhokT.append(cR.pressure_minRhokT + p_LRC)
                densities.append(cR.density)
                numberDensity = cR.numParts / (cR.containerLength ** 3)
                numberDensities.append(numberDensity)
                errorsIn_u.append(cR.intEnstdev)
                errorsIn_pressure.append(cR.pressurestdev)

        axU1.plot(1000*np.array(densities[:5]), 1000*np.array(intEns_perPart[:5]),
                 ls='', marker='o', ms=4, mfc='none', color=cpick.to_rgba(rT), label=f'T*={rT}')
        axP1.plot(1000*np.array(densities[:5]), 1000*np.array(pressures_minRhokT[:5]),
                 ls='', marker='o', ms=4, mfc='none', color=cpick.to_rgba(rT), label=f'T*={rT}')

        axU1.errorbar(1000*np.array(densities[:5]), 1000*np.array(intEns_perPart[:5]), yerr=1000*np.array(errorsIn_u[:5]),
                     marker='', ls='', lw=1, ecolor=cpick.to_rgba(rT), capsize=2)
        axP1.errorbar(1000*np.array(densities[:5]), 1000*np.array(pressures_minRhokT[:5]), yerr=1000*np.array(errorsIn_pressure[:5]),
                     marker='', ls='', lw=1, ecolor=cpick.to_rgba(rT), capsize=2)

        axU2.plot(densities[5:], intEns_perPart[5:],
                  ls='', marker='o', ms=4, mfc='none', color=cpick.to_rgba(rT), label=f'T*={rT}')
        axP2.plot(densities[5:], pressures_minRhokT[5:],
                  ls='', marker='o', ms=4, mfc='none', color=cpick.to_rgba(rT), label=f'T*={rT}')

        axU2.errorbar(densities[5:], intEns_perPart[5:], yerr=errorsIn_u[5:],
                      marker='', ls='', lw=1, ecolor=cpick.to_rgba(rT), capsize=2)
        axP2.errorbar(densities[5:], pressures_minRhokT[5:], yerr=errorsIn_pressure[5:],
                      marker='', ls='', lw=1, ecolor=cpick.to_rgba(rT), capsize=2)

        residualU = 0
        residualP = 0
        if rT == 0.85:
            residualU = np.array(intEns_perPart) - np.array(NIST_T085_u)
            residualP = np.array(pressures_minRhokT) - np.array(NIST_T085_p)
        elif rT == 0.90:
            residualU = np.array(intEns_perPart) - np.array(NIST_T09_u)
            residualP = np.array(pressures_minRhokT) - np.array(NIST_T09_p)

        axRU1.plot(1000 * np.array(densities[:5]), 1000*np.array(residualU[:5]),
                  ls='', marker='o', ms=4, color=cpick.to_rgba(rT), label=f'T*={rT}')
        axRP1.plot(1000 * np.array(densities[:5]), 1000*np.array(residualP[:5]),
                  ls='', marker='o', ms=4, color=cpick.to_rgba(rT), label=f'T*={rT}')
        axRU2.plot(np.array(densities[5:]), residualU[5:],
                  ls='', marker='o', ms=4, color=cpick.to_rgba(rT), label=f'T*={rT}')
        axRP2.plot(np.array(densities[5:]), residualP[5:],
                  ls='', marker='o', ms=4, color=cpick.to_rgba(rT), label=f'T*={rT}')

    figU.savefig(join(saveLocation, r'NIST_u.png'), format='png', dpi=1200)
    figP.savefig(join(saveLocation, r'NIST_p.png'), format='png', dpi=1200)


    # figU.show()
    figP.show()


def plotConvergence():
    # import convergence data for just a single rT

    cm1 = mcol.LinearSegmentedColormap.from_list("BlueToRed", ["b", "r"])
    cnorm = mcol.LogNorm(vmin=0.75, vmax=1.5)  # normalising the colourmap
    cpick = cm.ScalarMappable(norm=cnorm, cmap=cm1)  # object which maps redTemp to colour

    fig = plt.figure(figsize=[6.4, 4.8])
    gs = fig.add_gridspec(2, 1)
    gs.update(hspace=0.0)
    axUp = fig.add_subplot(gs[0])
    axLow = fig.add_subplot(gs[1])

    axUp.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True, labelbottom=False)
    axLow.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
    axUp.set(
        # xlabel='MC Moves (in millions)',
        ylabel=r'$g(r^{*}=1)$',
    )
    axLow.set(
        xlabel='MC moves (in millions)',
        ylabel=r'$g(r^{*}=1)$',
    )

    parent_folder = r'report/figures/convergence/compResult'
    saveFigureTo = r'report/figures/convergence'
    rT_folders = ['redTemp0.75', 'redTemp1.5']
    for rT_subfolder in rT_folders:
        gAtSigma_LIST = []
        densities = []
        redTemp = 0
        axis = None
        if '1.5' in rT_subfolder:
            redTemp = 1.5
            axis = axUp
        else:
            redTemp = 0.75
            axis = axLow

        rT_folder = join(parent_folder, rT_subfolder)
        data_files = [f for f in listdir(rT_folder) if isfile(join(rT_folder, f)) and ('gAtSigma' in f)]
        for j in range(len(data_files)):
            df = data_files[j]
            data_path = join(rT_folder, df)
            infile = bz2.BZ2File(data_path, 'r')
            gAS = pickle.load(infile)
            infile.close()
            gAtSigma_LIST.append(gAS)
            densities.append(float(df[-3:]))
            # print(f'density: {df[-3:]}')

        for i in range(len(gAtSigma_LIST)):
            gAS = gAtSigma_LIST[i]
            density = densities[i]
            numSamples = gAS.shape[0]
            numMoves = 2e7
            movesPerSample = numMoves/numSamples
            MCMoves = [j*movesPerSample for j in range(numSamples)]

            # gASMA = movingaverage(gAS.tolist(), 1000)
            gASCA = cumulativeaverage(gAS.tolist())  # cumulative average

            window = 500
            gASMA = movingaverage(gAS.tolist(), window)

            print(f"""
            list: {gAS}
            rA_list: {len(gASCA)}
            difference: {np.array(gASCA) - gAS[numSamples - len(gASCA):]}
            """)

            MCMoves_cropped = np.array(MCMoves[numSamples - len(gASMA):])

            if density == max(densities):
                linestyle = '-'
            elif density == min(densities):
                linestyle = ':'
            else:
                linestyle = '-.'

            baseColour = cpick.to_rgba(redTemp)
            lightenBy = 0.5
            try:
                c = mcol.cnames[baseColour]
            except:
                c = baseColour
            c = np.array(colorsys.rgb_to_hls(*mcol.to_rgb(c)))
            lighterColour = colorsys.hls_to_rgb(c[0], 1-lightenBy*(1-c[1]), c[2])


            axis.plot(MCMoves_cropped / 1e6, gASMA,
                      ls='-', marker='', lw=0.5, color=lighterColour,
                      label=r'$\rho{}$*' + f'={density}')

            axis.plot(np.array(MCMoves)/1e6, gASCA,
                      ls=linestyle, marker='', lw=1, color=cpick.to_rgba(redTemp), label=r'$\rho{}$*'+f'={density}')

    # ax.legend(loc='bottom right')
    fig.show()
    fig.savefig(join(saveFigureTo, r'gAtSigma_convergence.png'), format='png', dpi=1200)


def plotBinWidthEffect(BWFcRs):
    cm1 = mcol.LinearSegmentedColormap.from_list("BlueToRed", ["b", "r"])
    cnorm = mcol.LogNorm(vmin=0.75, vmax=1.5)  # normalising the colourmap
    cpick = cm.ScalarMappable(norm=cnorm, cmap=cm1)  # object which maps redTemp to colour

    fig = plt.figure(figsize=[6.4, 4.8])
    gs = fig.add_gridspec(2, 1)
    gs.update(hspace=0.0)
    axU = fig.add_subplot(gs[0])
    axP = fig.add_subplot(gs[1])

    axU.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True,
                     labelbottom=False)
    axP.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
    axU.set(
        # xlabel='MC Moves (in millions)',
        ylabel=r'$u*$',
    )
    axP.set(
        xlabel=r'$\Delta{}r*$',
        ylabel=r'$p*$',
    )

    baseBlue = cpick.to_rgba(0.75)
    baseRed = cpick.to_rgba(1.5)
    lightenBy = 0.5
    try:
        c1 = mcol.cnames[baseBlue]
        c2 = mcol.cnames[baseRed]
    except:
        c1 = baseBlue
        c2 = baseRed
    c1 = np.array(colorsys.rgb_to_hls(*mcol.to_rgb(c1)))
    c2 = np.array(colorsys.rgb_to_hls(*mcol.to_rgb(c2)))

    lighterBlue = colorsys.hls_to_rgb(c1[0], 1 - lightenBy * (1 - c1[1]), c1[2])
    lighterRed = colorsys.hls_to_rgb(c2[0], 1 - lightenBy * (1 - c2[1]), c2[2])



    BWFsHT = []
    intEnPlusErrorHT = []
    intEnMinusErrorHT = []
    pressurePlusErrorHT = []
    pressureMinusErrorHT = []
    
    BWFsLT = []
    intEnPlusErrorLT = []
    intEnMinusErrorLT = []
    pressurePlusErrorLT = []
    pressureMinusErrorLT = []

    for i in range(len(BWFcRs)):
        (BWF, cR) = BWFcRs[i]
        cutOffR = cR.containerLength / 2
        u_LRC = (8 / 9) * np.pi * cR.density * ((1 / cutOffR) ** 9 - 3 * (1 / cutOffR) ** 3)
        p_LRC = (32 / 9) * np.pi * (cR.density ** 2) * ((1 / cutOffR) ** 9 - (3 / 2) * (1 / cutOffR) ** 3)
        # print(f'uLRC: {u_LRC},\n'
        #       f'pLRC: {p_LRC}\n')

        if cR.redTemp == 1.5:
            BWFsHT.append(BWF)
            intEnPlusErrorHT.append(cR.intEn_perPart + u_LRC + cR.intEnstdev)
            intEnMinusErrorHT.append(cR.intEn_perPart + u_LRC - cR.intEnstdev)
            pressurePlusErrorHT.append(cR.pressure_minRhokT + p_LRC + cR.pressurestdev)
            pressureMinusErrorHT.append(cR.pressure_minRhokT + p_LRC - cR.pressurestdev)
        elif cR.redTemp == 0.75:
            BWFsLT.append(BWF)
            intEnPlusErrorLT.append(cR.intEn_perPart + u_LRC + cR.intEnstdev)
            intEnMinusErrorLT.append(cR.intEn_perPart + u_LRC - cR.intEnstdev)
            pressurePlusErrorLT.append(cR.pressure_minRhokT + p_LRC + cR.pressurestdev)
            pressureMinusErrorLT.append(cR.pressure_minRhokT + p_LRC - cR.pressurestdev)
            
        axU.plot(BWF, cR.intEn_perPart + u_LRC,
                 ls='', marker='.', ms=4, color=cpick.to_rgba(cR.redTemp))
        axP.plot(BWF, cR.pressure_minRhokT + p_LRC,
                 ls='', marker='.', ms=4, color=cpick.to_rgba(cR.redTemp))

        print(f'at redTemp: {cR.redTemp},    BWF: {BWF},\n'
              f'u* = {cR.intEn_perPart + u_LRC},\n'
              f'p* = {cR.pressure_minRhokT + p_LRC}')


        axU.errorbar(BWF, cR.intEn_perPart + u_LRC, yerr=cR.intEnstdev,
                     marker='', ls='', lw=1, ecolor=cpick.to_rgba(cR.redTemp), capsize=2)
        axP.errorbar(BWF, cR.pressure_minRhokT + p_LRC, yerr=cR.pressurestdev,
                     marker='', ls='', lw=1, ecolor=cpick.to_rgba(cR.redTemp), capsize=2)

    alpha = 0.1
    axU.fill_between(x=BWFsHT, y1=intEnPlusErrorHT, y2=intEnMinusErrorHT, color=cpick.to_rgba(1.5), alpha=alpha)
    axU.fill_between(x=BWFsLT, y1=intEnPlusErrorLT, y2=intEnMinusErrorLT, color=cpick.to_rgba(0.75), alpha=alpha)
    axP.fill_between(x=BWFsHT, y1=pressurePlusErrorHT, y2=pressureMinusErrorHT, color=cpick.to_rgba(1.5), alpha=alpha)
    axP.fill_between(x=BWFsLT, y1=pressurePlusErrorLT, y2=pressureMinusErrorLT, color=cpick.to_rgba(0.75), alpha=alpha)

    fig.show()
    parent_folder = r'report/figures/binWidth'
    # fig.savefig(join(parent_folder, r'BWFeffect.png'), format='png', dpi=1200)


def plotRDFs(cRs, saveLocation):
    plt.style.use('seaborn-deep')
    plt.style.use(r'PaperDoubleFig.mplstyle')

    fig = plt.figure(figsize=[6.4, 7])
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1])
    gs.update(hspace=0)
    axS = fig.add_subplot(gs[0])  # solid
    axL = fig.add_subplot(gs[1], sharex=axS)  # liquid
    axG = fig.add_subplot(gs[2], sharex=axS)  # gas

    axS.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True, labelbottom=False)
    axL.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True, labelbottom=False)
    axG.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True, labelbottom=True)

    axS.set(
        ylabel=r'$g(r*)$',
        xlim=[0, 4]
    )
    axL.set(
        ylabel=r'$g(r*)$',
    )
    axG.set(
        ylabel=r'$g(r*)$',
        xlabel=r'$r*$'
    )
    for cR in cRs:
        axis = None
        if cR.density == 0.1:
            axis = axG
        elif cR.density == 0.8:
            axis = axL
        elif cR.density == 1.5:
            axis = axS
        else:
            print('No axis selected, check whether the correct files are present (dens0.3, 0.8, 1.3)')
        axis.plot(cR.PCFradii, cR.PCF,
                  marker='', ls='-', lw=1.5, c='k')
    fig.savefig(join(saveLocation, r'RDFs.png'), format='png', dpi=1200)
    fig.show()


def movingaverage(values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma


# def moving_average(a, n=3):
#     ret = np.cumsum(a, dtype=float)
#     ret[n:] = ret[n:] - ret[:-n]
#     return ret[n - 1:] / n


def cumulativeaverage(values):
    cumsum = np.cumsum(values, dtype=float)
    cumaverage = [cumsum[i]/(i+1) for i in range(len(cumsum))]
    return cumaverage


def importCompResults(cR_parent_dir):
    compResults = []
    redTemp_dirs = [f for f in listdir(cR_parent_dir) if not isfile(join(cR_parent_dir, f)) and 'redTemp' in f]

    if len(redTemp_dirs) == 0:
        data_files = [f for f in listdir(cR_parent_dir) if isfile(join(cR_parent_dir, f)) and ('den' in f)]
        for j in range(len(data_files)):
            df = data_files[j]
            data_path = join(cR_parent_dir, df)
            infile = bz2.BZ2File(data_path, 'r')
            cR = pickle.load(infile)
            infile.close()
            compResults.append(cR)

    else:
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
    print(f'---{timedelta(seconds=startTime)}---')


    cR_parent_dir = r'report/figures/excessEnergyPressure'
    saveLocation = cR_parent_dir
    compResults = importCompResults(cR_parent_dir)
    plotCompResults(compResults, saveLocation)

    # cR_parent_dir = r'report/figures/Johnson'
    # saveLocation = cR_parent_dir
    # compResults = importCompResults(cR_parent_dir)
    # plotJohnsoncomparison(compResults, saveLocation)

    # cR_parent_dir = r'report/figures/NIST_comparison/compResult'
    # saveLocation = r'report/figures/NIST_comparison'
    # compResults = importCompResults(cR_parent_dir)
    # plotNISTcomparison(compResults, saveLocation)

    # plotConvergence()

    # BWF_parent_dir = r'report/figures/binWidth/compResult_VaryBWFs'
    # BWFcRs = importCompResults(BWF_parent_dir)
    # plotBinWidthEffect(BWFcRs)

    # cR_parent_dir = r'report\figures\rdfs'
    # saveLocation = cR_parent_dir
    # compResults = importCompResults(cR_parent_dir)
    # plotRDFs(compResults, saveLocation)


    endTime = time() - startTime
    print(f'---{timedelta(seconds=endTime)}---')
    quit()




