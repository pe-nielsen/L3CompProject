from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from time import time
from datetime import timedelta
import pickle
import bz2
from pathlib import Path
from multiprocessing import Pool
from resultClasses import simResult
from commonFunctions import getDistanceSquared, redLJ
from initial_parameters import simResultDir
import numpy as np
import numexpr as ne
np.random.seed(1)


class Container:
    # the LJ particles exist in the container
    def __init__(self, numParts, density, dim, sigma, redTemp, mode='LJ'):
        self.numParts = numParts  # number of particles in the system
        self.density = density  # density of particles
        self.dim = dim  # dimension of system (only dim=3 is supported)
        self.sigma = sigma  # LJ(sigma) = 0 -> taken to be the radius of the particles
        self.redTemp = redTemp  # the unit-free reduced temperature of the system = kT/epsilon
        self.mode = mode  # treat system as 'Hard' Sphere or under 'LJ' potential

        if mode == 'LJ':
            self.partRad = 1
            self.length = ((numParts*np.pi)/(6*density))**(1/3)  # container length in units of sigma

        # initialising some parameters
        self.partLocs = np.zeros(shape=(numParts, dim))
        self.displacement = 0
        self.validMoves_i = 0
        self.validMoves = 0

    def buildLattice(self):

        startPos = np.array([self.partRad]*3)  # x, y, z
        numIterationsPerDim = int(np.ceil((self.numParts/4)**(1/3)))
        lattice_parameter = np.ceil(self.length/((self.numParts/4)**(1/3)))

        config = np.zeros(shape=(2*self.numParts, self.dim))
        for i in range(numIterationsPerDim):
            for j in range(numIterationsPerDim):
                for k in range(numIterationsPerDim):
                    partIndex0 = 4*(numIterationsPerDim**2*i + numIterationsPerDim*j + k)
                    config[partIndex0, :] = lattice_parameter*np.array([i, j, k])
                    config[partIndex0 + 1, :] = lattice_parameter*np.array([i + 0.5, j + 0.5, k])
                    config[partIndex0 + 2, :] = lattice_parameter*np.array([i + 0.5, j, k + 0.5])
                    config[partIndex0 + 3, :] = lattice_parameter*np.array([i, j + 0.5, k + 0.5])
        self.partLocs = config[:self.numParts, :]
        # print(f'partLocs:\n{self.partLocs}')

    def moveParticles(self, moves):  # mode = 'Hard' for Hard Sphere or 'LJ' for Lennard-Jones
        # todo: Alder & Wainwright apply the moves randomly rather than sequentially

        self.validMoves_i = 0

        # todo: consider instead computing the total LJ of the system and storing to track the difference
        numIters = moves.shape[1]
        for m in range(numIters):
            partNum = m % self.numParts
            original_loc = np.copy(self.partLocs[partNum, :])  # (x, y, z)

            displ = moves[:, m]  # (dx, dy, dz)

            if self.mode == 'LJ':
                # ### reduced values ### - distances are already in reduced values
                diffsSq_i = getDistanceSquared(self.partLocs[partNum, :], self.partLocs, self.length, self.dim)  # getting the original distances
                self.partLocs[partNum, :] = (original_loc + displ) % self.length  # setting the new location
                diffsSq_f = getDistanceSquared(self.partLocs[partNum, :], self.partLocs, self.length, self.dim)  # getting the new distances

                redLJ_i = redLJ(diffsSq_i)  # computing the reduced potentials for each config
                redLJ_f = redLJ(diffsSq_f)

                red_sumLJ_i = np.nansum(redLJ_i)
                red_sumLJ_f = np.nansum(redLJ_f)
                red_dLJ = red_sumLJ_f - red_sumLJ_i

                if red_dLJ > 0:
                    # accept with some probability
                    rnd_criteria = np.random.rand()  # random number between 0 and 1 to compare Boltzmann factor to
                    redBoltzmannFactor = np.exp(-red_dLJ/self.redTemp)

                    if redBoltzmannFactor > rnd_criteria:
                        self.validMoves_i += 1
                        self.validMoves += 1
                    else:
                        self.partLocs[partNum, :] = original_loc  # revert the move if B < rnd
                else:
                    self.validMoves_i += 1
                    self.validMoves += 1

    def showState(self, plotTitle):
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.set(
            xlabel='x',
            ylabel='y',
            zlabel='z',
            title=plotTitle,
            xlim=(0, self.length),
            ylim=(0, self.length),
            zlim=(0, self.length),
        )
        xs = self.partLocs[:, 0]  # picking out the x, y, z values for the locations for plotting
        ys = self.partLocs[:, 1]
        zs = self.partLocs[:, 2]

        s = (2 * self.partRad) ** 4

        ax.scatter(xs, ys, zs, s=1,
                   color='grey')  # todo: verify the particle size is correct - def wrong as every point exactly the same size
        fig.show()


def genMoves_rndV(numMoves, rndBound):
    moves = 2*rndBound*(np.random.rand(3, int(numMoves))-0.5)  # random within rndBound centred on 0
    return moves


def initialiseSystem(numParts, density, dim, sigma, redTemp, mode):
    container = Container(numParts, density, dim, sigma, redTemp, mode)
    container.buildLattice()
    # container.showState(plotTitle=f'Initial configuration;\n'
    #                               f'density={density}')
    return container


def equilibrateSystem(container, numMoves, numEquilIter, initDisplacement):

    movesPerIter = numMoves/numEquilIter
    displacement = initDisplacement

    for i in range(numEquilIter):
        container.validMoves_i = 0
        moves = genMoves_rndV(movesPerIter, displacement)
        container.moveParticles(moves)
        acceptance_rate = container.validMoves_i / movesPerIter

        if acceptance_rate < 0.3:
            displacement *= 0.95
        elif acceptance_rate > 0.5:
            displacement *= 1.05

    container.displacement = displacement
    # container.showState(plotTitle=f'2. Equilibrated Config.;\n'
    #                               f'n_m={numMoves:.2E}, final displ.={displacement:.2E}, mode={container.mode}')
    print(f'Equilibration completed:\n'
          f'    final displacement: {container.displacement},\n'
          f'    acceptance: {container.validMoves/numMoves},\n'
          f'    redTemp: {container.redTemp},\n'
          f'    density: {container.density}\n')


def evolveSystem(container, numMoves, configSampleRate):
    container.validMoves = 0
    numPcfSamples = int(numMoves/configSampleRate)
    displacement = container.displacement
    configEquil = container.partLocs

    # for configs: [particle number, xyz, sample number]
    configs = np.zeros(shape=(configEquil.shape[0], configEquil.shape[1], numPcfSamples))  # all sample configurations

    for i in range(numPcfSamples):
        container.validMoves_i = 0  # setting current valid moves to 0 for optimising displacement
        moves = genMoves_rndV(configSampleRate, displacement)
        container.moveParticles(moves)

        configs[:, :, i] = container.partLocs

    print(f'\nSystem evolved:\n'
          f'    num moves = {numMoves:.2E}\n'
          f'    density = {container.density}\n'
          f'    total acceptance = {container.validMoves/numMoves}\n'
          f'    final displacement = {displacement}\n'
          f'    mode = {container.mode}\n')
    # container.showState(plotTitle=f'3. Post-Evolution Configuration;\n'
    #                               f'n_m={numMoves:.2E}, mode={container.mode}')#, displ.={displacement}')
    return configs


def runSim(tempDens):
    '''create single container with given redTemp and density, evolve it and save the configuration at an interval'''

    redTemp = tempDens[0]
    density = tempDens[1]

    numParts = 10**3
    dim = 3

    # Lennard-Jones Parameters: - http://stp.clarku.edu/simulations/lj/index.html
    sigma = 3.4e-10  # r_min = 2^(1/6)*sigma
    well_depth = 1.65e-21  # J
    Boltzmann_const = 1.38064852e-23  # m^2.kg.s^-2.K-1
    # system_temp = 273.15 + 55  # the initial argon test - Wood and Parker
    # redTemp = Boltzmann_const * system_temp / well_depth
    mode = 'LJ'

    numEquilMoves = int(1e6)  # must be integer
    # numEquilMoves = int(1e4)  # must be integer
    numEquilIter = 1000  # number of iterations of moves for equil - for varying displacement
    initDisplacement = 0.1  # todo: compare displacement and acceptance rate with Ben's work

    numEvolveMoves = 2e6
    configSampleRate = 1000

    displacement = initDisplacement

    print(f'current density: {density}')
    print(f'current redTemp: {redTemp}')

    container = initialiseSystem(numParts, density, dim, sigma, redTemp, mode)
    equilibrateSystem(container, numEquilMoves, numEquilIter, displacement)
    configs = evolveSystem(container, numEvolveMoves, configSampleRate)

    sR = simResult(container.numParts, container.sigma, well_depth, container.redTemp,
                   numEvolveMoves, container.displacement, container.mode, container.density,
                   container.partRad, container.length, configs)

    # Saving the simResult
    # parent_dir = r'C:\Users\Peter Nielsen\Documents\l3\comp_proj\simulationResults\testRun' + f'\\redTemp{redTemp*100}'
    # parent_dir = join(r'C:\Users\splb68\comp_proj\simulationResults\medRun_isotherms', f'redTemp{redTemp}')
    parent_dir = join(simResultDir, f'redTemp{redTemp}')
    Path(parent_dir).mkdir(parents=True, exist_ok=True)
    filename = f'den{sR.density}'
    filepath = join(parent_dir, filename)
    # outfile = open(filepath, 'wb')
    outfile = bz2.BZ2File(filepath, 'w')
    toPickle = sR
    pickle.dump(toPickle, outfile)
    outfile.close()

    print(f'''\n\n
        Configurations saved for:
        density: {sR.density},
        reduced temperature: {sR.redTemp},
        at time: ---{timedelta(seconds=(time()))}---''')


if __name__ == '__main__':
    startTime = time()
    print(f'---{timedelta(seconds=startTime)}---')



    # densities = [0.2, 0.4, 0.5, 0.6, 0.7, 0.8]
    # redTemps = [1, 5, 10]
    # densities = [0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8]
    # redTemps = [0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    densities = [1, 2, 3, 4, 5, 6]
    redTemps = [1, 11]

    tempDenss = []
    for rT in redTemps:
        for d in densities:
            tempDens = [rT, d]
            tempDenss.append(tempDens)

    numPools = 5  # number of cores used for calculations

    with Pool(numPools) as p:
        p.map(runSim, tempDenss)

    endTime = time() - startTime
    print(f'---{timedelta(seconds=endTime)}---')
    quit()
