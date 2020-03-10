import numpy as np
import numexpr as ne

def getDistanceSquared(testLoc, partLocs, containerLength, dim):
    diffs = np.absolute(ne.evaluate('partLocs - testLoc'))
    diffs_pbcs = np.zeros_like(diffs)

    for d in range(dim):
        # diffs_pbcs[:, d] = np.where(diffs[:, d] < 0.5 * containerLength, diffs[:, d], containerLength - diffs[:, d])
        diffs_indexed = diffs[:, d]
        diffs_pbcs[:, d] = ne.evaluate('where(diffs_indexed < 0.5 * containerLength, diffs_indexed, containerLength - diffs_indexed)')

    diffsSquared = np.sum(np.square(diffs_pbcs), axis=1)
    return diffsSquared


def redLJ(redSepSq):
    vdW = ne.evaluate('(1/redSepSq)**3')
    PEP = ne.evaluate('vdW**2')
    LJ_potn = ne.evaluate('4*(PEP - vdW)')
    return LJ_potn


def diff_redLJ(redSeps):
    PEP = ne.evaluate('-12 * (1/redSeps)**13')
    vdW = ne.evaluate('+6 * (1/redSeps)**7')
    LJ_potn = ne.evaluate('4*(PEP + vdW)')
    return LJ_potn
