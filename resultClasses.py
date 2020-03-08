class simResult:
    def __init__(self, numParts, sigma, wellDepth, redTemp, numEvolveMoves, displacement, mode, density,
                 partRad, containerLength, configs):
        self.numParts = numParts
        self.sigma = sigma
        self.wellDepth = wellDepth
        self.redTemp = redTemp
        self.numEvolveMoves = numEvolveMoves
        self.displacement = displacement
        self.mode = mode
        self.density = density
        self.partRad = partRad
        self.containerLength = containerLength
        self.configs = configs


class compResult:
    def __init__(self, numParts, sigma, wellDepth, redTemp, numEvolveMoves, displacement, mode, density,
                 partRad, containerLength, PCF, PCFstdev, PCFradii):
        self.numParts = numParts
        self.sigma = sigma
        self.wellDepth = wellDepth
        self.redTemp = redTemp
        self.numEvolveMoves = numEvolveMoves
        self.displacement = displacement
        self.mode = mode
        self.density = density
        self.partRad = partRad
        self.containerLength = containerLength

        self.PCF = PCF
        self.PCFstdev = PCFstdev
        self.PCFradii = PCFradii

        self.intEn_perPart = 0
        self.pressure_minRhokT = 0
