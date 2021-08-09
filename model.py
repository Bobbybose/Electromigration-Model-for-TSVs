#EM Model for TSVs in 3d-stacked DRAM architectures
#Author: Bobby Bose
#Worked conducted with the University of Kentucky

import math

#Constants
ratioOfCapturedVacancies = 1
ratioOfVacancyVolume = 0.4
atomicVolume = 1.18e-29
tsvVoidThickness = 5e-9
initialDiffusivity = 0.0047
activationEnergy = 1.30e-19
boltzmanConstant = 1.38e-23
temperature = 4.53e2
electronCharge = 1.6e-19
effectiveChargeConstant = 1
barrierResistivity = 3e-6
tsvRadius = 1.15e-6
effectiveVoidRadius = 1.15e-6
atomicConcentration = 1.53e28

IDD0 = 55e-3
numTSVs = 32
timeStep = 5e6
resGainSlope = 7.78
resGainInt = -8.73944
maxParallelism = 16

#Needed values
vacancyConcentration = atomicConcentration*math.exp(-activationEnergy/(boltzmanConstant*temperature))
vacancyDiffusivity = initialDiffusivity*math.exp(-activationEnergy/(boltzmanConstant*temperature))
tsvCrossSectionArea = math.pi*tsvRadius*tsvRadius

#TSV Class
class TSV:
    def __init__(self, parallelismLevel):
        self.parallelism = parallelismLevel
        self.totalCurrent = self.parallelism*IDD0
        self.currentPerTSV = self.totalCurrent/64
        self.currentDensity = self.currentPerTSV/tsvCrossSectionArea
        self.vacancyFlux = vacancyDiffusivity*vacancyConcentration*((electronCharge*effectiveChargeConstant)/(boltzmanConstant*temperature))*barrierResistivity*self.currentDensity
        self.dr = (ratioOfCapturedVacancies*ratioOfVacancyVolume*atomicVolume*effectiveVoidRadius*self.vacancyFlux*timeStep)/tsvVoidThickness


#Contains TSV data for different SA-level parallelism values
parallelismTSVs = []

#Filling in list for different allowable parallelism levels
i = maxParallelism
while i > 1:
   parallelismTSVs.append( TSV(i))
   i = int(i/2)
#parallelismTSVs.append(TSV(1))

resLimits = [2.79, 6.76, 14.7, 30.58]

timeDrops = []

time = 0
lastTime = 0

rVoid = 0
lastRVoid = 0

resGain = 0
lastResGain = 0

stop = 0
level = 0
currentParallelism = maxParallelism

print(maxParallelism, 'SA-level parallelism:')

#Applying EM Model
while currentParallelism > 1:
    time += timeStep
    rVoid += parallelismTSVs[level].dr
    resGain = resGainSlope*rVoid*1e6+resGainInt

    #print(time, rVoid, resGain)
    #print()

    if(resGain > resLimits[level]):
        timeDrops.append(lastTime)
        time = lastTime
        rVoid = lastRVoid
        resGain = lastResGain
        level += 1
        currentParallelism /= 2

        print('Lifetime =', "{:2e}".format(time), 'seconds.')
        print('Or', "{:.2f}".format(time/60/60/24/365.25), 'years.')

        if(currentParallelism == 1):
            break
        else:
            print();
            print(int(currentParallelism), 'SA-level parallelism:')
            continue

    lastTime = time
    lastRVoid = rVoid
    lastResGain = resGain

