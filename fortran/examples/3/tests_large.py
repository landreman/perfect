#!/usr/bin/env python

# This python script checks the perfectOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main PERFECT directory.

execfile('../testsCommon.py')

desiredRelativeTolerance = 0.001
desiredAbsoluteTolerance = 0.0001

numFailures = 0

numFailures += shouldBe("run  1/BHat[1,1;;;]", 0.833779, desiredRelativeTolerance)
numFailures += shouldBe("run  1/THat[25,0;;;]", 0.896523, desiredRelativeTolerance)
numFailures += shouldBe("run  1/d(PhiHat)d(psi)[25;;;]", 16.1246, desiredRelativeTolerance)
numFailures += shouldBe("run  1/d(nHat)d(psi)[24,1;;;]", -4.87363, desiredRelativeTolerance)
numFailures += shouldBe("run  1/FSABFlow[0,0;;;]", 37.5427, desiredRelativeTolerance)
numFailures += shouldBe("run  1/FSABFlow[0,1;;;]", -15.4288, desiredRelativeTolerance)
numFailures += shouldBe("run  1/densityPerturbation[2,8,0;;;]", 0.000362008, desiredRelativeTolerance)
numFailures += shouldBe("run  1/densityPerturbation[2,8,1;;;]", 6.39766e-05, desiredRelativeTolerance)
numFailures += shouldBe("run  1/flow[16,7,0;;;]", 33.6952, desiredRelativeTolerance)
numFailures += shouldBe("run  1/flow[16,7,1;;;]", 16.4715, desiredRelativeTolerance)
numFailures += shouldBe("run  1/poloidalFlow[7,2,0;;;]", 1.21959, desiredRelativeTolerance)
numFailures += shouldBe("run  1/poloidalFlow[24,8,1;;;]", -36.5724, desiredRelativeTolerance)
numFailures += shouldBe("run  1/toroidalFlow[0,9,0;;;]", 31.5534, desiredRelativeTolerance)
numFailures += shouldBe("run  1/toroidalFlow[25,8,1;;;]", -236.895, desiredRelativeTolerance)
numFailures += shouldBe("run  1/heatFlux[6,0;;;]", -0.123687, desiredRelativeTolerance)
numFailures += shouldBe("run  1/heatFlux[6,1;;;]", -0.0393983, desiredRelativeTolerance)
numFailures += diffAll(desiredAbsoluteTolerance)

exit(numFailures > 0)
