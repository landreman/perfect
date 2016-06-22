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
numFailures += shouldBe("run  1/d(nHat)d(psi)[24,1;;;]", -15.8659, desiredRelativeTolerance)
numFailures += shouldBe("run  1/FSABFlow[0,0;;;]", 37.5427, desiredRelativeTolerance)
numFailures += shouldBe("run  1/FSABFlow[0,1;;;]", -15.4288, desiredRelativeTolerance)
numFailures += shouldBe("run  1/densityPerturbation[2,8,0;;;]", 0.000296172, desiredRelativeTolerance)
numFailures += shouldBe("run  1/densityPerturbation[2,8,1;;;]", 6.39766e-05, desiredRelativeTolerance)
numFailures += shouldBe("run  1/flow[16,7,0;;;]", 31.6209, desiredRelativeTolerance)
numFailures += shouldBe("run  1/flow[16,7,1;;;]", 15.5161, desiredRelativeTolerance)
numFailures += shouldBe("run  1/heatFlux[6,0;;;]", -0.130092, desiredRelativeTolerance)
numFailures += shouldBe("run  1/heatFlux[6,1;;;]", -0.0392926, desiredRelativeTolerance)
numFailures += diffAll(desiredAbsoluteTolerance)

exit(numFailures > 0)
