#!/usr/bin/env python

# This python script checks the perfectOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main PERFECT directory.

execfile('../testsCommon.py')

desiredRelativeTolerance = 0.001
desiredAbsoluteTolerance = 0.001

numFailures = 0

numFailures += shouldBe("run  1/BHat[1,1;;;]", 0.874211, desiredRelativeTolerance)
numFailures += shouldBe("run  1/THat[25,0;;;]", 0.920468, desiredRelativeTolerance)
numFailures += shouldBe("run  1/d(PhiHat)d(psi)[1;;;]", 0.00178617, desiredRelativeTolerance)
numFailures += shouldBe("run  1/d(nHat)d(psi)[24,1;;;]", -0.828644, desiredRelativeTolerance)
numFailures += shouldBe("run  1/FSABFlow[0,0;;;]", 17.4789, desiredRelativeTolerance)
numFailures += shouldBe("run  1/FSABFlow[0,1;;;]", 12.523, desiredRelativeTolerance)
numFailures += shouldBe("run  1/densityPerturbation[2,1,0;;;]", -0.00148854, desiredRelativeTolerance)
numFailures += shouldBe("run  1/densityPerturbation[2,1,1;;;]", -0.000326652, desiredRelativeTolerance)
numFailures += shouldBe("run  1/flow[16,7,0;;;]", 19.1385, desiredRelativeTolerance)
numFailures += shouldBe("run  1/flow[16,7,1;;;]", 15.5835, desiredRelativeTolerance)
numFailures += shouldBe("run  1/heatFlux[6,0;;;]", -1.7565, desiredRelativeTolerance)
numFailures += shouldBe("run  1/heatFlux[6,1;;;]", -0.375692, desiredRelativeTolerance)
numFailures += diffAll(desiredAbsoluteTolerance)

exit(numFailures > 0)
