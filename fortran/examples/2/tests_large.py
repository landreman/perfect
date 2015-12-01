#!/usr/bin/env python

# This python script checks the perfectOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main PERFECT directory.

execfile('../testsCommon.py')

desiredRelativeTolerance = 0.001
desiredAbsoluteTolerance = 0.0001

numFailures = 0

numFailures += shouldBe("run  1/BHat[1,1;;;]", 0.874211, desiredRelativeTolerance)
numFailures += shouldBe("run  1/THat[25,0;;;]", 0.920468, desiredRelativeTolerance)
numFailures += shouldBe("run  1/d(PhiHat)d(psi)[1;;;]", 0.00178617, desiredRelativeTolerance)
numFailures += shouldBe("run  1/d(nHat)d(psi)[24,1;;;]", -0.827984, desiredRelativeTolerance)
numFailures += shouldBe("run  1/FSABFlow[0,0;;;]", 16.6699, desiredRelativeTolerance)
numFailures += shouldBe("run  1/FSABFlow[0,1;;;]", 11.8503, desiredRelativeTolerance)
numFailures += shouldBe("run  1/densityPerturbation[2,1,0;;;]", -0.00132442, desiredRelativeTolerance)
numFailures += shouldBe("run  1/densityPerturbation[2,1,1;;;]", -0.00032581, desiredRelativeTolerance)
numFailures += shouldBe("run  1/flow[16,7,0;;;]", 15.4467, desiredRelativeTolerance)
numFailures += shouldBe("run  1/flow[16,7,1;;;]", 13.3528, desiredRelativeTolerance)
numFailures += shouldBe("run  1/heatFlux[6,0;;;]", -1.5764, desiredRelativeTolerance)
numFailures += shouldBe("run  1/heatFlux[6,1;;;]", -0.37723, desiredRelativeTolerance)
numFailures += diffAll(desiredAbsoluteTolerance)

exit(numFailures > 0)
