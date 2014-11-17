#!/usr/bin/env python

# This python script checks the sfincsOutput.h5 file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main SFINCS directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

numFailures += shouldBe("run  1/BHat[1,1;;;]", 0.874211, desiredTolerance)
numFailures += shouldBe("run  1/THat[25,0;;;]", 0.920468, desiredTolerance)
numFailures += shouldBe("run  1/d(PhiHat)d(psi)[1;;;]", 0.00178617, desiredTolerance)
numFailures += shouldBe("run  1/d(nHat)d(psi)[24,1;;;]", -0.828644, desiredTolerance)
numFailures += shouldBe("run  1/FSABFlow[0,0;;;]", 17.417, desiredTolerance)
numFailures += shouldBe("run  1/FSABFlow[0,1;;;]", 12.4774, desiredTolerance)
numFailures += shouldBe("run  1/densityPerturbation[2,1,0;;;]", -0.0014995, desiredTolerance)
numFailures += shouldBe("run  1/densityPerturbation[2,1,1;;;]", -0.000326344, desiredTolerance)
numFailures += shouldBe("run  1/flow[16,7,0;;;]", 19.293, desiredTolerance)
numFailures += shouldBe("run  1/flow[16,7,1;;;]", 15.6281, desiredTolerance)
numFailures += shouldBe("run  1/heatFlux[6,0;;;]", -1.75734, desiredTolerance)
numFailures += shouldBe("run  1/heatFlux[6,1;;;]", -0.376159, desiredTolerance)

exit(numFailures > 0)
