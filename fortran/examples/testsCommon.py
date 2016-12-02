def readHDF5(variableName):
    # First, check whether perfectOutput.h5 exists:
    
    import os.path
    
    if not os.path.isfile("perfectOutput.h5"):
        print "Error! The file perfectOutput.h5 has not been created."
        exit(1)
        
    # Try executing h5dump with no arguments, just to see if h5dump is installed:
        
    import subprocess
        
    try:
        p = subprocess.Popen("h5dump", stdout=subprocess.PIPE)
    except:
        print "Error! Unable to execute h5dump."
        raise

    # If we made it this far, then the h5dump command is indeed accessible to the shell. Next call h5dump with arguments:

    try:
        p = subprocess.Popen(["h5dump", "-y", "-d", "/"+variableName, "perfectOutput.h5"], stdout=subprocess.PIPE)
    except:
        print "Error! Unable to read perfectOutput.h5 using h5dump. It is likely that perfectOutput.h5 is corrupted."
        raise

    (output, err) = p.communicate()

    # Pick out the lines with the variable contents.
    temp = output.splitlines()
    try:
        # If the correct [indices;;;] variableName was used, the desired value should be on line 10:
        return float(temp[10])
    except:
        print "Error! Unable to convert line 10 of the h5dump output to a single float."
        print "This may occur if the requested variableName lacks the proper [index1,index2,...;;;] indices."
        print "Here is the text that h5dump returned:"
        print output
        raise

def shouldBe(variableName, trueValue, relativeTolerance):
    latestValue = readHDF5(variableName)
    relativeDifference = abs((latestValue - trueValue) / trueValue)
    if relativeDifference > relativeTolerance:
        print "*** TEST FAILED!!  Variable "+variableName+" should be close to "+str(trueValue)+", but it is instead "+str(latestValue)
        print "Actual / correct = ",latestValue/trueValue
        return 1
    else:
        print "    Test passed:   Variable "+variableName+" should be close to "+str(trueValue)+", and it came out to be "+str(latestValue)+", which is within tolerance."
        print "Actual / correct = ",latestValue/trueValue
        return 0

def diffAll(absoluteTolerance):

    import h5py
    import numpy

    newFile = h5py.File("perfectOutput.h5")
    origFile = h5py.File("expectedResults_perfectOutput.h5")

    failureCount = 0

    # Differences in runtime and commit do not matter
    ignoredVariables = ["elapsed time (s)","gitCommit"]
    # these are intermediate variables used to calculate poloidal and toroidal flows
    # they will be practically zero and noisy in the local limit,
    # so we ignore them and only report errors if they actually affect
    # the resulting flows
    ignoredVariables += ["pPerpTermInVp", "pPerpTermInVpBeforePsiDerivative"]

    # these are control variables from older implementation of sources
    ignoredVariables += ["noChargeSource","noChargeSourceOption","sourcePoloidalVariation"]
    for varname in origFile["run  1"]:
        if varname in ignoredVariables:
            continue
        newData = newFile["run  1"][varname][...]
        origData = origFile["run  1"][varname][...]

        absoluteDifference = numpy.abs((newData-origData).max())
        if absoluteDifference > absoluteTolerance:
            print "*** TEST FAILED!!  Variable "+varname+" has a maximum difference of "+str(absoluteDifference)+" ***"
            failureCount += 1

    if failureCount == 0:
        print "    Test passed:   All variables are within "+str(absoluteTolerance)+" absolute difference"
    return failureCount
