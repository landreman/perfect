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
        return 1
    else:
        print "    Test passed:   Variable "+variableName+" should be close to "+str(trueValue)+", and it came out to be "+str(latestValue)+", which is within tolerance."
        return 0
  
