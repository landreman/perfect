import h5py
import numpy
from sys import exit,argv


#for if main:
from perfectInputFile import perfectInput


def create_globalMultiplier_of_Npsi(globalMultiplierFilename, Npsi, globalMultiplier):

    # function on the physical psiN grid that multiplies the global terms in the DKE
        
    if globalMultiplierFilename==None:
        print "Error: globalMultiplierFilename must be a filename"
        exit(1)
    globalMultiplierFile = h5py.File(globalMultiplierFilename,'w')
    groupname = "Npsi"+str(Npsi)
    if groupname in globalMultiplierFile:
        print "Error: profiles already added for Npsi="+str(Npsi)
        exit(1)
    profilesgroup = globalMultiplierFile.create_group(groupname)


    globalMultiplier = numpy.array(globalMultiplier)
    if globalMultiplier.shape != (Npsi,):
        print "Error: globalTermMultiplier has dimensions " + str(globalMultiplier.shape) + " instead of (Npsi,) = " + str((Npsi,))  
    profilesgroup.create_dataset("globalTermMultiplier",data=globalMultiplier)
    
        
if __name__=="__main__":
    inputfilename = argv[1]
    inputfile = perfectInput(inputfilename)
    Npsi = inputfile.Npsi
    globalMultiplier= 1
    globalMultiplier=[1]*Npsi
    create_globalMultiplier_for_Npsi("globalMultiplier.h5",Npsi,globalMultiplier)
    
