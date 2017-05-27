import h5py
import numpy
from sys import exit,argv


#for if main:
from perfectInputFile import perfectInput


def create_psiAHat_of_Npsi(psiAHatFilename, Npsi, psiAHatArray,psiArray=None):

    # psiAHatArray: relates derivatives on uniform (psiN) and non-uniform (psi) grid
    # psiArray: optional array of f^{-1}(psiN[i]), psiN=f(psi). Needed to intepret simulation outputs on nonuniform grid
    
    if psiAHatFilename is None:
        print "Error: psiAHatFilename must be a filename"
        exit(1)
    psiAHatFile = h5py.File(psiAHatFilename,'w')
    groupname = "Npsi"+str(Npsi)
    if groupname in psiAHatFile:
        print "Error: profiles already added for Npsi="+str(Npsi)
        exit(1)
    profilesgroup = psiAHatFile.create_group(groupname)


    psiAHatArray = numpy.array(psiAHatArray)
    if psiAHatArray.shape != (Npsi,):
        print "Error: psiAHatArray has dimensions " + str(psiAHatArray.shape) + " instead of (Npsi,) = " + str((Npsi,))  
    profilesgroup.create_dataset("psiAHatArray",data=psiAHatArray)
    
    if psiArray is not None:
        psiArray = numpy.array(psiArray)
        if psiArray.shape != (Npsi,):
            print "Error: psiArray has dimensions " + str(psiArray.shape) + " instead of (Npsi,) = " + str((Npsi,))  
        profilesgroup.create_dataset("psiArray",data=psiArray)
    
        
if __name__=="__main__":
    inputfilename = argv[1]
    inputfile = perfectInput(inputfilename)
    Npsi = inputfile.Npsi
    psiAHat= inputfile.psiAHat
    psiAHatArray=[psiAHat]*Npsi
    create_psiAHat_for_Npsi("psiAHat.h5",Npsi,psiAHatArray)
