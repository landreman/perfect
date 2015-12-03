import h5py
import numpy
from sys import exit

####################################################
# Class for creating profile input files for PERFECT
####################################################
class perfectProfiles:
    
    def __init__(self,profilesFilename=None):
        if profilesFilename==None:
            print "Error: perfectProfiles must be initialized with a name for the profiles file"
            exit(1)
        self.profilesFile = h5py.File(profilesFilename,'w')

    def __del__(self):
        try:
            self.profilesFile.close()
        except ValueError: # If file is not open, don't need to do anything
            pass

    def close(self):
        self.profilesFile.close()

    def create_profiles_for_Npsi(self, Npsi, numSpecies, PhiHat, dPhiHatdpsi, THats, dTHatdpsis, nHats, dnHatdpsis, etaHats, detaHatdpsis):
        # Make sure inputs are numpy arrays
        PhiHat = numpy.array(PhiHat)
        dPhiHatdpsi = numpy.array(dPhiHatdpsi)
        THats = numpy.array(THats)
        dTHatdpsis = numpy.array(dTHatdpsis)
        nHats = numpy.array(nHats)
        dnHatdpsis = numpy.array(dnHatdpsis)
        etaHats = numpy.array(etaHats)
        detaHatdpsis = numpy.array(detaHatdpsis)

        # Check arrays are the right size
        for var,name in ((PhiHat,"PhiHat"), (dPhiHatdpsi,"dPhiHatdpsi")):
            if var.shape != (Npsi,):
                print "Error: "+name+" has dimensions "+str(var.shape)+" instead of (Npsi) =",str((Npsi,))
                exit(1)
        for var,name in ((THats,"THats"), (dTHatdpsis,"dTHatdpsis"), (nHats,"nHats"), (dnHatdpsis,"dnHatdpsis"), (etaHats,"etaHats"), (detaHatdpsis,"detaHatdpsis")):
            if var.shape != (Npsi,numSpecies):
                print "Error: "+name+" has dimensions "+str(var.shape)+" instead of (Npsi,numSpecies) =",str((Npsi,numSpecies))
                exit(1)

        # Write profiles to HDF5 file
        groupname = "Npsi"+str(Npsi)
        if groupname in self.profilesFile:
            print "Error: profiles already added for Npsi="+str(Npsi)
            exit(1)
        profilesgroup = self.profilesFile.create_group(groupname)

        profilesgroup.create_dataset("PhiHat",data=PhiHat)
        profilesgroup.create_dataset("dPhiHatdpsi",data=dPhiHatdpsi)
        profilesgroup.create_dataset("THats",data=THats)
        profilesgroup.create_dataset("dTHatdpsis",data=dTHatdpsis)
        profilesgroup.create_dataset("nHats",data=nHats)
        profilesgroup.create_dataset("dnHatdpsis",data=dnHatdpsis)
        profilesgroup.create_dataset("etaHats",data=etaHats)
        profilesgroup.create_dataset("detaHatdpsis",data=detaHatdpsis)

    def create_profiles_for_Npsi_Ntheta(self, Npsi, Ntheta, nHatNeutral, dnHatNeutraldpsi):
        groupname = "Npsi"+str(Npsi)+"Ntheta"+str(Ntheta)
        if groupname in self.profilesFile:
            print "Error: profiles already added for Npsi="+str(Npsi)+" Ntheta="+str(Ntheta)
            exit(1)
        profilesgroup = self.profilesFile.create_group(groupname)

        profilesgroup.create_dataset("nHatNeutral",data=nHatNeutral)
        profilesgroup.create_dataset("dnHatNeutraldpsi",data=dnHatNeutraldpsi)

