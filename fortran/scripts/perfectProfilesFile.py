import h5py
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
        self.profilesFile.close()

    def create_profiles_for_Npsi(self, Npsi, PhiHat, dPhiHatdpsi, THats, dTHatdpsis, nHats, dnHatdpsis, etaHats, detaHatdpsis):
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

