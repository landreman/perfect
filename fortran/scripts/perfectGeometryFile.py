import h5py
from sys import exit

####################################################
# Class for creating geometry input files for PERFECT
####################################################
class perfectGeometry:
    
    def __init__(self,geometryFilename=None):
        if geometryFilename==None:
            print "Error: perfectGeometry must be initialized with a name for the geometry file"
            exit(1)
        self.geometryFile = h5py.File(geometryFilename,'w')

    def __del__(self):
        self.geometryFile.close()

    def create_geometry_for_Npsi_Ntheta(self, Npsi, Ntheta,BHat,dBHatdpsi,dBHatdtheta,JHat,IHat,dIHatdpsi):
        # Write geometry to HDF5 file
        groupname = "Npsi"+str(Npsi)+"Ntheta"+str(Ntheta)
        if groupname in self.geometryFile:
            print "Error: profiles already added for Npsi="+str(Npsi)+", Ntheta="+str(Ntheta)
            exit(1)
        geometrygroup = self.geometryFile.create_group(groupname)

        geometrygroup.create_dataset("BHat",data=BHat)
        geometrygroup.create_dataset("dBHatdpsi",data=dBHatdpsi)
        geometrygroup.create_dataset("dBHatdtheta",data=dBHatdtheta)
        geometrygroup.create_dataset("JHat",data=JHat)
        geometrygroup.create_dataset("IHat",data=IHat)
        geometrygroup.create_dataset("dIHatdpsi",data=dIHatdpsi)


