import h5py
import numpy
from sys import exit

####################################################
# Class for creating geometry input files for PERFECT
####################################################
class perfectGeometry:
    
    def __init__(self,geometryFilename=None,append=False):
        if geometryFilename==None:
            print "Error: perfectGeometry must be initialized with a name for the geometry file"
            exit(1)
        if append:
            self.geometryFile = h5py.File(geometryFilename,'a')
        else:
            self.geometryFile = h5py.File(geometryFilename,'w')

    def __del__(self):
        try:
            self.geometryFile.close()
        except ValueError:
            # presumably file is not open
            pass

    def close(self):
        self.geometryFile.close()

    def create_geometry_for_Npsi_Ntheta(self, Npsi, Ntheta, psiMin, psiMax, BHat, dBHatdpsi, dBHatdtheta, RHat, JHat, IHat, dIHatdpsi,replace=False):
        # Make sure inputs are numpy arrays
        BHat = numpy.array(BHat)
        dBHatdpsi = numpy.array(dBHatdpsi)
        dBHatdtheta = numpy.array(dBHatdtheta)
        RHat = numpy.array(RHat)
        JHat = numpy.array(JHat)
        IHat = numpy.array(IHat)
        dIHatdpsi = numpy.array(dIHatdpsi)

        # Check arrays are the right size
        for var,name in [[BHat,"BHat"], [dBHatdpsi,"dBHatdpsi"], [dBHatdtheta,"dBHatdtheta"], [RHat,"RHat"],[JHat,"JHat"]]:
            if var.shape != (Npsi,Ntheta):
                print "Error: "+name+" has dimensions "+str(var.shape)+" instead of (Npsi,Ntheta) =",str((Npsi,Ntheta))
                exit(1)
        for var,name in [[IHat,"IHat"],[dIHatdpsi,"dIHatdpsi"]]:
            if var.shape != (Npsi,):
                print "Error: "+name+" has dimensions "+str(var.shape)+" instead of (Npsi) =",str((Npsi,))
                exit(1)

        # Write geometry to HDF5 file
        groupname = "Npsi"+str(Npsi)+"Ntheta"+str(Ntheta)
        if groupname in self.geometryFile:
            if replace:
                del self.geometryFile[groupname]
            else:
                print "Error: profiles already added for Npsi="+str(Npsi)+", Ntheta="+str(Ntheta)
                exit(1)
        geometrygroup = self.geometryFile.create_group(groupname)

        geometrygroup.create_dataset("psiMin",data=psiMin) # psiMin and psiMax are used to check that the geometryFile is consistent with the input.namelist file when PERFECT is run
        geometrygroup.create_dataset("psiMax",data=psiMax)
        geometrygroup.create_dataset("BHat",data=BHat)
        geometrygroup.create_dataset("dBHatdpsi",data=dBHatdpsi)
        geometrygroup.create_dataset("dBHatdtheta",data=dBHatdtheta)
        geometrygroup.create_dataset("RHat",data=RHat)
        geometrygroup.create_dataset("JHat",data=JHat)
        geometrygroup.create_dataset("IHat",data=IHat)
        geometrygroup.create_dataset("dIHatdpsi",data=dIHatdpsi)

    def add_field(self,name,var):
        if name in self.geometryFile:
            if (not isinstance(var,numpy.ndarray)) or var.shape == self.geometryFile[name].shape:
                self.geometryFile[name][...] = var
            else:
                del self.geometryFile[name]
                self.geometryFile.create_dataset(name,data=var)
        else:
            self.geometryFile.create_dataset(name,data=var)
