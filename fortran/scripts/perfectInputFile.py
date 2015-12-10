import f90nml # You may need to install this module if you do not have it already, e.g. "pip install --user f90nml"
import numpy
import scipy.special
import subprocess
from sys import exit

########################################################
# Class for reading and manipulating PERFECT input files
########################################################
class perfectInput:
    """Read and parse an input file for PERFECT

    This class reads the options from a PERFECT input file and creates some derived quantities as well (such as an array of psi values)

    WARNING: The class is currently incomplete (it does not read all options, or calculate all possible derived quantities). Please feel free to add more if they are useful to you ;-)

    Keyword arguments:
    inputfilename -- the name of the PERFECT input file
    """
    def computeNpsi(self,NpsiPerDiameter0, psiDiameter0, widthExtender0, leftBoundaryShift, rightBoundaryShift):
        return int(round(NpsiPerDiameter0 * (psiDiameter0 + 2*widthExtender0 - leftBoundaryShift + rightBoundaryShift))+1)

    def RHat(self,theta):
        return 1. + 1./self.Miller_A*numpy.cos(theta + self.Miller_x*numpy.sin(theta))

    def BPoloidal(self,theta):
        return self.Miller_QQ/(self.Miller_kappa*self.Miller_q)*numpy.sqrt((numpy.sin(theta+self.Miller_x*numpy.sin(theta)) * (1+self.Miller_x*numpy.cos(theta)))**2 + (self.Miller_kappa*numpy.cos(theta))**2) / (self.RHat(theta) * ( numpy.cos(self.Miller_x*numpy.sin(theta)) + self.Miller_dRdr*numpy.cos(theta) + (self.Miller_s_kappa-self.Miller_s_delta*numpy.cos(theta) + (1+self.Miller_s_kappa)*self.Miller_x*numpy.cos(theta)) * numpy.sin(theta) * numpy.sin(theta + self.Miller_x*numpy.sin(theta)) ))

    def computeBs_1D(self,theta):
        if self.geometryToUse==0:
            return 1./(1.+self.epsil*numpy.cos(theta))
        elif self.geometryToUse==1:
            return numpy.sqrt(self.BPoloidal(theta)**2 + 1./self.RHat(theta)**2)
        elif self.geometryToUse==2:
            return 1. + self.epsil*numpy.cos(theta)
        else:
            print "Error! Invalid geometry"
            exit(1)

    def computedBdthetas_1D(self,theta):
        if self.geometryToUse==0:
            return self.epsil * numpy.sin(theta) / (1.+self.epsil*numpy.cos(theta))**2
        elif self.geometryToUse==1:
            print "Error, case not implemented yet"
            
            
            exit(2)
        elif self.geometryToUse==2:
            return -self.epsil * numpy.sin(theta)
        else:
            print "Error! Invalid geometry"
            exit(1)

    def computeOneOverqRbDotGradThetas_1D(self,theta):
        if self.geometryToUse==0:
            return numpy.ones(len(theta))
        elif self.geometryToUse==1:
            bs = computeBs_1D(theta)
            return bs/Miller_q/BDotGradTheta(theta)
        elif self.geometryToUse==2:
            return 1. / (1.+epsil*numpy.cos(theta))
        else:
            print "Error! Invalid geometry"
            exit(1)

    def ddpsi_accurate(self,x):
        if self.Npsi<5:
            print "ddpsi_accurate cannot be calculated, Npsi<5"
            return;
        # uniformDiffMatrices, scheme 12 (5-point centered difference) 
        dxdpsi = numpy.zeros(x.shape)

        # Central part
        dxdpsi[2:-2] = (x[:-4]-8.*x[1:-3]+8.*x[3:-1]-x[4:])/12./self.dpsi

        # Lower boundary
        dxdpsi[:2] = (-25./12.*x[:2]+4.*x[1:3]-3.*x[2:4]+4./3.*x[3:5]-1./4.*x[4:6])/self.dpsi

        # Upper boundary
        dxdpsi[-2:] = (25./12.*x[-2:]-4.*x[-3:-1]+3.*x[-4:-2]-4./3.*x[-5:-3]+1./4.*x[-6:-4])/self.dpsi

        return dxdpsi

    def ddpsi_spectral(self,x):
        # uniformDiffMatrices, scheme 20 (spectral diffenentation matrices)
        print "Spectral diff matrix unimplemented."
        exit(2)
        dBdthetaResolutionMultiplier = 10
        N=dBdthetaResolutionMultiplier*self.Ntheta
        dxdpsi = numpy.zeros(x.shape)


    def read(self,inputfilename):
        self.inputfilename = inputfilename
        self.inputfile = f90nml.read(inputfilename)
        flowControl = self.inputfile["flowControl"]
        physicsParameters = self.inputfile["physicsParameters"]
        speciesParameters = self.inputfile["speciesParameters"]
        resolutionParameters = self.inputfile["resolutionParameters"]
        otherNumericalParameters = self.inputfile["otherNumericalParameters"]
        #if not physicsParameters["profilesScheme"] == 7: # Check that we actually need to create profiles for this input file
        #!!!! To me, this does not seem like a relevant check for a class to read perfect inputs.
        #    print "This input file does not require a profiles file (profilesScheme!=7)."
        #    exit(1)
        if resolutionParameters["NpsiNumRuns"]>0 or resolutionParameters["psiDiameterNumRuns"]>0 or resolutionParameters["widthExtenderNumRuns"]>0:
            print "Scans that change Npsi are not supported"
            exit(1)
        self.outputFilename=flowControl["outputFilename"]
            
        self.profilesScheme = physicsParameters["profilesScheme"]
        self.profilesFilename = physicsParameters["profilesFilename"]
        self.NpsiPerDiameter = resolutionParameters["NpsiPerDiameter"]
        self.psiDiameter = resolutionParameters["psiDiameter"]
        self.widthExtender = resolutionParameters["widthExtender"]

        self.charges=speciesParameters["charges"]
        self.masses=speciesParameters["masses"]

        self.leftBoundaryShift = physicsParameters["leftBoundaryShift"]
        self.rightBoundaryShift = physicsParameters["rightBoundaryShift"]
        self.psiMid = physicsParameters["psiMid"]
        self.desiredU = physicsParameters["desiredU"]
        if physicsParameters["desiredUNumRuns"]>0:
            print "Error, scans in desiredU not supported"
            exit(2)
        self.desiredFWHMInRhoTheta = physicsParameters["desiredFWHMInRhoTheta"]
        self.Delta= physicsParameters["Delta"]
        self.omega= physicsParameters["omega"]
        self.nu_r= physicsParameters["nu_r"]
        self.dTHatdpsiScalar = physicsParameters["dTHatdpsiScalar"]
        self.detaHatdpsiScalar = physicsParameters["detaHatdpsiScalar"]

        self.Npsi = self.computeNpsi(self.NpsiPerDiameter, self.psiDiameter, self.widthExtender, self.leftBoundaryShift, self.rightBoundaryShift)
        self.makeLocalApproximation=physicsParameters["makeLocalApproximation"]

        # psi grid is uniformly spaced and should always include psiMin and psiMax points (options for psiDerivative Scheme are 1 and 2)
        self.psiMin = self.psiMid - self.psiDiameter/2. - self.widthExtender + self.leftBoundaryShift
        self.psiMax = self.psiMid + self.psiDiameter/2. + self.widthExtender + self.rightBoundaryShift
        self.psi = numpy.linspace(self.psiMin, self.psiMax, self.Npsi)
        #We might have a psi grid with only one point
        if self.Npsi>1:
            self.dpsi = self.psi[1]-self.psi[0]
        elif self.makeLocalApproximation==False:
            print "Npsi is less than 2; this only makes sense when makeLocalApproximation is true, but it is false"
            exit(1)
            
        
        # theta grid is uniform and periodic, needs grid point at 0 but not 2pi
        self.Ntheta = resolutionParameters["Ntheta"]
        self.theta = numpy.linspace(0.,2.*numpy.pi,self.Ntheta,endpoint=False)
        dtheta = self.theta[1] - self.theta[0]
        self.thetaWeights = numpy.full(self.Ntheta,dtheta)
        #self.thetaDerivativeScheme = otherNumericalParameters["thetaDerivativeScheme"]
        #if self.thetaDerivativeScheme==0:
        #elif self.thetaDerivative in [1,2]:
        #else:
        #    print "Error, invalid thetaDerivativeScheme"
        #    exit(1)

        try:
            self.numSpecies = len(self.inputfile["speciesParameters"]["charges"])
        except TypeError:
            self.numSpecies = 1

        # Geometrical stuff
        # Should perhaps be moved if we rewrite so that geometry
        # is not defined by inputfile
        geometryParameters = self.inputfile["geometryParameters"]
        self.geometryToUse = geometryParameters["geometryToUse"]
        self.geometryFilename = geometryParameters["geometryFilename"]

        self.epsil = geometryParameters["epsil"]
        self.Miller_kappa = geometryParameters["Miller_kappa"]
        self.Miller_delta = geometryParameters["Miller_delta"]
        self.Miller_s_delta = geometryParameters["Miller_s_delta"]
        self.Miller_s_kappa = geometryParameters["Miller_s_kappa"]
        self.Miller_dRdr = geometryParameters["Miller_dRdr"]
        self.Miller_q = geometryParameters["Miller_q"]
        if self.geometryToUse==1:
            self.Miller_x = numpy.arcsin(self.Miller_delta)
            self.Miller_A = 1./self.epsil
#            NThetaIntegral = 100
#            QQIntegrand = ((1.+self.Miller_s_kappa)*numpy.sin(self.theta + self.Miller_x*numpy.sin(self.theta)) * (1.+self.Miller_x*numpy.cos(self.theta)) * numpy.sin(self.theta) + numpy.cos(self.theta) * (self.Miller_dRdr + numpy.cos(self.theta + self.Miller_x*numpy.sin(self.theta)) - self.Miller_s_delta*numpy.sin(self.theta + self.Miller_x*numpy.sin(self.theta)) * numpy.sin(self.theta))) / self.RHat(self.theta)
#            self.Miller_QQ = self.Miller_kappa / (2.*numpy.pi*self.Miller_A) * QQIntegrand.sum() * 2.*numpy.pi/NThetaIntegral
#        if self.geometryToUse in [0,1,2]:
#            bs_1D = self.computeBs_1D(self.theta)
#            dbdthetas_1D = self.computedBdthetas_1D(self.theta)
#            self.oneOverqRbDotGradThetas_1D = self.computeOneOverqRbDotGradThetas_1D(self.theta)
#            self.BHat = numpy.repeat(bs_1D[numpy.newaxis,:],self.Npsi,axis=0) # Note, using C-ordering instead of Fortran-ordering, so indices are transposed compared to PERFECT, i.e. here the first index is for psi and the second is for theta
#            self.dBHatdtheta = numpy.repeat(dbdthetas_1D[numpy.newaxis,:],self.Npsi,axis=0)
#            self.JHat = numpy.repeat(bs_1D[numpy.newaxis,:] / self.oneOverqRbDotGradThetas_1D / self.Miller_q, self.Npsi, axis=0)
#            self.dBHatdpsi = 0.
#            self.IHat = 1.
#            self.dIHatdpsi = 0.
#        elif self.geometryToUse==3:
#            print "EFIT interface not implemented yet"
#        elif self.geometryToUse==4:
#            #May want to calculate things here, or not.
#            print "How would the input file parser calculate this??"
#        else:
#            print "Error! Invalid setting for geometryToUse"
#            exit(1)
#
#        self.VPrimeHat = (self.thetaWeights[numpy.newaxis,:] / self.JHat).sum(axis=1)
#        self.FSABHat2 = ( self.thetaWeights[numpy.newaxis,:] * self.BHat**2 / self.JHat / self.VPrimeHat[:,numpy.newaxis] ).sum(axis=1)
#        self.typicalB = numpy.sqrt(self.FSABHat2)

    #def setvar(self,group,var,value):
    #    print group,var,value
    #    print self.inputfile[group][var]
    #    self.inputfile[group][var] = value
    #    print self.inputfile[group][var]
    #    print ""

    #def write(self):
    #    print self.inputfilename
    #    self.inputfile.write(self.inputfilename,force=True)

    def changevar(self,group,var,value):
        # Warning: this command will fail silently if the pattern is not found. Sorry about that.
        # Warning: case insensitive
        subprocess.call("sed -i -e '/\&"+group+"/,/\&/{ s/^  "+var+" =.*/  "+var+" = "+str(value)+"/I } ' "+self.inputfilename, shell=True)

    def __init__(self,inputfilename=None):
        if inputfilename:
            self.read(inputfilename)


