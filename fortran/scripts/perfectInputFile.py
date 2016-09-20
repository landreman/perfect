import f90nml # You may need to install this module if you do not have it already, e.g. "pip install --user f90nml"
import numpy
import scipy.special
import subprocess
from sys import exit

########################################################
# Class for reading and manipulating PERFECT input files
########################################################
class perfectInput(object):
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


    def groupSetter(group,attr):
        def setter(self, value):
            try:
                #just to see if it exists
                if attr not in self.inputfile[group]:
                    raise KeyError('Attribute ' + attr + 'does not exists in group ' + group +'.')
            except KeyError:
                print "PerfectInputFile: Error: group: " + group + ", attribute: " + attr + " combination does not exist!"
                return None
            self.changevar(group,attr,value)
            self.inputfile[group][attr]=value
        return setter
    
    def groupGetter(group,attr):
        def getter(self):
            try:
                return self.inputfile[group][attr]
            except KeyError:
                print "PerfectInputFile: Error: group: " + group + ", attribute: " + attr + " combination does not exist!"
                return None
        return getter

    programMode = property(groupGetter('flowControl','programMode'), groupSetter('flowControl','programMode'))
    outputFilename = property(groupGetter('flowControl','outputFilename'), groupSetter('flowControl','outputFilename'))
    solveSystem = property(groupGetter('flowControl','solveSystem'), groupSetter('flowControl','solveSystem'))

    profilesScheme = property(groupGetter('physicsParameters','profilesScheme'), groupSetter('physicsParameters','profilesScheme'))
    profilesFilename = property(groupGetter('physicsParameters','profilesFilename'), groupSetter('physicsParameters','profilesFilename'))
    psiAHat = property(groupGetter('physicsParameters','psiAHat'), groupSetter('physicsParameters','psiAHat'))
    leftBoundaryShift = property(groupGetter('physicsParameters','leftBoundaryShift'), groupSetter('physicsParameters','leftBoundaryShift'))
    rightBoundaryShift = property(groupGetter('physicsParameters','rightBoundaryShift'), groupSetter('physicsParameters','rightBoundaryShift'))
    psiMid = property(groupGetter('physicsParameters','psiMid'), groupSetter('physicsParameters','psiMid'))
    desiredU = property(groupGetter('physicsParameters','desiredU'), groupSetter('physicsParameters','desiredU'))
    desiredUNumRuns = property(groupGetter('physicsParameters','desiredUNumRuns'), groupSetter('physicsParameters','desiredUNumRuns'))
    desiredFWHMInRhoTheta = property(groupGetter('physicsParameters','desiredFWHMInRhoTheta'), groupSetter('physicsParameters','desiredFWHMInRhoTheta'))
    Delta =  property(groupGetter('physicsParameters','Delta'), groupSetter('physicsParameters','Delta'))
    omega =  property(groupGetter('physicsParameters','omega'), groupSetter('physicsParameters','omega'))
    nu_r =  property(groupGetter('physicsParameters','nu_r'), groupSetter('physicsParameters','nu_r'))
    dTHatdpsiScalar =  property(groupGetter('physicsParameters','dTHatdpsiScalar'), groupSetter('physicsParameters','dTHatdpsiScalar'))
    detaHatdpsiScalar =  property(groupGetter('physicsParameters','detaHatdpsiScalar'), groupSetter('physicsParameters','detaHatdpsiScalar'))
    makeLocalApproximation =  property(groupGetter('physicsParameters','makeLocalApproximation'), groupSetter('physicsParameters','makeLocalApproximation'))
    includeddpsiTerm =  property(groupGetter('physicsParameters','includeddpsiTerm'), groupSetter('physicsParameters','includeddpsiTerm'))
    
    NpsiPerDiameter =  property(groupGetter('resolutionParameters','NpsiPerDiameter'), groupSetter('resolutionParameters','NpsiPerDiameter'))
    psiDiameter =  property(groupGetter('resolutionParameters','psiDiameter'), groupSetter('resolutionParameters','psiDiameter'))
    widthExtender =  property(groupGetter('resolutionParameters','widthExtender'), groupSetter('resolutionParameters','widthExtender'))
    Nxi = property(groupGetter('resolutionParameters','Nxi'), groupSetter('resolutionParameters','Nxi'))
    Ntheta = property(groupGetter('resolutionParameters','Ntheta'), groupSetter('resolutionParameters','Ntheta'))
    NpsiNumRuns = property(groupGetter('resolutionParameters','NpsiNumRuns'), groupSetter('resolutionParameters','NpsiNumRuns'))
    psiDiameterNumRuns = property(groupGetter('resolutionParameters','psiDiameterNumRuns'), groupSetter('resolutionParameters','psiDiameterNumRuns'))
    widthExtenderNumRuns = property(groupGetter('resolutionParameters','widthExtenderNumRuns'), groupSetter('resolutionParameters','widthExtenderNumRuns'))
    

    charges = property(groupGetter('speciesParameters','charges'), groupSetter('speciesParameters','charges'))
    masses = property(groupGetter('speciesParameters','masses'), groupSetter('speciesParameters','masses'))

    geometryToUse = property(groupGetter('geometryParameters','geometryToUse'), groupSetter('geometryParameters','geometryToUse'))
    geometryFilename = property(groupGetter('geometryParameters','geometryFilename'), groupSetter('geometryParameters','geometryFilename'))
    epsil = property(groupGetter('geometryParameters','epsil'), groupSetter('geometryParameters','epsil'))
    Miller_kappa = property(groupGetter('geometryParameters','Miller_kappa'), groupSetter('geometryParameters','Miller_kappa'))
    Miller_delta = property(groupGetter('geometryParameters','Miller_delta'), groupSetter('geometryParameters','Miller_delta'))
    Miller_s_delta = property(groupGetter('geometryParameters','Miller_s_delta'), groupSetter('geometryParameters','Miller_s_delta'))
    Miller_s_kappa = property(groupGetter('geometryParameters','Miller_s_kappa'), groupSetter('geometryParameters','Miller_s_kappa'))
    Miller_dRdr = property(groupGetter('geometryParameters','Miller_dRdr'), groupSetter('geometryParameters','Miller_dRdr'))
    Miller_q = property(groupGetter('geometryParameters','Miller_q'), groupSetter('geometryParameters','Miller_q'))

    psiGridType = property(groupGetter('otherNumericalParameters','psiGridType'), groupSetter('otherNumericalParameters','psiGridType'))

    useGlobalTermMultiplier = property(groupGetter('otherNumericalParameters','useGlobalTermMultiplier'), groupSetter('otherNumericalParameters','useGlobalTermMultiplier'))

    def read(self,inputfilename):
        self.inputfilename = inputfilename
        self.inputfile = f90nml.read(inputfilename)
        flowControl = self.inputfile["flowControl"]
        physicsParameters = self.inputfile["physicsParameters"]
        speciesParameters = self.inputfile["speciesParameters"]
        resolutionParameters = self.inputfile["resolutionParameters"]
        otherNumericalParameters = self.inputfile["otherNumericalParameters"]
       
        #self.programMode=flowControl["programMode"]

        if (self.NpsiNumRuns>0 or self.psiDiameterNumRuns>0 or self.widthExtenderNumRuns>0) and ((self.programMode == 2) or (self.programMode == 3)):
            print "perfectInputFile: Error: Scans that change Npsi are not supported"
            exit(1)
        
        #self.outputFilename=flowControl["outputFilename"]
        #self.solveSystem=flowControl["solveSystem"]
        #self.profilesScheme = physicsParameters["profilesScheme"]
        if self.profilesScheme == 7:
            if self.profilesFilename == None:
                print "PerfectInputFile: Error: profilesScheme == 7, but no profilesFilename specified."
                exit(1)
        #self.psiAHat = physicsParameters["psiAHat"]
        #self.NpsiPerDiameter = resolutionParameters["NpsiPerDiameter"]
        #self.psiDiameter = resolutionParameters["psiDiameter"]
        #self.widthExtender = resolutionParameters["widthExtender"]
        #self.Nxi = resolutionParameters["Nxi"]
        
        
        #self.charges=speciesParameters["charges"]
        #self.masses=speciesParameters["masses"]

        #self.leftBoundaryShift = physicsParameters["leftBoundaryShift"]
        #self.rightBoundaryShift = physicsParameters["rightBoundaryShift"]
        #self.psiMid = physicsParameters["psiMid"]
        #self.desiredU = physicsParameters["desiredU"]
        if (self.desiredUNumRuns>0) and (self.programMode == 4):
            print "PerfectInputFile: Error: scans in desiredU not supported"
            exit(2)
        #self.desiredFWHMInRhoTheta = physicsParameters["desiredFWHMInRhoTheta"]
        #self.Delta= physicsParameters["Delta"]
        #self.omega= physicsParameters["omega"]
        #self.nu_r= physicsParameters["nu_r"]
        #self.dTHatdpsiScalar = physicsParameters["dTHatdpsiScalar"]
        #self.detaHatdpsiScalar = physicsParameters["detaHatdpsiScalar"]

        self.Npsi = self.computeNpsi(self.NpsiPerDiameter, self.psiDiameter, self.widthExtender, self.leftBoundaryShift, self.rightBoundaryShift)
        #self.makeLocalApproximation=physicsParameters["makeLocalApproximation"]

        # psi grid is uniformly spaced and should always include psiMin and psiMax points (options for psiDerivative Scheme are 1 and 2)
        #self.psiGridType = otherNumericalParameters["psiGridType"]
        if self.psiGridType == None:
            #this will give the Getter something to return,
            #the sed part of the Setter should silently fail
            self.inputfile["otherNumericalParameters"]["psiGridType"] = 0

        if self.useGlobalTermMultiplier == None:
            #this will give the Getter something to return,
            #the sed part of the Setter should silently fail
            self.inputfile["otherNumericalParameters"]["useGlobalTermMultiplier"] = 0
        
            
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
        #self.Ntheta = resolutionParameters["Ntheta"]
        self.theta = numpy.linspace(0.,2.*numpy.pi,self.Ntheta,endpoint=False)
        self.dtheta = self.theta[1] - self.theta[0]
        self.thetaWeights = numpy.full(self.Ntheta,self.dtheta)
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
        #self.geometryToUse = geometryParameters["geometryToUse"]
        if self.geometryToUse == 4:
            if self.geometryFilename == None:
                print "PerfectInputFile: Error: geometryToUse == 4, but no geometryFilename specified."
                exit(1)
                
        #self.epsil = geometryParameters["epsil"]
        #self.Miller_kappa = geometryParameters["Miller_kappa"]
        #self.Miller_delta = geometryParameters["Miller_delta"]
        #self.Miller_s_delta = geometryParameters["Miller_s_delta"]
        #self.Miller_s_kappa = geometryParameters["Miller_s_kappa"]
        #self.Miller_dRdr = geometryParameters["Miller_dRdr"]
        #self.Miller_q = geometryParameters["Miller_q"]
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
        if type(value) == str:
            #strings must be enclosed in "" in namelists
            #may be wise to see if the string contains citation marks...
            if (value.find("'") != -1) or (value.find('"') != -1):
                print "Warning! String to changevar contains a ' or \" character." 
            value = '"' + value + '"'
        if (type(value) == list) or (type(value) == numpy.ndarray):
            # arrays are space seperated
            delimiter=' '
            value_temp = '' 
            for val in value:
                value_temp =  value_temp + str(val) + delimiter
            value = value_temp.rsplit(delimiter,1)[0]
            
        subprocess.call("sed -i -e '/\&"+group+"/,/\&/{ s/^  "+var+" =.*/  "+var+" = "+str(value)+"/I } ' "+self.inputfilename, shell=True)

    def __init__(self,inputfilename=None):
        if inputfilename:
            self.read(inputfilename)
