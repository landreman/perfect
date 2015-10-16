#!/usr/bin/env python

#roughly simulates T to get constant heat flux.

from perfectInputFile import perfectInput
from perfectProfilesFile import perfectProfiles
from perfectGeometryFile import perfectGeometry

import numpy
import scipy.special
from sys import exit, argv

inputfilename = argv[1]
inputfile = perfectInput(inputfilename)
Npsi = inputfile.Npsi
Ntheta = inputfile.Ntheta
numSpecies = inputfile.numSpecies
psi = inputfile.psi
psiMid = inputfile.psiMid
theta = inputfile.theta
epsil = inputfile.epsil
desiredFWHMInRhoTheta = inputfile.desiredFWHMInRhoTheta
Miller_q = inputfile.Miller_q
desiredU = inputfile.desiredU
dTHatdpsiScalar = inputfile.dTHatdpsiScalar
detaHatdpsiScalar = inputfile.detaHatdpsiScalar
ddpsi_accurate = inputfile.ddpsi_accurate

# Create geometry
# This creates circular concentric flux surfaces from the parameters in the input file. No psi dependence.
psi_dependence=numpy.ones((Npsi,1))
BHat=psi_dependence*(1/(1 + epsil * numpy.cos(theta)))
dBHatdpsi=numpy.zeros((Npsi,Ntheta))
dBHatdtheta=psi_dependence*(epsil*numpy.sin(theta)/((1 + epsil * numpy.cos(theta))**2))
JHat=BHat/Miller_q
IHat=numpy.ones((Npsi))
dIHatdpsi=numpy.zeros((Npsi))

#transpose test
#BHat=numpy.transpose(BHat)
#dBHatdpsi=numpy.transpose(dBHatdpsi)
#dBHatdtheta=numpy.transpose(dBHatdtheta)
#JHat=numpy.transpose(JHat)
#IHat=numpy.transpose(IHat)
#dIHatdpsi=numpy.transpose(dIHatdpsi)



# Output geometry
geometry_outputfile = perfectGeometry(inputfile.geometryFilename)
geometry_outputfile.create_geometry_for_Npsi_Ntheta(Npsi,Ntheta,BHat,dBHatdpsi,dBHatdtheta,JHat,IHat,dIHatdpsi)



exit(0)
