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

# Create profiles

#PhiHat = numpy.zeros(Npsi)
#dPhiHatdpsi = numpy.zeros(Npsi)
THats = numpy.zeros((numSpecies,Npsi))
dTHatdpsis = numpy.zeros((numSpecies,Npsi))
nHats = numpy.zeros((numSpecies,Npsi))
dnHatdpsis = numpy.zeros((numSpecies,Npsi))
etaHats = numpy.zeros((numSpecies,Npsi))
detaHatdpsis = numpy.zeros((numSpecies,Npsi))

delta = (0.0006+0.)*((epsil/(0.01+0.))**2)/(desiredFWHMInRhoTheta/(0.24+0.))/30.
omega = delta
psiAHat = epsil**2/Miller_q
#FSABHat2Fine = inputfile.FSABHat2[0]

# Write delta, omega, psiAHat back to the input file (to match their being overridden in profiles.F90)
#inputfile.setvar("physicsParameters","delta",delta)
#inputfile.setvar("physicsParameters","omega",omega)
#inputfile.setvar("physicsParameters","psiAHat",psiAHat)
#inputfile.write()
inputfile.changevar("physicsParameters","delta",delta)
inputfile.changevar("physicsParameters","omega",omega)
inputfile.changevar("physicsParameters","psiAHat",psiAHat)

ss = 25
rampAmplitude = desiredU * psiAHat / omega

PhiHat = rampAmplitude*scipy.special.erf((psi-psiMid)*ss)*numpy.sqrt(numpy.pi)/2./ss
dPhiHatdpsi = rampAmplitude*numpy.exp(-((psi-psiMid)*ss)**2)

for species in range(numSpecies):
    #THats=[0]*numSpecies
    #dTHatdpsis=[0]*numSpecies
 
    THats[species] = 1. + dTHatdpsiScalar * (psi-psiMid)
    dTHatdpsis[species] = dTHatdpsiScalar

    
    
    etaHats[species] = 1. + (psi-psiMid)*detaHatdpsiScalar
    detaHatdpsis[species] = ddpsi_accurate(etaHats[species])

    nHats[species] = etaHats[species] * numpy.exp(-2.*omega/delta*PhiHat/THats[species])
    dnHatdpsis[species] = ddpsi_accurate(nHats[species])

    
# Output profiles
profile_outputfile = perfectProfiles(inputfile.profilesFilename)
profile_outputfile.create_profiles_for_Npsi(Npsi,PhiHat,dPhiHatdpsi,THats,dTHatdpsis,nHats,dnHatdpsis,etaHats,detaHatdpsis)

# Create geometry
# This creates circular concentric flux surfaces from the parameters in the input file. No psi dependence.
psi_dependence=numpy.ones((Npsi,1))
BHat=psi_dependence*(1/(1 + inputfile.epsil * numpy.cos(theta)))
dBHatdpsi=numpy.zeros((Npsi,Ntheta))
dBHatdtheta=psi_dependence*(inputfile.epsil*numpy.sin(theta)/((1 + inputfile.epsil * numpy.cos(theta))**2))
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
