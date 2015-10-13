#!/usr/bin/env python

from perfectInputFile import perfectInput
from perfectProfilesFile import perfectProfiles
import numpy
import scipy.special
from sys import exit, argv

inputfilename = argv[1]
inputfile = perfectInput(inputfilename)
Npsi = inputfile.Npsi
numSpecies = inputfile.numSpecies
psi = inputfile.psi
psiMid = inputfile.psiMid
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
FSABHat2Fine = inputfile.FSABHat2[0]

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
    THats[species] = 1. + dTHatdpsiScalar * (psi-psiMid)
    dTHatdpsis[species] = dTHatdpsiScalar
    
    etaHats[species] = 1. + (psi-psiMid)*detaHatdpsiScalar
    detaHatdpsis[species] = ddpsi_accurate(etaHats[species])

    nHats[species] = etaHats[species] * numpy.exp(-2.*omega/delta*PhiHat/THats[species])
    dnHatdpsis[species] = ddpsi_accurate(nHats[species])

# Output profiles
outputfile = perfectProfiles(inputfile.profilesFilename)
outputfile.create_profiles_for_Npsi(Npsi,PhiHat,dPhiHatdpsi,THats,dTHatdpsis,nHats,dnHatdpsis,etaHats,detaHatdpsis)

exit(0)
