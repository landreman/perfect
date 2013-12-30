module globalVariables

  implicit none

#include <finclude/petscsysdef.h>

  integer, parameter :: integerToRepresentTrue  =  1
  integer, parameter :: integerToRepresentFalse = -1

  PetscScalar :: one = 1., oneHalf = 0.5d+0
  PetscScalar :: zero = 0., two = 2., three = 3., four = 4., five = 5.

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0

  ! ********************************************************
  ! ********************************************************
  !
  ! Options for program flow control:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: programMode

  character(len=200) :: outputFilename

  integer :: outputScheme

  logical :: saveMatlabOutput

  character(len=200) :: MatlabOutputFilename

  logical :: parallelizeOverScan

  logical :: solveSystem

  ! ********************************************************
  ! ********************************************************
  !
  ! Geometry input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: geometryToUse = 0

  PetscScalar :: epsil

  ! Miller parameters: (these are only used when geometryToUse = 1.)
  PetscScalar :: Miller_kappa
  PetscScalar :: Miller_delta
  PetscScalar :: Miller_s_delta
  PetscScalar :: Miller_s_kappa
  PetscScalar :: Miller_dRdr
  PetscScalar :: Miller_q
  ! The inverse aspect ratio epsil is also used for Miller geometry.

  ! ********************************************************
  ! ********************************************************
  !
  ! Physics input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  PetscScalar :: desiredU
  PetscScalar :: desiredUMin
  PetscScalar :: desiredUMax
  integer :: desiredUNumRuns
  PetscScalar :: desiredFWHMInRhoTheta
  PetscScalar :: dTHatdpsiScalar
  PetscScalar :: detaHatdpsiScalar

  PetscScalar :: delta
  PetscScalar :: omega
  PetscScalar :: psiAHat

  PetscScalar :: psiMid, psiMin, psiMax

  PetscScalar :: nuPrime
  PetscScalar :: nuStar
  PetscScalar :: nu_r

  integer :: profilesScheme

  PetscScalar :: exponent

  logical :: setTPrimeToBalanceHeatFlux

  integer :: sourcePoloidalVariation

  logical :: makeLocalApproximation

  logical :: includeddpsiTerm

  PetscScalar :: leftBoundaryShift, rightBoundaryShift
  integer :: leftBoundaryScheme, rightBoundaryScheme

  ! ********************************************************
  ! ********************************************************
  !
  ! Species quantities
  !
  ! ********************************************************
  ! ********************************************************

  integer, parameter :: maxNumSpecies = 100

  integer, parameter :: speciesNotInitialized = -9999

  integer :: numSpecies

  PetscScalar, dimension(maxNumSpecies) :: charges, masses, scalarNHats, scalarTHats

  ! ********************************************************
  ! ********************************************************
  !
  ! Numerical resolution parameters:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: Npsi
  PetscScalar :: NpsiPerDiameter
  PetscScalar :: NpsiMinFactor, NpsiMaxFactor
  integer :: NpsiNumRuns

  PetscScalar :: psiDiameter
  PetscScalar :: psiDiameterMinFactor, psiDiameterMaxFactor
  integer :: psiDiameterNumRuns

  PetscScalar :: widthExtender
  PetscScalar :: widthExtenderMin, widthExtenderMax
  integer :: widthExtenderNumRuns

  integer :: Ntheta
  PetscScalar :: NthetaMinFactor, NthetaMaxFactor
  integer :: NthetaNumRuns

  integer :: Nxi
  PetscScalar :: NxiMinFactor, NxiMaxFactor
  integer :: NxiNumRuns

  integer :: NL
  PetscScalar :: NLMinFactor, NLMaxFactor
  integer :: NLNumRuns

  integer :: Nx 
  PetscScalar :: NxMinFactor, NxMaxFactor
  integer :: NxNumRuns 

  PetscScalar  :: NxPotentialsPerVth
  PetscScalar :: NxPotentialsPerVthMinFactor, NxPotentialsPerVthMaxFactor
  integer :: NxPotentialsPerVthNumRuns

  PetscScalar :: xMax
  PetscScalar :: xMaxMinFactor, xMaxMaxFactor
  integer :: xMaxNumRuns

  PetscScalar :: xMaxForDistribution

  PetscScalar :: solverTolerance
  PetscScalar :: solverToleranceMinFactor
  PetscScalar :: solverToleranceMaxFactor
  integer :: solverToleranceNumRuns

  logical :: forceOddNtheta

  ! ********************************************************
  ! ********************************************************
  !
  ! Other numerical parameters:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: psiDerivativeScheme
  integer :: thetaDerivativeScheme
  integer :: xDerivativeScheme

  PetscScalar :: thresh

  PetscScalar :: xScaleFactor

  logical :: useIterativeSolver

  integer :: whichParallelSolverToFactorPreconditioner

  logical :: isAParallelDirectSolverInstalled

  integer :: layout
  ! layout is not presently used.

  ! ********************************************************
  ! ********************************************************
  ! 
  ! Preconditioner options:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: preconditioner_x, preconditioner_psi, preconditioner_species
  integer :: preconditioner_theta, preconditioner_xi, preconditioner_x_min_L

  ! ********************************************************
  ! ********************************************************
  !
  !  Outputs and numerical data which will be saved in the output file
  !
  ! ********************************************************
  ! ********************************************************

  PetscScalar, dimension(:), allocatable :: psi, theta
  PetscScalar, dimension(:,:), allocatable :: BHat, JHat, dBHatdtheta, dBHatdpsi
  PetscScalar, dimension(:), allocatable :: IHat, dIHatdpsi, dPhiHatdpsi, PhiHat
  PetscScalar, dimension(:,:), allocatable :: THats, dTHatdpsis, nHats, dnHatdpsis, etaHats, detaHatdpsis
  PetscScalar, dimension(:,:), allocatable :: particleSourceProfile, heatSourceProfile
!  PetscScalar, dimension(:,:), allocatable :: LHSOfKParEquation
  PetscScalar, dimension(:), allocatable :: VPrimeHat, FSABHat2, typicalB
  PetscScalar, dimension(:,:,:), allocatable :: flow, kPar, densityPerturbation, pressurePerturbation
  PetscScalar, dimension(:,:,:), allocatable :: particleFluxBeforeThetaIntegral
  PetscScalar, dimension(:,:,:), allocatable :: momentumFluxBeforeThetaIntegral
  PetscScalar, dimension(:,:,:), allocatable :: heatFluxBeforeThetaIntegral
  PetscScalar, dimension(:,:), allocatable :: kParOutboard, kParInboard, FSAKPar
  PetscScalar, dimension(:,:), allocatable :: flowOutboard, flowInboard, FSABFlow
  PetscScalar, dimension(:,:), allocatable :: FSADensityPerturbation, FSAPressurePerturbation
!  PetscScalar, dimension(:), allocatable :: kThetaOutboardWith3PointStencil, kThetaInboardWith3PointStencil
!  PetscScalar, dimension(:), allocatable :: kThetaOutboardWith5PointStencil, kThetaInboardWith5PointStencil
!  PetscScalar, dimension(:,:), allocatable :: kThetaWith3PointStencil, kThetaWith5PointStencil
  PetscScalar, dimension(:,:), allocatable :: particleFlux, heatFlux, momentumFlux
!  PetscScalar, dimension(:,:), allocatable :: potentialTermInPoloidalFlow
!  PetscScalar, dimension(:,:), allocatable :: PhiTermInKTheta
!  PetscScalar, dimension(:,:), allocatable :: pPerpTermInKThetaWith3PointStencil
!  PetscScalar, dimension(:,:), allocatable :: pPerpTermInKThetaWith5PointStencil
!  PetscScalar, dimension(:,:), allocatable :: pPerpTermInKThetaBeforePsiDerivative
  PetscScalar, dimension(:,:), allocatable :: nuPrimeProfile, nuStarProfile
  PetscScalar, dimension(:,:), allocatable :: deltaN, deltaT, deltaEta, U, r

  PetscLogDouble :: elapsedTime
  integer :: didItConverge

  ! ********************************************************
  !
  !  Variables related to parallelization:
  !
  ! ********************************************************

  integer :: numProcs, myRank 
  logical :: masterProc
  ! The above quantities refer to the PETSC_COMM_WORLD communicator, not to the smaller communicators
  ! used for parameter scans.

  ! The quantities below refer to the sub-communicator:
  integer :: numCommunicators
  integer, dimension(:), allocatable :: commMinProcs
  integer, dimension(:), allocatable :: minUnits, maxUnits
  MPI_Comm :: MPIComm
  integer :: myRankInSubComm, numProcsInSubComm
  logical :: masterProcInSubComm
  integer :: myCommunicatorIndex

contains

  ! ------------------------------------------------------------------------

  subroutine deallocateArrays()

    implicit none

    deallocate(psi)
    deallocate(theta)
    deallocate(JHat)
    deallocate(BHat)
    deallocate(dBHatdpsi)
    deallocate(dBHatdtheta)
    deallocate(IHat)
    deallocate(dIHatdpsi)
    deallocate(PhiHat)
    deallocate(dPhiHatdpsi)

    deallocate(THats)
    deallocate(dTHatdpsis)
    deallocate(nHats)
    deallocate(dnHatdpsis)
    deallocate(etaHats)
    deallocate(detaHatdpsis)
    deallocate(nuPrimeProfile)
    deallocate(nuStarProfile)
    deallocate(deltaN)
    deallocate(deltaT)
    deallocate(deltaEta)
    deallocate(VPrimeHat)
    deallocate(FSABHat2)
    deallocate(typicalB)
    deallocate(U)
    deallocate(r)

    if (masterProcInSubComm) then
       deallocate(particleSourceProfile)
       deallocate(heatSourceProfile)

       deallocate(densityPerturbation)
       deallocate(flow)
       deallocate(kPar)
       deallocate(pressurePerturbation)
       deallocate(particleFluxBeforeThetaIntegral)
       deallocate(momentumFluxBeforeThetaIntegral)
       deallocate(heatFluxBeforeThetaIntegral)

       deallocate(FSADensityPerturbation)
       deallocate(kParOutboard)
       deallocate(kParInboard)
       deallocate(FSAKPar)
       deallocate(flowOutboard)
       deallocate(flowInboard)
       deallocate(FSABFlow)
       deallocate(FSAPressurePerturbation)
!       deallocate(LHSOfKParEquation)
       deallocate(particleFlux)
       deallocate(momentumFlux)
       deallocate(heatFlux)

!       deallocate(kThetaWith3PointStencil)
!       deallocate(kThetaWith5PointStencil)
!       deallocate(kThetaOutboardWith3PointStencil)
!       deallocate(kThetaInboardWith3PointStencil)
!       deallocate(kThetaOutboardWith5PointStencil)
!       deallocate(kThetaInboardWith5PointStencil)
!       deallocate(PhiTermInKTheta)
!       deallocate(pPerpTermInKThetaWith3PointStencil)
!       deallocate(pPerpTermInKThetaWith5PointStencil)
!       deallocate(pPerpTermInKThetaBeforePsiDerivative)
    end if

  end subroutine deallocateArrays

end module globalVariables

