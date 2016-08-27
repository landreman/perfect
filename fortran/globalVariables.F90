module globalVariables

  implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

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

  integer :: programMode = 1

  character(len=200) :: outputFilename = 'PERFECTOutput.h5'

  integer :: outputScheme = 1

  logical :: saveMatlabOutput = .false.

  character(len=200) :: MatlabOutputFilename = 'PERFECT.m'

  logical :: parallelizeOverScan = .true.

  logical :: solveSystem = .true.

  ! ********************************************************
  ! ********************************************************
  !
  ! Geometry input parameters:
  !
  ! ********************************************************
  ! ********************************************************

  integer :: geometryToUse = 0
  character(len=100) :: geometryFilename


  PetscScalar :: epsil = 0.1d+0

  ! Miller parameters: (these are only used when geometryToUse = 1.)
  PetscScalar :: Miller_kappa = 1.66d+0
  PetscScalar :: Miller_delta = 0.416d+0
  PetscScalar :: Miller_s_delta = 1.37d+0
  PetscScalar :: Miller_s_kappa = 0.7d+0
  PetscScalar :: Miller_dRdr = -0.354d+0
  PetscScalar :: Miller_q = 3.0d+0
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

  PetscScalar :: delta = 0.0011d+0
  PetscScalar :: omega = 0.0014d+0
  PetscScalar :: psiAHat
  PetscScalar, dimension(:), allocatable :: psiAHatArray

  PetscScalar :: psiMid, psiMin, psiMax

  PetscScalar :: nuPrime
  PetscScalar :: nuStar
  PetscScalar :: nu_r

  integer :: profilesScheme
  character(len=100) :: profilesFilename

  PetscScalar :: exponent

  logical :: setTPrimeToBalanceHeatFlux

  integer :: Nsources = 2
  
  integer :: sourcePoloidalVariation
  PetscScalar :: sourcePoloidalVariationStrength
  PetscScalar :: sourcePoloidalVariationPhase


  logical :: makeLocalApproximation

  logical :: includeCollisionOperator

  logical :: includeddpsiTerm

  PetscScalar :: leftBoundaryShift=0, rightBoundaryShift=0
  integer :: leftBoundaryScheme=1, rightBoundaryScheme=1

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
  PetscScalar :: thetaGridShift, scaledThetaGridShift

  integer :: Nxi
  PetscScalar :: NxiMinFactor, NxiMaxFactor
  integer :: NxiNumRuns
  integer, dimension(:), allocatable :: Nxi_for_x, min_x_for_L

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

  logical :: forceOddNtheta = .true.

  integer :: NxUniform = 25, NxiUniform = 31
  PetscScalar :: xUniformMax = 3d+0

  integer :: matrixSize, localMatrixSize, localDKEMatrixSize

  ! ********************************************************
  ! ********************************************************
  !
  ! Other numerical parameters:
  !
  ! ********************************************************
  ! ********************************************************

  ! for non-uniform grid
  integer :: psiGridType
  character(len=100) :: psiAHatFilename

  integer :: psiDerivativeScheme
  integer :: thetaDerivativeScheme
  integer :: xDerivativeScheme=2

  ! control treatment of constraints and sources at boundary

  integer :: NpsiSourcelessRight, NpsiSourcelessLeft

  ! lowest/highest psi indices where constraint are enforced
  integer :: lowestEnforcedIpsi, highestEnforcedIpsi

  PetscScalar :: thresh

  PetscScalar :: xScaleFactor = 1

  logical :: useIterativeSolver = .true.
  logical :: useIterativeBoundarySolver = .true.

  integer :: whichParallelSolverToFactorPreconditioner

  logical :: isAParallelDirectSolverInstalled

  integer :: layout = 0
  ! layout is not presently used.

  integer :: PETSCPreallocationStrategy = 1
  integer :: Nxi_for_x_option = 1

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
  PetscScalar, dimension(:,:), allocatable :: BHat, BPHat, BTHat, JHat, RHat, dBHatdtheta, dBHatdpsi
  PetscScalar, dimension(:), allocatable :: IHat, dIHatdpsi, dPhiHatdpsi, PhiHat
  PetscScalar, dimension(:,:), allocatable :: THats, dTHatdpsis, nHats, dnHatdpsis, etaHats, detaHatdpsis
  PetscScalar, dimension(:,:,:), allocatable :: sourceProfile
!  PetscScalar, dimension(:,:), allocatable :: LHSOfKParEquation
  PetscScalar, dimension(:), allocatable :: VPrimeHat, FSABHat2, typicalB
  PetscScalar, dimension(:,:,:), allocatable :: flow, kPar, densityPerturbation, pressurePerturbation
  PetscScalar, dimension(:,:,:), allocatable :: pPerpTermInVp,pPerpTermInVpBeforePsiDerivative
  PetscScalar, dimension(:,:,:), allocatable :: toroidalFlow,poloidalFlow
  PetscScalar, dimension(:,:,:), allocatable :: particleFluxBeforeThetaIntegral
  PetscScalar, dimension(:,:,:), allocatable :: momentumFluxBeforeThetaIntegral
  PetscScalar, dimension(:,:,:), allocatable :: heatFluxBeforeThetaIntegral
  PetscScalar, dimension(:,:), allocatable :: kParOutboard, kParInboard, FSAKPar
  PetscScalar, dimension(:,:), allocatable :: flowOutboard, flowInboard, FSAFlow, FSABFlow
  PetscScalar, dimension(:,:), allocatable :: FSAToroidalFlow,FSAPoloidalFlow
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
  PetscScalar, dimension(:), allocatable :: xUniform, xiUniform
  integer :: thetaIndexForOutboard, thetaIndexForInboard
  PetscScalar, dimension(:,:,:,:), allocatable :: deltaFOutboard, fullFOutboard

  PetscLogDouble :: elapsedTime
  integer :: didItConverge

  ! ********************************************************
  !
  !  Useful intermediate variables
  !
  ! ********************************************************

  PetscScalar, dimension(:,:), allocatable :: sqrtTHats
  PetscScalar, dimension(:,:,:,:,:,:), allocatable :: RosenbluthPotentialTerms

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
    deallocate(psiAHatArray)
    deallocate(xUniform)
    deallocate(xiUniform)
    deallocate(JHat)
    deallocate(RHat)
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

    deallocate(sqrtTHats)

    if (masterProcInSubComm) then
       deallocate(sourceProfile)
       deallocate(densityPerturbation)
       deallocate(flow)
       deallocate(toroidalFlow)
       deallocate(poloidalFlow)
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
       deallocate(FSAFlow)
       deallocate(FSABFlow)
       deallocate(FSAToroidalFlow)
       deallocate(FSAPoloidalFlow)
       deallocate(FSAPressurePerturbation)
!       deallocate(LHSOfKParEquation)
       deallocate(particleFlux)
       deallocate(momentumFlux)
       deallocate(heatFlux)
       deallocate(deltaFOutboard)
       deallocate(fullFOutboard)

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

