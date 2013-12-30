#define HAVE_PARALLEL_HDF5

module writeHDF5Output

  use globalVariables
  use scan
  use petscsysdef
  use HDF5

  implicit none

#include <finclude/petscsysdef.h>

  integer, private :: HDF5Error
  integer(HID_T), private :: HDF5FileID, parallelID, dspaceIDForScalar
  integer(HID_T), dimension(:), allocatable, private :: dspaceIDForSpecies
  integer(HID_T), dimension(:), allocatable, private :: dspaceIDForProfile
  integer(HID_T), dimension(:), allocatable, private :: dspaceIDForSpeciesProfile
  integer(HID_T), dimension(:), allocatable, private :: dspaceIDForTheta
  integer(HID_T), dimension(:), allocatable, private :: dspaceIDForThetaPsi
  integer(HID_T), dimension(:), allocatable, private :: dspaceIDForSpeciesThetaPsi
  integer(HID_T), dimension(:), allocatable, private :: groupIDs

  integer(HID_T), private :: dsetID_programMode

  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_charges
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_masses
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_scalarTHats
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_scalarNHats
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_NpsiPerDiameter
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Npsi
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_psiDiameter
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_widthExtender
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Ntheta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Nxi
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_NL
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Nx
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_NxPotentialsPerVth
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_xMax
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_xMaxForDistribution
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_solverTolerance
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Delta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_omega
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_psiAHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_nu_r
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_Miller_q
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_epsil
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_psi
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_theta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_JHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_BHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dBHatdpsi
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dBHatdtheta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_IHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dIHatdpsi
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_PhiHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dPhiHatdpsi
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_THats
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dTHatdpsis
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_nHats
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dnHatdpsis
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_etaHats
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_detaHatdpsis
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_elapsedTime
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_sourcePoloidalVariation
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_particleSourceProfile
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_heatSourceProfile
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_VPrimeHat
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_FSABHat2
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_U
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_r
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_densityPerturbation
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_flow
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_kPar
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_pressurePerturbation
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_particleFluxBeforeThetaIntegral
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_momentumFluxBeforeThetaIntegral
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_heatFluxBeforeThetaIntegral
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_FSADensityPerturbation
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_flowOutboard
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_flowInboard
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_FSABFlow
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_kParOutboard
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_kParInboard
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_FSAKPar
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_FSAPressurePerturbation
!  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_LHSOfKParEquation
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_particleFlux
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_momentumFlux
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_heatFlux
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_makeLocalApproximation
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_useIterativeSolver
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_didItConverge
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_psiDerivativeScheme
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_thetaDerivativeScheme
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_xDerivativeScheme
!!$  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_kThetaWith3PointStencil
!!$  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_kThetaWith5PointStencil
!!$  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_kThetaOutboardWith3PointStencil
!!$  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_kThetaOutboardWith5PointStencil
!!$  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_kThetaInboardWith3PointStencil
!!$  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_kThetaInboardWith5PointStencil
!!$  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_PhiTermInKTheta
!!$  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_pPerpTermInKThetaBeforePsiDerivative
!!$  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_pPerpTermInKThetaWith3PointStencil
!!$  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_pPerpTermInKThetaWith5PointStencil
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_nuPrimeProfile
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_nuStarProfile
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_deltaN
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_deltaT
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_deltaEta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_desiredU
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_desiredFWHMInRhoTheta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_exponent
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_dTHatdpsiScalar
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_detaHatdpsiScalar
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_setTPrimeToBalanceHeatFlux
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_includeddpsiTerm
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_species
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_x
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_x_min_L
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_xi
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_theta
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_preconditioner_psi
  integer(HID_T), dimension(:), allocatable, private :: dsetIDs_layout


  !  integer(HSIZE_T), parameter, private :: dimForScalar = 1
  integer(HSIZE_T), dimension(1), parameter, private :: dimForScalar = 1
  integer(HSIZE_T), dimension(:,:), allocatable, private :: dimForSpecies
  integer(HSIZE_T), dimension(:,:), allocatable, private :: dimForProfile
  integer(HSIZE_T), dimension(:,:), allocatable, private :: dimForSpeciesProfile
  integer(HSIZE_T), dimension(:,:), allocatable, private :: dimForTheta
  integer(HSIZE_T), dimension(:,:), allocatable, private :: dimForThetaPsi
  integer(HSIZE_T), dimension(:,:), allocatable, private :: dimForSpeciesThetaPsi

contains

  ! -----------------------------------------------------------------------------------

  subroutine openOutputFile()

    implicit none

#ifdef HAVE_PARALLEL_HDF5
    if (outputScheme > 0) then
#else
    if (outputScheme > 0 .and. masterProcInSubComm) then
#endif
       call h5open_f(HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error initializing HDF5."
          stop
       end if

#ifdef HAVE_PARALLEL_HDF5
       ! Initialize some stuff related to parallel file access
       call h5pcreate_f(H5P_FILE_ACCESS_F, parallelID, HDF5Error)
       call h5pset_fapl_mpio_f(parallelID, MPI_COMM_WORLD, MPI_INFO_NULL, HDF5Error)

       call h5fcreate_f(outputFilename, H5F_ACC_TRUNC_F, HDF5FileID, HDF5Error, access_prp=parallelID)
#else
       call h5fcreate_f(outputFilename, H5F_ACC_TRUNC_F, HDF5FileID, HDF5Error)
#endif
       if (HDF5Error < 0) then
          print *,"Error opening HDF5 output file."
          stop
       end if

    end if

  end subroutine openOutputFile

  ! -----------------------------------------------------------------------------------

  subroutine createHDF5Structures()

    implicit none

    integer :: i, rank
    character(20) :: groupName

#ifdef HAVE_PARALLEL_HDF5
    if (outputScheme > 0) then
#else
    if (outputScheme > 0 .and. masterProcInSubComm) then
#endif
       allocate(groupIDs(numRunsInScan))

       allocate(dsetIDs_charges(numRunsInScan))
       allocate(dsetIDs_masses(numRunsInScan))
       allocate(dsetIDs_scalarTHats(numRunsInScan))
       allocate(dsetIDs_scalarNHats(numRunsInScan))
       allocate(dsetIDs_NpsiPerDiameter(numRunsInScan))
       allocate(dsetIDs_Npsi(numRunsInScan))
       allocate(dsetIDs_psiDiameter(numRunsInScan))
       allocate(dsetIDs_widthExtender(numRunsInScan))
       allocate(dsetIDs_Ntheta(numRunsInScan))
       allocate(dsetIDs_Nxi(numRunsInScan))
       allocate(dsetIDs_NL(numRunsInScan))
       allocate(dsetIDs_Nx(numRunsInScan))
       allocate(dsetIDs_NxPotentialsPerVth(numRunsInScan))
       allocate(dsetIDs_xMax(numRunsInScan))
       allocate(dsetIDs_xMaxForDistribution(numRunsInScan))
       allocate(dsetIDs_solverTolerance(numRunsInScan))
       allocate(dsetIDs_Delta(numRunsInScan))
       allocate(dsetIDs_omega(numRunsInScan))
       allocate(dsetIDs_psiAHat(numRunsInScan))
       allocate(dsetIDs_nu_r(numRunsInScan))
       allocate(dsetIDs_Miller_q(numRunsInScan))
       allocate(dsetIDs_epsil(numRunsInScan))
       allocate(dsetIDs_psi(numRunsInScan))
       allocate(dsetIDs_theta(numRunsInScan))
       allocate(dsetIDs_JHat(numRunsInScan))
       allocate(dsetIDs_BHat(numRunsInScan))
       allocate(dsetIDs_dBHatdpsi(numRunsInScan))
       allocate(dsetIDs_dBHatdtheta(numRunsInScan))
       allocate(dsetIDs_IHat(numRunsInScan))
       allocate(dsetIDs_dIHatdpsi(numRunsInScan))
       allocate(dsetIDs_PhiHat(numRunsInScan))
       allocate(dsetIDs_dPhiHatdpsi(numRunsInScan))
       allocate(dsetIDs_THats(numRunsInScan))
       allocate(dsetIDs_dTHatdpsis(numRunsInScan))
       allocate(dsetIDs_nHats(numRunsInScan))
       allocate(dsetIDs_dnHatdpsis(numRunsInScan))
       allocate(dsetIDs_etaHats(numRunsInScan))
       allocate(dsetIDs_detaHatdpsis(numRunsInScan))
       allocate(dsetIDs_elapsedTime(numRunsInScan))
       allocate(dsetIDs_sourcePoloidalVariation(numRunsInScan))
       allocate(dsetIDs_particleSourceProfile(numRunsInScan))
       allocate(dsetIDs_heatSourceProfile(numRunsInScan))
       allocate(dsetIDs_VPrimeHat(numRunsInScan))
       allocate(dsetIDs_FSABHat2(numRunsInScan))
       allocate(dsetIDs_U(numRunsInScan))
       allocate(dsetIDs_r(numRunsInScan))
       allocate(dsetIDs_densityPerturbation(numRunsInScan))
       allocate(dsetIDs_flow(numRunsInScan))
       allocate(dsetIDs_kPar(numRunsInScan))
       allocate(dsetIDs_pressurePerturbation(numRunsInScan))
       allocate(dsetIDs_particleFluxBeforeThetaIntegral(numRunsInScan))
       allocate(dsetIDs_momentumFluxBeforeThetaIntegral(numRunsInScan))
       allocate(dsetIDs_heatFluxBeforeThetaIntegral(numRunsInScan))
       allocate(dsetIDs_FSADensityPerturbation(numRunsInScan))
       allocate(dsetIDs_flowOutboard(numRunsInScan))
       allocate(dsetIDs_flowInboard(numRunsInScan))
       allocate(dsetIDs_FSABFlow(numRunsInScan))
       allocate(dsetIDs_kParOutboard(numRunsInScan))
       allocate(dsetIDs_kParInboard(numRunsInScan))
       allocate(dsetIDs_FSAKPar(numRunsInScan))
!       allocate(dsetIDs_LHSOfKParEquation(numRunsInScan))
       allocate(dsetIDs_FSAPressurePerturbation(numRunsInScan))
       allocate(dsetIDs_particleFlux(numRunsInScan))
       allocate(dsetIDs_momentumFlux(numRunsInScan))
       allocate(dsetIDs_heatFlux(numRunsInScan))
       allocate(dsetIDs_makeLocalApproximation(numRunsInScan))
       allocate(dsetIDs_useIterativeSolver(numRunsInScan))
       allocate(dsetIDs_didItConverge(numRunsInScan))
       allocate(dsetIDs_psiDerivativeScheme(numRunsInScan))
       allocate(dsetIDs_thetaDerivativeScheme(numRunsInScan))
       allocate(dsetIDs_xDerivativeScheme(numRunsInScan))
!!$       allocate(dsetIDs_kThetaWith3PointStencil(numRunsInScan))
!!$       allocate(dsetIDs_kThetaWith5PointStencil(numRunsInScan))
!!$       allocate(dsetIDs_kThetaOutboardWith3PointStencil(numRunsInScan))
!!$       allocate(dsetIDs_kThetaOutboardWith5PointStencil(numRunsInScan))
!!$       allocate(dsetIDs_kThetaInboardWith3PointStencil(numRunsInScan))
!!$       allocate(dsetIDs_kThetaInboardWith5PointStencil(numRunsInScan))
!!$       allocate(dsetIDs_PhiTermInKTheta(numRunsInScan))
!!$       allocate(dsetIDs_pPerpTermInKThetaBeforePsiDerivative(numRunsInScan))
!!$       allocate(dsetIDs_pPerpTermInKThetaWith3PointStencil(numRunsInScan))
!!$       allocate(dsetIDs_pPerpTermInKThetaWith5PointStencil(numRunsInScan))
       allocate(dsetIDs_nuPrimeProfile(numRunsInScan))
       allocate(dsetIDs_nuStarProfile(numRunsInScan))
       allocate(dsetIDs_deltaN(numRunsInScan))
       allocate(dsetIDs_deltaT(numRunsInScan))
       allocate(dsetIDs_deltaEta(numRunsInScan))
       allocate(dsetIDs_desiredU(numRunsInScan))
       allocate(dsetIDs_desiredFWHMInRhoTheta(numRunsInScan))
       allocate(dsetIDs_exponent(numRunsInScan))
       allocate(dsetIDs_dTHatdpsiScalar(numRunsInScan))
       allocate(dsetIDs_detaHatdpsiScalar(numRunsInScan))
       allocate(dsetIDs_setTPrimeToBalanceHeatFlux(numRunsInScan))
       allocate(dsetIDs_includeddpsiTerm(numRunsInScan))
       allocate(dsetIDs_preconditioner_species(numRunsInScan))
       allocate(dsetIDs_preconditioner_x(numRunsInScan))
       allocate(dsetIDs_preconditioner_x_min_L(numRunsInScan))
       allocate(dsetIDs_preconditioner_xi(numRunsInScan))
       allocate(dsetIDs_preconditioner_theta(numRunsInScan))
       allocate(dsetIDs_preconditioner_psi(numRunsInScan))
       allocate(dsetIDs_layout(numRunsInScan))

       allocate(dspaceIDForSpecies(numRunsInScan))
       allocate(dspaceIDForProfile(numRunsInScan))
       allocate(dspaceIDForSpeciesProfile(numRunsInScan))
       allocate(dspaceIDForTheta(numRunsInScan))
       allocate(dspaceIDForThetaPsi(numRunsInScan))
       allocate(dspaceIDForSpeciesThetaPsi(numRunsInScan))

       allocate(dimForSpecies(numRunsInScan,1))
       allocate(dimForProfile(numRunsInScan,1))
       allocate(dimForSpeciesProfile(numRunsInScan,2))
       allocate(dimForTheta(numRunsInScan,1))
       allocate(dimForThetaPsi(numRunsInScan,2))
       allocate(dimForSpeciesThetaPsi(numRunsInScan,3))

       ! Create a dataspace for storing single numbers:
       rank = 0
       call h5screate_simple_f(rank, dimForScalar, dspaceIDForScalar, HDF5Error)

       ! Save programMode in the file:
       call h5dcreate_f(HDF5FileID, "programMode", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
            dsetID_programMode, HDF5Error)
       if (masterProc) then
          call h5dwrite_f(dsetID_programMode, H5T_NATIVE_INTEGER, programMode, dimForScalar, HDF5Error)
       end if
       call h5dclose_f(dsetID_programMode, HDF5Error)

       do i=1,numRunsInScan
          ! Create a group to hold all data for the run:
          write (groupName, "(a, i3)") "run",i
          call h5gcreate_f(HDF5FileID, groupName, groupIDs(i), HDF5Error)

          ! Create dataspaces that depend on resolution parameters:
          rank = 1
          dimForProfile(i,1)=NpsisForScan(i)
          call h5screate_simple_f(rank, dimForProfile(i,:), dspaceIDForProfile(i), HDF5Error)

          dimForTheta(i,1)=NthetasForScan(i)
          call h5screate_simple_f(rank, dimForTheta(i,:), dspaceIDForTheta(i), HDF5Error)

          dimForSpecies(i,1)=numSpecies
          call h5screate_simple_f(rank, dimForSpecies(i,:), dspaceIDForSpecies(i), HDF5Error)

          rank = 2
          dimForSpeciesProfile(i,1)=numSpecies
          dimForSpeciesProfile(i,2)=NpsisForScan(i)
          call h5screate_simple_f(rank, dimForSpeciesProfile(i,:), dspaceIDForSpeciesProfile(i), HDF5Error)

          rank = 2
          dimForThetaPsi(i,1)=NthetasForScan(i)
          dimForThetaPsi(i,2)=NpsisForScan(i)
          call h5screate_simple_f(rank, dimForThetaPsi(i,:), dspaceIDForThetaPsi(i), HDF5Error)

          rank = 3
          dimForSpeciesThetaPsi(i,1)=numSpecies
          dimForSpeciesThetaPsi(i,2)=NthetasForScan(i)
          dimForSpeciesThetaPsi(i,3)=NpsisForScan(i)
          call h5screate_simple_f(rank, dimForSpeciesThetaPsi(i,:), dspaceIDForSpeciesThetaPsi(i), HDF5Error)

          ! Create datasets for each quantity in each run:
          call h5dcreate_f(groupIDs(i), "charges", H5T_NATIVE_DOUBLE, dspaceIDForSpecies(i), &
               dsetIDs_charges(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "masses", H5T_NATIVE_DOUBLE, dspaceIDForSpecies(i), &
               dsetIDs_masses(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "scalarTHats", H5T_NATIVE_DOUBLE, dspaceIDForSpecies(i), &
               dsetIDs_scalarTHats(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "scalarNHats", H5T_NATIVE_DOUBLE, dspaceIDForSpecies(i), &
               dsetIDs_scalarNHats(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "NpsiPerDiameter", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_NpsiPerDiameter(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Npsi", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_Npsi(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "psiDiameter", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_psiDiameter(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "widthExtender", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_widthExtender(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Ntheta", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_Ntheta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Nxi", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_Nxi(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "NL", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_NL(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Nx", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_Nx(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "NxPotentialsPerVth", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_NxPotentialsPerVth(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "xMax", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_xMax(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "xMaxForDistribution", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_xMaxForDistribution(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "solverTolerance", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_solverTolerance(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Delta", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_Delta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "omega", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_omega(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "psiAHat", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_psiAHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "nu_r", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_nu_r(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "Miller_q", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_Miller_q(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "epsil", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_epsil(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "psi", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
               dsetIDs_psi(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "theta", H5T_NATIVE_DOUBLE, dspaceIDForTheta(i), &
               dsetIDs_theta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "JHat", H5T_NATIVE_DOUBLE, dspaceIDForThetaPsi(i), &
               dsetIDs_JHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "BHat", H5T_NATIVE_DOUBLE, dspaceIDForThetaPsi(i), &
               dsetIDs_BHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(BHat)d(psi)", H5T_NATIVE_DOUBLE, dspaceIDForThetaPsi(i), &
               dsetIDs_dBHatdpsi(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(BHat)d(theta)", H5T_NATIVE_DOUBLE, dspaceIDForThetaPsi(i), &
               dsetIDs_dBHatdtheta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "IHat", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
               dsetIDs_IHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(IHat)d(psi)", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
               dsetIDs_dIHatdpsi(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "PhiHat", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
               dsetIDs_PhiHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(PhiHat)d(psi)", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
               dsetIDs_dPhiHatdpsi(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "THat", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_THats(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(THat)d(psi)", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_dTHatdpsis(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "nHat", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_nHats(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(nHat)d(psi)", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_dnHatdpsis(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "etaHat", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_etaHats(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "d(etaHat)d(psi)", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_detaHatdpsis(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "elapsed time (s)", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_elapsedTime(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "sourcePoloidalVariation", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_sourcePoloidalVariation(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "particleSourceProfile", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_particleSourceProfile(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "heatSourceProfile", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_heatSourceProfile(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "VPrimeHat", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
               dsetIDs_VPrimeHat(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "FSABHat2", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
               dsetIDs_FSABHat2(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "U", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_U(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "r", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_r(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "densityPerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaPsi(i), &
               dsetIDs_densityPerturbation(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "flow", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaPsi(i), &
               dsetIDs_flow(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "kPar", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaPsi(i), &
               dsetIDs_kPar(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "pressurePerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaPsi(i), &
               dsetIDs_pressurePerturbation(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "particleFluxBeforeThetaIntegral", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaPsi(i), &
               dsetIDs_particleFluxBeforeThetaIntegral(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "momentumFluxBeforeThetaIntegral", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaPsi(i), &
               dsetIDs_momentumFluxBeforeThetaIntegral(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "heatFluxBeforeThetaIntegral", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesThetaPsi(i), &
               dsetIDs_heatFluxBeforeThetaIntegral(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "FSADensityPerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_FSADensityPerturbation(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "flowOutboard", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_flowOutboard(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "flowInboard", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_flowInboard(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "FSABFlow", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_FSABFlow(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "kParOutboard", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_kParOutboard(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "kParInboard", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_kParInboard(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "FSAKPar", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_FSAKPar(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "FSAPressurePerturbation", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_FSAPressurePerturbation(i), HDF5Error)

!!$          call h5dcreate_f(groupIDs(i), "LHSOfKParEquation", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
!!$               dsetIDs_LHSOfKParEquation(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "particleFlux", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_particleFlux(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "momentumFlux", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_momentumFlux(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "heatFlux", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_heatFlux(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "makeLocalApproximation", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_makeLocalApproximation(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "useIterativeSolver", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_useIterativeSolver(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "didItConverge", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_didItConverge(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "psiDerivativeScheme", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_psiDerivativeScheme(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "thetaDerivativeScheme", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_thetaDerivativeScheme(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "xDerivativeScheme", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_xDerivativeScheme(i), HDF5Error)

!!$          call h5dcreate_f(groupIDs(i), "kThetaWith3PointStencil", H5T_NATIVE_DOUBLE, dspaceIDForThetaPsi(i), &
!!$               dsetIDs_kThetaWith3PointStencil(i), HDF5Error)
!!$
!!$          call h5dcreate_f(groupIDs(i), "kThetaWith5PointStencil", H5T_NATIVE_DOUBLE, dspaceIDForThetaPsi(i), &
!!$               dsetIDs_kThetaWith5PointStencil(i), HDF5Error)
!!$
!!$          call h5dcreate_f(groupIDs(i), "kThetaOutboardWith3PointStencil", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
!!$               dsetIDs_kThetaOutboardWith3PointStencil(i), HDF5Error)
!!$
!!$          call h5dcreate_f(groupIDs(i), "kThetaOutboardWith5PointStencil", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
!!$               dsetIDs_kThetaOutboardWith5PointStencil(i), HDF5Error)
!!$
!!$          call h5dcreate_f(groupIDs(i), "kThetaInboardWith3PointStencil", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
!!$               dsetIDs_kThetaInboardWith3PointStencil(i), HDF5Error)
!!$
!!$          call h5dcreate_f(groupIDs(i), "kThetaInboardWith5PointStencil", H5T_NATIVE_DOUBLE, dspaceIDForProfile(i), &
!!$               dsetIDs_kThetaInboardWith5PointStencil(i), HDF5Error)
!!$
!!$          call h5dcreate_f(groupIDs(i), "PhiTermInKTheta", H5T_NATIVE_DOUBLE, dspaceIDForThetaPsi(i), &
!!$               dsetIDs_PhiTermInKTheta(i), HDF5Error)
!!$
!!$          call h5dcreate_f(groupIDs(i), "pPerpTermInKThetaBeforePsiDerivative", H5T_NATIVE_DOUBLE, dspaceIDForThetaPsi(i), &
!!$               dsetIDs_pPerpTermInKThetaBeforePsiDerivative(i), HDF5Error)
!!$
!!$          call h5dcreate_f(groupIDs(i), "pPerpTermInKThetaWith3PointStencil", H5T_NATIVE_DOUBLE, dspaceIDForThetaPsi(i), &
!!$               dsetIDs_pPerpTermInKThetaWith3PointStencil(i), HDF5Error)
!!$
!!$          call h5dcreate_f(groupIDs(i), "pPerpTermInKThetaWith5PointStencil", H5T_NATIVE_DOUBLE, dspaceIDForThetaPsi(i), &
!!$               dsetIDs_pPerpTermInKThetaWith5PointStencil(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "nuPrimeProfile", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_nuPrimeProfile(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "nuStarProfile", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_nuStarProfile(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "deltaN", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_deltaN(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "deltaT", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_deltaT(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "deltaEta", H5T_NATIVE_DOUBLE, dspaceIDForSpeciesProfile(i), &
               dsetIDs_deltaEta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "desiredU", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_desiredU(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "desiredFWHMInRhoTheta", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_desiredFWHMInRhoTheta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "dTHatdpsiScalar", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_dTHatdpsiScalar(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "detaHatdpsiScalar", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_detaHatdpsiScalar(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "exponent", H5T_NATIVE_DOUBLE, dspaceIDForScalar, &
               dsetIDs_exponent(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "setTPrimeToBalanceHeatFlux", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_setTPrimeToBalanceHeatFlux(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "includeddpsiTerm", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_includeddpsiTerm(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_species", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_species(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_x", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_x(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_x_min_L", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_x_min_L(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_xi", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_xi(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_theta", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_theta(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "preconditioner_psi", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_preconditioner_psi(i), HDF5Error)

          call h5dcreate_f(groupIDs(i), "layout", H5T_NATIVE_INTEGER, dspaceIDForScalar, &
               dsetIDs_layout(i), HDF5Error)


       end do
    end if

  end subroutine createHDF5Structures

  ! -----------------------------------------------------------------------------------

  subroutine writeRunToOutputFile(runNum)

    implicit none

    integer, intent(in) :: runNum
    integer :: temp

    if (outputScheme > 0 .and. masterProcInSubComm) then

       call h5dwrite_f(dsetIDs_charges(runNum), H5T_NATIVE_DOUBLE, &
            charges, dimForSpecies(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_masses(runNum), H5T_NATIVE_DOUBLE, &
            masses, dimForSpecies(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_scalarTHats(runNum), H5T_NATIVE_DOUBLE, &
            scalarTHats, dimForSpecies(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_scalarNHats(runNum), H5T_NATIVE_DOUBLE, &
            scalarNHats, dimForSpecies(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_NpsiPerDiameter(runNum), H5T_NATIVE_DOUBLE, &
            NpsiPerDiameter, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Npsi(runNum), H5T_NATIVE_INTEGER, &
            Npsi, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_psiDiameter(runNum), H5T_NATIVE_DOUBLE, &
            psiDiameter, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_widthExtender(runNum), H5T_NATIVE_DOUBLE, &
            widthExtender, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Ntheta(runNum), H5T_NATIVE_INTEGER, &
            Ntheta, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Nxi(runNum), H5T_NATIVE_INTEGER, &
            Nxi, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_NL(runNum), H5T_NATIVE_INTEGER, &
            NL, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Nx(runNum), H5T_NATIVE_INTEGER, &
            Nx, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_NxPotentialsPerVth(runNum), H5T_NATIVE_DOUBLE, &
            NxPotentialsPerVth, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_xMax(runNum), H5T_NATIVE_DOUBLE, &
            xMax, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_xMaxForDistribution(runNum), H5T_NATIVE_DOUBLE, &
            xMaxForDistribution, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_solverTolerance(runNum), H5T_NATIVE_DOUBLE, &
            solverTolerance, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Delta(runNum), H5T_NATIVE_DOUBLE, &
            Delta, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_omega(runNum), H5T_NATIVE_DOUBLE, &
            omega, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_psiAHat(runNum), H5T_NATIVE_DOUBLE, &
            psiAHat, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_nu_r(runNum), H5T_NATIVE_DOUBLE, &
            nu_r, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_Miller_q(runNum), H5T_NATIVE_DOUBLE, &
            Miller_q, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_epsil(runNum), H5T_NATIVE_DOUBLE, &
            epsil, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_psi(runNum), H5T_NATIVE_DOUBLE, &
            psi, dimForProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_theta(runNum), H5T_NATIVE_DOUBLE, &
            theta, dimForTheta(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_JHat(runNum), H5T_NATIVE_DOUBLE, &
            JHat, dimForThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_BHat(runNum), H5T_NATIVE_DOUBLE, &
            BHat, dimForThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_dBHatdpsi(runNum), H5T_NATIVE_DOUBLE, &
            dBHatdpsi, dimForThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_dBHatdtheta(runNum), H5T_NATIVE_DOUBLE, &
            dBHatdtheta, dimForThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_IHat(runNum), H5T_NATIVE_DOUBLE, &
            IHat, dimForProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_dIHatdpsi(runNum), H5T_NATIVE_DOUBLE, &
            dIHatdpsi, dimForProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_PhiHat(runNum), H5T_NATIVE_DOUBLE, &
            phiHat, dimForProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_dPhiHatdpsi(runNum), H5T_NATIVE_DOUBLE, &
            dPhiHatdpsi, dimForProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_THats(runNum), H5T_NATIVE_DOUBLE, &
            THats, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_dTHatdpsis(runNum), H5T_NATIVE_DOUBLE, &
            dTHatdpsis, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_nHats(runNum), H5T_NATIVE_DOUBLE, &
            nHats, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_dnHatdpsis(runNum), H5T_NATIVE_DOUBLE, &
            dnHatdpsis, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_etaHats(runNum), H5T_NATIVE_DOUBLE, &
            etaHats, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_detaHatdpsis(runNum), H5T_NATIVE_DOUBLE, &
            detaHatdpsis, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_elapsedTime(runNum), H5T_NATIVE_DOUBLE, &
            elapsedTime, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_sourcePoloidalVariation(runNum), H5T_NATIVE_INTEGER, &
            sourcePoloidalVariation, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_particleSourceProfile(runNum), H5T_NATIVE_DOUBLE, &
            particleSourceProfile, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_heatSourceProfile(runNum), H5T_NATIVE_DOUBLE, &
            heatSourceProfile, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_VPrimeHat(runNum), H5T_NATIVE_DOUBLE, &
            VPrimeHat, dimForProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_FSABHat2(runNum), H5T_NATIVE_DOUBLE, &
            FSABHat2, dimForProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_U(runNum), H5T_NATIVE_DOUBLE, &
            U, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_r(runNum), H5T_NATIVE_DOUBLE, &
            r, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_densityPerturbation(runNum), H5T_NATIVE_DOUBLE, &
            densityPerturbation, dimForSpeciesThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_flow(runNum), H5T_NATIVE_DOUBLE, &
            flow, dimForSpeciesThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_kPar(runNum), H5T_NATIVE_DOUBLE, &
            kPar, dimForSpeciesThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_pressurePerturbation(runNum), H5T_NATIVE_DOUBLE, &
            pressurePerturbation, dimForSpeciesThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_particleFluxBeforeThetaIntegral(runNum), H5T_NATIVE_DOUBLE, &
            particleFluxBeforeThetaIntegral, dimForSpeciesThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_momentumFluxBeforeThetaIntegral(runNum), H5T_NATIVE_DOUBLE, &
            momentumFluxBeforeThetaIntegral, dimForSpeciesThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_heatFluxBeforeThetaIntegral(runNum), H5T_NATIVE_DOUBLE, &
            heatFluxBeforeThetaIntegral, dimForSpeciesThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_FSADensityPerturbation(runNum), H5T_NATIVE_DOUBLE, &
            FSADensityPerturbation, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_flowOutboard(runNum), H5T_NATIVE_DOUBLE, &
            flowOutboard, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_flowInboard(runNum), H5T_NATIVE_DOUBLE, &
            flowInboard, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_FSABFlow(runNum), H5T_NATIVE_DOUBLE, &
            FSABFlow, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_kParOutboard(runNum), H5T_NATIVE_DOUBLE, &
            kParOutboard, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_kParInboard(runNum), H5T_NATIVE_DOUBLE, &
            kParInboard, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_FSAKPar(runNum), H5T_NATIVE_DOUBLE, &
            FSAKPar, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_FSAPressurePerturbation(runNum), H5T_NATIVE_DOUBLE, &
            FSAPressurePerturbation, dimForSpeciesProfile(runNum,:), HDF5Error)

!!$       call h5dwrite_f(dsetIDs_LHSOfKParEquation(runNum), H5T_NATIVE_DOUBLE, &
!!$            LHSOfKParEquation, dimForProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_particleFlux(runNum), H5T_NATIVE_DOUBLE, &
            particleFlux, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_momentumFlux(runNum), H5T_NATIVE_DOUBLE, &
            momentumFlux, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_heatFlux(runNum), H5T_NATIVE_DOUBLE, &
            heatFlux, dimForSpeciesProfile(runNum,:), HDF5Error)

       if (makeLocalApproximation) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetIDs_makeLocalApproximation(runNum), H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       if (useIterativeSolver) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetIDs_useIterativeSolver(runNum), H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_didItConverge(runNum), H5T_NATIVE_INTEGER, &
            didItConverge, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_psiDerivativeScheme(runNum), H5T_NATIVE_INTEGER, &
            psiDerivativeScheme, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_thetaDerivativeScheme(runNum), H5T_NATIVE_INTEGER, &
            thetaDerivativeScheme, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_xDerivativeScheme(runNum), H5T_NATIVE_INTEGER, &
            xDerivativeScheme, dimForScalar, HDF5Error)

!!$       call h5dwrite_f(dsetIDs_kThetaWith3PointStencil(runNum), H5T_NATIVE_DOUBLE, &
!!$            kThetaWith3PointStencil, dimForThetaPsi(runNum,:), HDF5Error)
!!$
!!$       call h5dwrite_f(dsetIDs_kThetaWith5PointStencil(runNum), H5T_NATIVE_DOUBLE, &
!!$            kThetaWith5PointStencil, dimForThetaPsi(runNum,:), HDF5Error)
!!$
!!$       call h5dwrite_f(dsetIDs_kThetaOutboardWith3PointStencil(runNum), H5T_NATIVE_DOUBLE, &
!!$            kThetaOutboardWith3PointStencil, dimForProfile(runNum,:), HDF5Error)
!!$
!!$       call h5dwrite_f(dsetIDs_kThetaOutboardWith5PointStencil(runNum), H5T_NATIVE_DOUBLE, &
!!$            kThetaOutboardWith5PointStencil, dimForProfile(runNum,:), HDF5Error)
!!$
!!$       call h5dwrite_f(dsetIDs_kThetaInboardWith3PointStencil(runNum), H5T_NATIVE_DOUBLE, &
!!$            kThetaInboardWith3PointStencil, dimForProfile(runNum,:), HDF5Error)
!!$
!!$       call h5dwrite_f(dsetIDs_kThetaInboardWith5PointStencil(runNum), H5T_NATIVE_DOUBLE, &
!!$            kThetaInboardWith5PointStencil, dimForProfile(runNum,:), HDF5Error)
!!$
!!$       call h5dwrite_f(dsetIDs_PhiTermInKTheta(runNum), H5T_NATIVE_DOUBLE, &
!!$            PhiTermInKTheta, dimForThetaPsi(runNum,:), HDF5Error)
!!$
!!$       call h5dwrite_f(dsetIDs_pPerpTermInKThetaBeforePsiDerivative(runNum), H5T_NATIVE_DOUBLE, &
!!$            pPerpTermInKThetaBeforePsiDerivative, dimForThetaPsi(runNum,:), HDF5Error)
!!$
!!$       call h5dwrite_f(dsetIDs_pPerpTermInKThetaWith3PointStencil(runNum), H5T_NATIVE_DOUBLE, &
!!$            pPerpTermInKThetaWith3PointStencil, dimForThetaPsi(runNum,:), HDF5Error)
!!$
!!$       call h5dwrite_f(dsetIDs_pPerpTermInKThetaWith5PointStencil(runNum), H5T_NATIVE_DOUBLE, &
!!$            pPerpTermInKThetaWith5PointStencil, dimForThetaPsi(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_nuPrimeProfile(runNum), H5T_NATIVE_DOUBLE, &
            nuPrimeProfile, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_nuStarProfile(runNum), H5T_NATIVE_DOUBLE, &
            nuStarProfile, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_deltaN(runNum), H5T_NATIVE_DOUBLE, &
            deltaN, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_deltaT(runNum), H5T_NATIVE_DOUBLE, &
            deltaT, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_deltaEta(runNum), H5T_NATIVE_DOUBLE, &
            deltaEta, dimForSpeciesProfile(runNum,:), HDF5Error)

       call h5dwrite_f(dsetIDs_desiredU(runNum), H5T_NATIVE_DOUBLE, &
            desiredU, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_desiredFWHMInRhoTheta(runNum), H5T_NATIVE_DOUBLE, &
            desiredFWHMInRhoTheta, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_exponent(runNum), H5T_NATIVE_DOUBLE, &
            exponent, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_dTHatdpsiScalar(runNum), H5T_NATIVE_DOUBLE, &
            dTHatdpsiScalar, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_detaHatdpsiScalar(runNum), H5T_NATIVE_DOUBLE, &
            detaHatdpsiScalar, dimForScalar, HDF5Error)

       if (setTPrimeToBalanceHeatFlux) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetIDs_setTPrimeToBalanceHeatFlux(runNum), H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       if (includeddpsiTerm) then
          temp = integerToRepresentTrue
       else
          temp = integerToRepresentFalse
       end if
       call h5dwrite_f(dsetIDs_includeddpsiTerm(runNum), H5T_NATIVE_INTEGER, &
            temp, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_species(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_species, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_x(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_x, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_x_min_L(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_x_min_L, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_xi(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_xi, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_theta(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_theta, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_preconditioner_psi(runNum), H5T_NATIVE_INTEGER, &
            preconditioner_psi, dimForScalar, HDF5Error)

       call h5dwrite_f(dsetIDs_layout(runNum), H5T_NATIVE_INTEGER, &
            layout, dimForScalar, HDF5Error)

    end if

  end subroutine writeRunToOutputFile

  ! -----------------------------------------------------------------------------------

  subroutine closeOutputFile()

    implicit none

    integer :: i

    
#ifdef HAVE_PARALLEL_HDF5
    if (outputScheme > 0) then
#else
       if (outputScheme > 0 .and. masterProcInSubComm) then
#endif
       do i=1,numRunsInScan
          call h5dclose_f(dsetIDs_charges(i), HDF5Error)
          call h5dclose_f(dsetIDs_masses(i), HDF5Error)
          call h5dclose_f(dsetIDs_scalarTHats(i), HDF5Error)
          call h5dclose_f(dsetIDs_scalarNHats(i), HDF5Error)
          call h5dclose_f(dsetIDs_NpsiPerDiameter(i), HDF5Error)
          call h5dclose_f(dsetIDs_Npsi(i), HDF5Error)
          call h5dclose_f(dsetIDs_psiDiameter(i), HDF5Error)
          call h5dclose_f(dsetIDs_widthExtender(i), HDF5Error)
          call h5dclose_f(dsetIDs_Ntheta(i), HDF5Error)
          call h5dclose_f(dsetIDs_Nxi(i), HDF5Error)
          call h5dclose_f(dsetIDs_NL(i), HDF5Error)
          call h5dclose_f(dsetIDs_Nx(i), HDF5Error)
          call h5dclose_f(dsetIDs_NxPotentialsPerVth(i), HDF5Error)
          call h5dclose_f(dsetIDs_xMax(i), HDF5Error)
          call h5dclose_f(dsetIDs_xMaxForDistribution(i), HDF5Error)
          call h5dclose_f(dsetIDs_solverTolerance(i), HDF5Error)
          call h5dclose_f(dsetIDs_Delta(i), HDF5Error)
          call h5dclose_f(dsetIDs_omega(i), HDF5Error)
          call h5dclose_f(dsetIDs_psiAHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_nu_r(i), HDF5Error)
          call h5dclose_f(dsetIDs_Miller_q(i), HDF5Error)
          call h5dclose_f(dsetIDs_epsil(i), HDF5Error)
          call h5dclose_f(dsetIDs_psi(i), HDF5Error)
          call h5dclose_f(dsetIDs_theta(i), HDF5Error)
          call h5dclose_f(dsetIDs_JHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_BHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_dBHatdpsi(i), HDF5Error)
          call h5dclose_f(dsetIDs_dBHatdtheta(i), HDF5Error)
          call h5dclose_f(dsetIDs_IHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_dIHatdpsi(i), HDF5Error)
          call h5dclose_f(dsetIDs_PhiHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_dPhiHatdpsi(i), HDF5Error)
          call h5dclose_f(dsetIDs_THats(i), HDF5Error)
          call h5dclose_f(dsetIDs_dTHatdpsis(i), HDF5Error)
          call h5dclose_f(dsetIDs_nHats(i), HDF5Error)
          call h5dclose_f(dsetIDs_dnHatdpsis(i), HDF5Error)
          call h5dclose_f(dsetIDs_etaHats(i), HDF5Error)
          call h5dclose_f(dsetIDs_detaHatdpsis(i), HDF5Error)
          call h5dclose_f(dsetIDs_elapsedTime(i), HDF5Error)
          call h5dclose_f(dsetIDs_sourcePoloidalVariation(i), HDF5Error)
          call h5dclose_f(dsetIDs_particleSourceProfile(i), HDF5Error)
          call h5dclose_f(dsetIDs_heatSourceProfile(i), HDF5Error)
          call h5dclose_f(dsetIDs_VPrimeHat(i), HDF5Error)
          call h5dclose_f(dsetIDs_FSABHat2(i), HDF5Error)
          call h5dclose_f(dsetIDs_U(i), HDF5Error)
          call h5dclose_f(dsetIDs_r(i), HDF5Error)
          call h5dclose_f(dsetIDs_densityPerturbation(i), HDF5Error)
          call h5dclose_f(dsetIDs_flow(i), HDF5Error)
          call h5dclose_f(dsetIDs_kPar(i), HDF5Error)
          call h5dclose_f(dsetIDs_pressurePerturbation(i), HDF5Error)
          call h5dclose_f(dsetIDs_particleFluxBeforeThetaIntegral(i), HDF5Error)
          call h5dclose_f(dsetIDs_momentumFluxBeforeThetaIntegral(i), HDF5Error)
          call h5dclose_f(dsetIDs_heatFluxBeforeThetaIntegral(i), HDF5Error)
          call h5dclose_f(dsetIDs_FSADensityPerturbation(i), HDF5Error)
          call h5dclose_f(dsetIDs_flowOutboard(i), HDF5Error)
          call h5dclose_f(dsetIDs_flowInboard(i), HDF5Error)
          call h5dclose_f(dsetIDs_FSABFlow(i), HDF5Error)
          call h5dclose_f(dsetIDs_kParOutboard(i), HDF5Error)
          call h5dclose_f(dsetIDs_kParInboard(i), HDF5Error)
          call h5dclose_f(dsetIDs_FSAKPar(i), HDF5Error)
          call h5dclose_f(dsetIDs_FSAPressurePerturbation(i), HDF5Error)
!          call h5dclose_f(dsetIDs_LHSOfKParEquation(i), HDF5Error)
          call h5dclose_f(dsetIDs_particleFlux(i), HDF5Error)
          call h5dclose_f(dsetIDs_momentumFlux(i), HDF5Error)
          call h5dclose_f(dsetIDs_heatFlux(i), HDF5Error)
          call h5dclose_f(dsetIDs_makeLocalApproximation(i), HDF5Error)
          call h5dclose_f(dsetIDs_useIterativeSolver(i), HDF5Error)
          call h5dclose_f(dsetIDs_didItConverge(i), HDF5Error)
          call h5dclose_f(dsetIDs_psiDerivativeScheme(i), HDF5Error)
          call h5dclose_f(dsetIDs_thetaDerivativeScheme(i), HDF5Error)
          call h5dclose_f(dsetIDs_xDerivativeScheme(i), HDF5Error)
!!$          call h5dclose_f(dsetIDs_kThetaWith3PointStencil(i), HDF5Error)
!!$          call h5dclose_f(dsetIDs_kThetaWith5PointStencil(i), HDF5Error)
!!$          call h5dclose_f(dsetIDs_kThetaOutboardWith3PointStencil(i), HDF5Error)
!!$          call h5dclose_f(dsetIDs_kThetaOutboardWith5PointStencil(i), HDF5Error)
!!$          call h5dclose_f(dsetIDs_kThetaInboardWith3PointStencil(i), HDF5Error)
!!$          call h5dclose_f(dsetIDs_kThetaInboardWith5PointStencil(i), HDF5Error)
!!$          call h5dclose_f(dsetIDs_PhiTermInKTheta(i), HDF5Error)
!!$          call h5dclose_f(dsetIDs_pPerpTermInKThetaBeforePsiDerivative(i), HDF5Error)
!!$          call h5dclose_f(dsetIDs_pPerpTermInKThetaWith3PointStencil(i), HDF5Error)
!!$          call h5dclose_f(dsetIDs_pPerpTermInKThetaWith5PointStencil(i), HDF5Error)
          call h5dclose_f(dsetIDs_nuPrimeProfile(i), HDF5Error)
          call h5dclose_f(dsetIDs_nuStarProfile(i), HDF5Error)
          call h5dclose_f(dsetIDs_deltaN(i), HDF5Error)
          call h5dclose_f(dsetIDs_deltaT(i), HDF5Error)
          call h5dclose_f(dsetIDs_deltaEta(i), HDF5Error)
          call h5dclose_f(dsetIDs_desiredU(i), HDF5Error)
          call h5dclose_f(dsetIDs_desiredFWHMInRhoTheta(i), HDF5Error)
          call h5dclose_f(dsetIDs_setTPrimeToBalanceHeatFlux(i), HDF5Error)
          call h5dclose_f(dsetIDs_exponent(i), HDF5Error)
          call h5dclose_f(dsetIDs_dTHatdpsiScalar(i), HDF5Error)
          call h5dclose_f(dsetIDs_detaHatdpsiScalar(i), HDF5Error)
          call h5dclose_f(dsetIDs_includeddpsiTerm(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_species(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_x(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_x_min_L(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_xi(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_theta(i), HDF5Error)
          call h5dclose_f(dsetIDs_preconditioner_psi(i), HDF5Error)
          call h5dclose_f(dsetIDs_layout(i), HDF5Error)

          call h5gclose_f(groupIDs(i), HDF5Error)

!          call h5sclose_f(dspaceIDForSpecies(i), HDF5Error)
!          call h5sclose_f(dspaceIDForProfile(i), HDF5Error)
!          call h5sclose_f(dspaceIDForSpeciesProfile(i), HDF5Error)
!          call h5sclose_f(dspaceIDForTheta(i), HDF5Error)
!          call h5sclose_f(dspaceIDForThetaPsi(i), HDF5Error)
!          call h5sclose_f(dspaceIDForSpeciesThetaPsi(i), HDF5Error)
       end do

       call h5sclose_f(dspaceIDForScalar, HDF5Error)
       call h5pclose_f(parallelID, HDF5Error)
       call h5fclose_f(HDF5FileID, HDF5Error)
       call h5close_f(HDF5Error)
    end if

  end subroutine closeOutputFile

end module writeHDF5Output

