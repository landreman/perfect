#define HAVE_PARALLEL_HDF5

module writeHDF5Output

use globalVariables
use scan
use petscsysdef
use HDF5

implicit none

#include <finclude/petscsysdef.h>

integer, private :: HDF5Error
integer(HID_T), private :: HDF5FileID, parallelID
integer(HID_T), dimension(:), allocatable, private :: groupIDs
interface writeVariable
  module procedure writeVariable_integer
  module procedure writeVariable_scalar
  module procedure writeVariable_1d
  module procedure writeVariable_2d
  module procedure writeVariable_3d
  module procedure writeVariable_4d
  module procedure writeVariable_5d
end interface writeVariable

! Stuff for writing debugging output:
integer(HID_T), private :: HDF5DebugFileID
interface writeDebugArray
  module procedure writeDebugArray_scalar
  module procedure writeDebugArray_1d
  module procedure writeDebugArray_2d
  module procedure writeDebugArray_3d
  module procedure writeDebugArray_4d
  module procedure writeDebugArray_5d
end interface writeDebugArray

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

  subroutine setupOutput()

    implicit none

    integer :: i, rank
    character(20) :: groupName

#ifdef HAVE_PARALLEL_HDF5
    if (outputScheme > 0) then
#else
    if (outputScheme > 0 .and. masterProcInSubComm) then
#endif
       allocate(groupIDs(numRunsInScan))

       ! Save programMode in the file:
       call writeIntegerNoGroup(programMode, "programMode")

       do i=1,numRunsInScan
          ! Create a group to hold all data for the run:
          write (groupName, "(a, i3)") "run",i
          call h5gcreate_f(HDF5FileID, groupName, groupIDs(i), HDF5Error)
       end do
    end if

  end subroutine setupOutput

  ! -----------------------------------------------------------------------------------

  subroutine writeRunToOutputFile(runNum)

    implicit none

    integer, intent(in) :: runNum
    integer :: temp

    if (outputScheme > 0) then

      call writeVariable_1d_nonalloc(charges,"charges",runNum)
      call writeVariable_1d_nonalloc(masses,"masses",runNum)
      call writeVariable_1d_nonalloc(scalarTHats,"scalarTHats",runNum)
      call writeVariable_1d_nonalloc(scalarNHats,"scalarNHats",runNum)
      call writeVariable(NpsiPerDiameter,"NpsiPerDiameter",runNum)
      call writeVariable(Npsi,"Npsi",runNum)
      call writeVariable(psiDiameter,"psiDiameter",runNum)
      call writeVariable(widthExtender,"widthExtender",runNum)
      call writeVariable(numSpecies,"Nspecies",runNum)
      call writeVariable(Ntheta,"Ntheta",runNum)
      call writeVariable(Nxi,"Nxi",runNum)
      call writeVariable(NL,"NL",runNum)
      call writeVariable(Nx,"Nx",runNum)
      call writeVariable(NxPotentialsPerVth,"NxPotentialsPerVth",runNum)
      call writeVariable(xMax,"xMax",runNum)
      call writeVariable(xMaxForDistribution,"xMaxForDistribution",runNum)
      call writeVariable(solverTolerance,"solverTolerance",runNum)
      call writeVariable(Delta,"Delta",runNum)
      call writeVariable(omega,"omega",runNum)
      call writeVariable(psiAHat,"psiAHat",runNum)
      call writeVariable(nu_r,"nu_r",runNum)
      call writeVariable(Miller_q,"Miller_q",runNum)
      call writeVariable(epsil,"epsil",runNum)
      call writeVariable(psi,"psi",runNum)
      call writeVariable(theta,"theta",runNum)
      call writeVariable(JHat,"JHat",runNum)
      call writeVariable(BHat,"BHat",runNum)
      call writeVariable(dBHatdpsi,"d(BHat)d(psi)",runNum)
      call writeVariable(dBHatdtheta,"d(BHat)d(theta)",runNum)
      call writeVariable(IHat,"IHat",runNum)
      call writeVariable(dIHatdpsi,"d(IHat)d(psi)",runNum)
      call writeVariable(PhiHat,"PhiHat",runNum)
      call writeVariable(dPhiHatdpsi,"d(PhiHat)d(psi)",runNum)
      call writeVariable(THats,"THat",runNum)
      call writeVariable(dTHatdpsis,"d(THat)d(psi)",runNum)
      call writeVariable(nHats,"nHat",runNum)
      call writeVariable(dnHatdpsis,"d(nHat)d(psi)",runNum)
      call writeVariable(etaHats,"etaHat",runNum)
      call writeVariable(detaHatdpsis,"d(etaHat)d(psi)",runNum)
      call writeVariable(elapsedTime,"elapsed time (s)",runNum)
      call writeVariable(sourcePoloidalVariation,"sourcePoloidalVariation",runNum)
      call writeVariable(particleSourceProfile,"particleSourceProfile",runNum)
      call writeVariable(heatSourceProfile,"heatSourceProfile",runNum)
      call writeVariable(VPrimeHat,"VPrimeHat",runNum)
      call writeVariable(FSABHat2,"FSABHat2",runNum)
      call writeVariable(U,"U",runNum)
      call writeVariable(r,"r",runNum)
      call writeVariable(densityPerturbation,"densityPerturbation",runNum)
      call writeVariable(flow,"flow",runNum)
      call writeVariable(kPar,"kPar",runNum)
      call writeVariable(pressurePerturbation,"pressurePerturbation",runNum)
      call writeVariable(particleFluxBeforeThetaIntegral,"particleFluxBeforeThetaIntegral",runNum)
      call writeVariable(momentumFluxBeforeThetaIntegral,"momentumFluxBeforeThetaIntegral",runNum)
      call writeVariable(heatFluxBeforeThetaIntegral,"heatFluxBeforeThetaIntegral",runNum)
      call writeVariable(FSADensityPerturbation,"FSADensityPerturbation",runNum)
      call writeVariable(flowOutboard,"flowOutboard",runNum)
      call writeVariable(flowInboard,"flowInboard",runNum)
      call writeVariable(FSABFlow,"FSABFlow",runNum)
      call writeVariable(kParOutboard,"kParOutboard",runNum)
      call writeVariable(kParInboard,"kParInboard",runNum)
      call writeVariable(FSAKPar,"FSAKPar",runNum)
      call writeVariable(FSAPressurePerturbation,"FSAPressurePerturbation",runNum)
  !!$    call writeVariable(LHSOfKParEquation,"LHSOfKParEquation",runNum)
      call writeVariable(particleFlux,"particleFlux",runNum)
      call writeVariable(momentumFlux,"momentumFlux",runNum)
      call writeVariable(heatFlux,"heatFlux",runNum)
      if (makeLocalApproximation) then
         temp = integerToRepresentTrue
      else
         temp = integerToRepresentFalse
      end if
      call writeVariable(temp,"makeLocalApproximation",runNum)
      if (useIterativeSolver) then
         temp = integerToRepresentTrue
      else
         temp = integerToRepresentFalse
      end if
      call writeVariable(temp,"useIterativeSolver",runNum)
      call writeVariable(didItConverge,"didItConverge",runNum)
      call writeVariable(psiDerivativeScheme,"psiDerivativeScheme",runNum)
      call writeVariable(thetaDerivativeScheme,"thetaDerivativeScheme",runNum)
      call writeVariable(xDerivativeScheme,"xDerivativeScheme",runNum)
  !!$    call writeVariable(kThetaWith3PointStencil,"kThetaWith3PointStencil",runNum)
  !!$    call writeVariable(kThetaWith5PointStencil,"kThetaWith5PointStencil",runNum)
  !!$    call writeVariable(kThetaOutboardWith3PointStencil,"kThetaOutboardWith3PointStencil",runNum)
  !!$    call writeVariable(kThetaOutboardWith5PointStencil,"kThetaOutboardWith5PointStencil",runNum)
  !!$    call writeVariable(kThetaInboardWith3PointStencil,"kThetaInboardWith3PointStencil",runNum)
  !!$    call writeVariable(kThetaInboardWith5PointStencil,"kThetaInboardWith5PointStencil",runNum)
  !!$    call writeVariable(PhiTermInKTheta,"PhiTermInKTheta",runNum)
  !!$    call writeVariable(pPerpTermInKThetaBeforePsiDerivative,"pPerpTermInKThetaBeforePsiDerivative",runNum)
  !!$    call writeVariable(pPerpTermInKThetaWith3PointStencil,"pPerpTermInKThetaWith3PointStencil",runNum)
  !!$    call writeVariable(pPerpTermInKThetaWith5PointStencil,"pPerpTermInKThetaWith5PointStencil",runNum)
      call writeVariable(nuPrimeProfile,"nuPrimeProfile",runNum)
      call writeVariable(nuStarProfile,"nuStarProfile",runNum)
      call writeVariable(deltaN,"deltaN",runNum)
      call writeVariable(deltaT,"deltaT",runNum)
      call writeVariable(deltaEta,"deltaEta",runNum)
      call writeVariable(desiredU,"desiredU",runNum)
      call writeVariable(desiredFWHMInRhoTheta,"desiredFWHMInRhoTheta",runNum)
      call writeVariable(dTHatdpsiScalar,"dTHatdpsiScalar",runNum)
      call writeVariable(detaHatdpsiScalar,"detaHatdpsiScalar",runNum)
      call writeVariable(exponent,"exponent",runNum)
      if (setTPrimeToBalanceHeatFlux) then
         temp = integerToRepresentTrue
      else
         temp = integerToRepresentFalse
      end if
      call writeVariable(temp,"setTPrimeToBalanceHeatFlux",runNum)
      if (includeddpsiTerm) then
         temp = integerToRepresentTrue
      else
         temp = integerToRepresentFalse
      end if
      call writeVariable(temp,"includeddpsiTerm",runNum)
      call writeVariable(preconditioner_species,"preconditioner_species",runNum)
      call writeVariable(preconditioner_x,"preconditioner_x",runNum)
      call writeVariable(preconditioner_x_min_L,"preconditioner_x_min_L",runNum)
      call writeVariable(preconditioner_xi,"preconditioner_xi",runNum)
      call writeVariable(preconditioner_theta,"preconditioner_theta",runNum)
      call writeVariable(preconditioner_psi,"preconditioner_psi",runNum)
      call writeVariable(layout,"layout",runNum)
      call writeVariable(NxUniform,"NxUniform",runNum)
      call writeVariable(NxiUniform,"NxiUniform",runNum)
      call writeVariable(xUniform,"xUniform",runNum)
      call writeVariable(xiUniform,"xiUniform",runNum)
      call writeVariable(thetaIndexForOutboard,"thetaIndexForOutboard",runNum)
      if (outputScheme > 1) then
        call writeVariable(deltaFOutboard,"deltaFOutboard",runNum)
        call writeVariable(fullFOutboard,"fullFOutboard",runNum)
      end if
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
          call h5gclose_f(groupIDs(i), HDF5Error)
       end do

       call h5pclose_f(parallelID, HDF5Error)
       call h5fclose_f(HDF5FileID, HDF5Error)
       call h5close_f(HDF5Error)
    end if

  end subroutine closeOutputFile

  subroutine writeVariable_scalar(var, varname, i)

    PetscScalar, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(0) :: dimensions

    dimensions = shape(var)

#ifdef HAVE_PARALLEL_HDF5
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (masterProcInSubComm) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProcInSubComm) then
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeVariable_scalar

  subroutine writeVariable_integer(var, varname, i)

    integer, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(0) :: dimensions

    dimensions = shape(var)

#ifdef HAVE_PARALLEL_HDF5
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (masterProcInSubComm) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProcInSubComm) then
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_INTEGER dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeVariable_integer

  subroutine writeVariable_1d(var, varname, i)

    PetscScalar, dimension(:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(1) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (masterProcInSubComm) then
      dimensions = shape(var)
    end if
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPIComm,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (masterProcInSubComm) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProcInSubComm) then
      dimensions = shape(var)
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeVariable_1d

  subroutine writeVariable_1d_nonalloc(var, varname, i)

    PetscScalar, dimension(:), intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(1) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (masterProcInSubComm) then
      dimensions = shape(var)
    end if
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPIComm,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (masterProcInSubComm) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProcInSubComm) then
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeVariable_1d_nonalloc

  subroutine writeVariable_2d(var, varname, i)

    PetscScalar, dimension(:,:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(2) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (masterProcInSubComm) then
      dimensions = shape(var)
    end if
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPIComm,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (masterProcInSubComm) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProcInSubComm) then
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeVariable_2d

  subroutine writeVariable_3d(var, varname, i)

    PetscScalar, dimension(:,:,:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(3) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (masterProcInSubComm) then
      dimensions = shape(var)
    end if
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPIComm,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (masterProcInSubComm) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProcInSubComm) then
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeVariable_3d

  subroutine writeVariable_4d(var, varname, i)

    PetscScalar, dimension(:,:,:,:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(4) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (masterProcInSubComm) then
      dimensions = shape(var)
    end if
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPIComm,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (masterProcInSubComm) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProcInSubComm) then
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeVariable_4d

  subroutine writeVariable_5d(var, varname, i)

    PetscScalar, dimension(:,:,:,:,:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(5) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (masterProcInSubComm) then
      dimensions = shape(var)
    end if
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPIComm,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (masterProcInSubComm) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProcInSubComm) then
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeVariable_5d

  subroutine writeIntegerNoGroup(var,varname)
    ! Only writes on masterProc

    Integer, intent(in) :: var
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(0) :: dimensions
    integer :: ierror

    dimensions = shape(var)

#ifdef HAVE_PARALLEL_HDF5
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
    if (masterProc) then
      call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProc) then
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeIntegerNoGroup

  subroutine openDebugOutputFile(filename)

    character(len=*), intent(in) :: filename

  !#ifdef HAVE_PARALLEL_HDF5
  !    if (outputScheme > 0) then
  !#else
  !    if (outputScheme > 0 .and. masterProcInSubComm) then
  !#endif
      !print *, "opening",filename
      call h5open_f(HDF5Error)
      if (HDF5Error < 0) then
         print *,"Error initializing HDF5."
         stop
      end if

  !#ifdef HAVE_PARALLEL_HDF5
  !       ! Initialize some stuff related to parallel file access
  !       call h5pcreate_f(H5P_FILE_ACCESS_F, parallelID, HDF5Error)
  !       call h5pset_fapl_mpio_f(parallelID, MPI_COMM_WORLD, MPI_INFO_NULL, HDF5Error)
  !
  !       call h5fcreate_f(filename, H5F_ACC_TRUNC_F, HDF5DebugFileID, HDF5Error, access_prp=parallelID)
  !#else
       call h5fcreate_f(filename, H5F_ACC_TRUNC_F, HDF5DebugFileID, HDF5Error)
  !#endif
       if (HDF5Error < 0) then
          print *,"Error opening HDF5 output file."
          stop
       end if

  !    end if

  end subroutine openDebugOutputFile

  subroutine closeDebugOutputFile()

    !print *,"closing"
    !call h5pclose_f(parallelID,HDF5Error)
    call h5fclose_f(HDF5DebugFileID,HDF5Error)
    !! call h5close_f(HDF5Error) ! breaks writing more output later???

  end subroutine closeDebugOutputFile

  subroutine writeDebugArray_scalar(var,varname)

    PetscScalar, intent(in) :: var
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(:), allocatable :: dimensions

    dimensions = shape(var)

    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5DebugFileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)

  end subroutine writeDebugArray_scalar

  subroutine writeDebugArray_1d(var,varname)

    PetscScalar, dimension(:), allocatable, intent(in) :: var
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(:), allocatable :: dimensions

    dimensions = shape(var)

    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5DebugFileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)

  end subroutine writeDebugArray_1d

  subroutine writeDebugArray_2d(var,varname)

    PetscScalar, dimension(:,:), allocatable, intent(in) :: var
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(:), allocatable :: dimensions

    dimensions = shape(var)

    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5DebugFileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)

  end subroutine writeDebugArray_2d

  subroutine writeDebugArray_3d(var,varname)

    PetscScalar, dimension(:,:,:), allocatable, intent(in) :: var
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(:), allocatable :: dimensions

    dimensions = shape(var)

    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5DebugFileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)

  end subroutine writeDebugArray_3d

  subroutine writeDebugArray_4d(var,varname)

    PetscScalar, dimension(:,:,:,:), allocatable, intent(in) :: var
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(:), allocatable :: dimensions

    dimensions = shape(var)

    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5DebugFileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)

  end subroutine writeDebugArray_4d

  subroutine writeDebugArray_5d(var,varname)

    PetscScalar, dimension(:,:,:,:,:), allocatable, intent(in) :: var
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(:), allocatable :: dimensions

    dimensions = shape(var)

    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5DebugFileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)

  end subroutine writeDebugArray_5d

end module writeHDF5Output

