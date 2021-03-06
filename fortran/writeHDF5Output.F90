#define HAVE_PARALLEL_HDF5

module writeHDF5Output

use globalVariables
use scan
use petscsysdef
use HDF5

implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

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

interface rank
   module procedure rank_real
   module procedure rank_integer
   module procedure rank_character
   module procedure rank_scalar
   module procedure rank_1d
   module procedure rank_1d_integer
   !module procedure rank_1d_nonalloc
   module procedure rank_2d
   module procedure rank_3d
   module procedure rank_4d
   module procedure rank_5d
end interface rank


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

  subroutine setupOutput(gitCommit)

    implicit none

    integer :: i
    character(len=*), intent(in) :: gitCommit
    character(len=32) :: arg
    character(:),allocatable:: argString
    character(20) :: groupName

#ifdef HAVE_PARALLEL_HDF5
    if (outputScheme > 0) then
#else
    if (outputScheme > 0 .and. masterProcInSubComm) then
#endif
       allocate(groupIDs(numRunsInScan))

       ! Save programMode in the file:
       call writeIntegerNoGroup(programMode, "programMode")

       ! Save the git commit hash of the PERFECT binary in the file:
       call writeStringNoGroup(gitCommit,"gitCommit")

       ! Save the extra command line arguments in the file:
       argString=""
       do i = 1, iargc()
          call getarg(i, arg)
          if (i == 1) then
                argString=trim(arg)
          else
             argString=trim(argString)//" "//trim(arg)
          end if        
       end do
       if (len(argString)==0) then
          argString=" "
       end if
       call writeStringNoGroup(argString,"cmdFlags")

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
    ! to make the writeVariable interface recognize slices of our source array
    PetscScalar, dimension(:,:), allocatable :: genericSourceProfile
    allocate(genericSourceProfile(numSpecies,Npsi-NpsiSourcelessLeft-NpsiSourcelessRight))
    
    if ((psiGridType == 1) .and. masterProcInSubComm) then
       call convertToPsiNDerivatives()
       print *,"Converted output derivatives to psiN"
    end if

    if (outputScheme > 0) then      
      call writeVariable_1d_nonalloc(charges,numSpecies,"charges",runNum)
      call writeVariable_1d_nonalloc(masses,numSpecies,"masses",runNum)
      call writeVariable_1d_nonalloc(scalarTHats,numSpecies,"scalarTHats",runNum)
      call writeVariable_1d_nonalloc(scalarNHats,numSpecies,"scalarNHats",runNum)
      call writeVariable(NpsiPerDiameter,"NpsiPerDiameter",runNum)
      call writeVariable(Npsi,"Npsi",runNum)
      call writeVariable(psiDiameter,"psiDiameter",runNum)
      call writeVariable(widthExtender,"widthExtender",runNum)
      call writeVariable(numSpecies,"Nspecies",runNum)
      call writeVariable(Ntheta,"Ntheta",runNum)
      call writeVariable(thetaGridShift,"thetaGridShift",runNum)
      call writeVariable(Nxi,"Nxi",runNum)
      call writeIntegerVariable_1d(Nxi_for_x,"Nxi_for_x",runNum)
      call writeIntegerVariable_1d(min_x_for_L,"min_x_for_L",runNum)
      call writeVariable(NL,"NL",runNum)
      call writeVariable(Nx,"Nx",runNum)
      call writeVariable(NxPotentialsPerVth,"NxPotentialsPerVth",runNum)
      call writeVariable(xMax,"xMax",runNum)
      call writeVariable(xMaxForDistribution,"xMaxForDistribution",runNum)
      call writeVariable(solverTolerance,"solverTolerance",runNum)
      call writeVariable(Delta,"Delta",runNum)
      call writeVariable(omega,"omega",runNum)
      call writeVariable(psiAHat,"psiAHat",runNum)
      call writeVariable(psiGridType,"psiGridType",runNum)
      call writeVariable(psiAHatArray,"psiAHatArray",runNum)
      call writeVariable(nu_r,"nu_r",runNum)
      call writeVariable(Miller_q,"Miller_q",runNum)
      call writeVariable(epsil,"epsil",runNum)
      call writeVariable(psi,"psi",runNum)
      call writeVariable(theta,"theta",runNum)
      call writeVariable(JHat,"JHat",runNum)
      call writeVariable(BHat,"BHat",runNum)
      call writeVariable(dBHatdpsi,"d(BHat)d(psi)",runNum)
      call writeVariable(dBHatdtheta,"d(BHat)d(theta)",runNum)
      call writeVariable(BPHat,"BPHat",runNum)
      call writeVariable(BTHat,"BTHat",runNum)
      call writeVariable(RHat,"RHat",runNum)
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
      call writeVariable(NpsiSourcelessLeft,"NpsiSourcelessLeft",runNum)
      call writeVariable(NpsiSourcelessRight,"NpsiSourcelessRight",runNum)
      call writeVariable(sourcePoloidalVariation,"sourcePoloidalVariation",runNum)
      call writeVariable(sourcePoloidalVariationStrength,"sourcePoloidalVariationStrength",runNum)
      call writeVariable(sourcePoloidalVariationPhase,"sourcePoloidalVariationPhase",runNum)
      ! todo: implement an array of strings with source names
      ! in parts of code specifying velocity space structure of sources
      do temp = 1,Nsources
         select case(temp)
         case(1)
            ! particle source
            if (masterProcInSubComm) then
               genericSourceProfile = sourceProfile(temp,:,:)
            end if
            call writeVariable(genericSourceProfile,"particleSourceProfile",runNum)
         case(2)
            ! heat source
            if (masterProcInSubComm) then
               genericSourceProfile = sourceProfile(temp,:,:)
            end if
            call writeVariable(genericSourceProfile,"heatSourceProfile",runNum)
         case default
            ! as of now unspecified source
            if (masterProcInSubComm) then
               print *,"Writing source with isource > 2. This should not happen!"
               genericSourceProfile = sourceProfile(temp,:,:)
            end if
            call writeVariable(genericSourceProfile,"unknownSourceProfile",runNum)
         end select
      end do
      call writeVariable(VPrimeHat,"VPrimeHat",runNum)
      call writeVariable(FSABHat2,"FSABHat2",runNum)
      call writeVariable(U,"U",runNum)
      call writeVariable(r,"r",runNum)
      call writeVariable(densityPerturbation,"densityPerturbation",runNum)
      call writeVariable(flow,"flow",runNum)
      call writeVariable(kPar,"kPar",runNum)
      call writeVariable(pPerpTermInVp,"pPerpTermInVp",runNum)
      call writeVariable(pPerpTermInVpBeforePsiDerivative,"pPerpTermInVpBeforePsiDerivative",runNum)
      call writeVariable(toroidalFlow,"toroidalFlow",runNum)
      call writeVariable(poloidalFlow,"poloidalFlow",runNum)
      call writeVariable(pressurePerturbation,"pressurePerturbation",runNum)
      call writeVariable(particleFluxBeforeThetaIntegral,"particleFluxBeforeThetaIntegral",runNum)
      call writeVariable(momentumFluxBeforeThetaIntegral,"momentumFluxBeforeThetaIntegral",runNum)
      call writeVariable(heatFluxBeforeThetaIntegral,"heatFluxBeforeThetaIntegral",runNum)
      call writeVariable(FSADensityPerturbation,"FSADensityPerturbation",runNum)
      call writeVariable(flowOutboard,"flowOutboard",runNum)
      call writeVariable(flowInboard,"flowInboard",runNum)
      call writeVariable(FSAFlow,"FSAFlow",runNum)
      call writeVariable(FSAToroidalFlow,"FSAToroidalFlow",runNum)
      call writeVariable(FSAPoloidalFlow,"FSAPoloidalFlow",runNum)
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
      if (includeCollisionOperator) then
         temp = integerToRepresentTrue
      else
         temp = integerToRepresentFalse
      end if
      call writeVariable(temp,"includeCollisionOperator",runNum)
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
      call writeVariable(thetaIndexForInboard,"thetaIndexForInboard",runNum)
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
    call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
    if (masterProcInSubComm) then
      call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProcInSubComm) then
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, var, dimensions, HDF5Error)
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
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
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
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      dimensions = shape(var)
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeVariable_1d

  subroutine writeVariable_1d_nonalloc(var, varsize, varname, i)

    PetscScalar, dimension(:), intent(in) :: var
    integer, intent(in) :: varsize, i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(1) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (masterProcInSubComm) then
      dimensions = varsize
    end if
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPIComm,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(size(shape(var)), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (masterProcInSubComm) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProcInSubComm) then
      call h5screate_simple_f(size(shape(var)), dimensions, dspaceID, HDF5Error)
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
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
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
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
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
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
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
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
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
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
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
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
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
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
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
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeVariable_5d

  subroutine writeIntegerVariable_1d(var, varname, i)

    integer, dimension(:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(1) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (masterProcInSubComm) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      dimensions = shape(var)
    end if
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPIComm,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
    if (masterProcInSubComm) then
      call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProcInSubComm) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      dimensions = shape(var)
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(groupIDs(i), varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeIntegerVariable_1d

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

  subroutine writeStringNoGroup(var,varname)
    ! Only writes on masterProc

    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID, stringType
    integer(HSIZE_T), dimension(0) :: dimensions
    integer(HSIZE_T) :: stringLength
    integer :: ierror

    stringLength = len(var)
    dimensions = shape(var)

#ifdef HAVE_PARALLEL_HDF5
    call h5tcopy_f(H5T_FORTRAN_S1, stringType, HDF5Error)
    call h5tset_size_f(stringtype, stringLength, HDF5Error)
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, stringType, dspaceID, dsetID, HDF5Error)
    if (masterProc) then
      call h5dwrite_f(dsetID, stringType, var, dimensions, HDF5Error)
    end if
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (masterProc) then
      call h5tcopy_f(H5T_FORTRAN_S1, stringType, HDF5Error)
      call h5tset_size_f(stringtype, stringLength, HDF5Error)
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, stringType, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, stringType, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeStringNoGroup

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


  integer function rank_real(A)
    real, intent(in) :: A(:)
    rank_real=size(shape(A))
    return
  end function rank_real

  integer function rank_character(A)
    character(len=*), intent(in) :: A
    rank_character=size(shape(A))
    return
  end function rank_character

  integer function rank_integer(A)
    integer, intent(in) :: A
    rank_integer=size(shape(A))
    return
  end function rank_integer

  integer function rank_scalar(A)
    PetscScalar, intent(in) :: A
    rank_scalar=size(shape(A))
    return
  end function rank_scalar

  integer function rank_1d(A)
    PetscScalar, dimension(:), allocatable, intent(in) :: A
    rank_1d=size(shape(A))
    return
  end function rank_1d

  integer function rank_1d_integer(A)
    integer, dimension(:), allocatable, intent(in) :: A
    rank_1d_integer=size(shape(A))
    return
  end function rank_1d_integer

  
  ! conflicts with rank_1d
  !integer function rank_1d_nonalloc(A)
  !  PetscScalar, dimension(:), intent(in) :: A
  !  rank_1d_nonalloc=size(shape(A))
  !  return
  !end function rank_1d_nonalloc

  integer function rank_2d(A)
    PetscScalar, dimension(:,:), allocatable, intent(in) :: A
    rank_2d=size(shape(A))
    return
  end function rank_2d

  integer function rank_3d(A)
    PetscScalar, dimension(:,:,:), allocatable, intent(in) :: A
    rank_3d=size(shape(A))
    return
  end function rank_3d

  integer function rank_4d(A)
    PetscScalar, dimension(:,:,:,:), allocatable, intent(in) :: A
    rank_4d=size(shape(A))
    return
  end function rank_4d


  integer function rank_5d(A)
    PetscScalar, dimension(:,:,:,:,:), allocatable, intent(in) :: A
    rank_5d=size(shape(A))
    return
  end function rank_5d


 
  ! -------------------------------------------------------------------------------------
  
  subroutine convertToPsiNDerivatives()
    implicit none
    integer :: ipsi,itheta,ispecies
    PetscScalar :: scaleFactor
    
    do ipsi=1,Npsi
       scaleFactor = psiAHat/psiAHatArray(ipsi)
       !print *,"ipsi: ",ipsi
       ! psi
       dIHatdpsi(ipsi)=scaleFactor*dIHatdpsi(ipsi)
       dPhiHatdpsi(ipsi)=scaleFactor*dPhiHatdpsi(ipsi)

       !psi,theta
       do itheta=1,Ntheta
          dBHatdpsi(itheta,ipsi) = scaleFactor*dBHatdpsi(itheta,ipsi)
       end do
    
       !psi,species
       do ispecies=1,numSpecies
          dTHatdpsis(ispecies,ipsi)=scaleFactor*dTHatdpsis(ispecies,ipsi)
          dnHatdpsis(ispecies,ipsi)=scaleFactor*dnHatdpsis(ispecies,ipsi)
          detaHatdpsis(ispecies,ipsi)=scaleFactor*detaHatdpsis(ispecies,ipsi)
       end do
    end do
  end subroutine convertToPsiNDerivatives
  
end module writeHDF5Output

