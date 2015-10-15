module geometry
  ! Subroutines and functions related to determining B(theta).

  use globalVariables
  use grids
  use HDF5
#include <finclude/petscsysdef.h>


  implicit none

  logical, private :: initializedYet = .false.
  PetscScalar :: Miller_A, Miller_x, Miller_QQ
  integer, parameter :: NThetaIntegral = 100
  integer :: HDF5Error
  integer(HID_T) :: HDF5FileID, HDF5GroupID, HDF5DatasetID, parallelID
  character(len=100) :: HDF5Groupname
  integer(HSIZE_T), dimension(1) :: arraydims1d
  integer(HSIZE_T), dimension(2) :: arraydims2d


contains

  !-----------------------------------------------------------------------------------------

  subroutine initializeGeometry()
    ! In case there are any calculations related to the magnetic geometry which
    ! only need to be done once, even for a parameter scan, they can go in this subroutine.
    ! This subroutine is called only once at the beginning of perfect.F90.

    implicit none

    integer :: i

    select case (geometryToUse)
    case (0)
       ! Circular concentric flux surfaces

    case (1)
       ! Miller geometry
       Miller_x = asin(Miller_delta)
       Miller_A = 1/epsil
       Miller_QQ = 0
       ! Integrate QQIntegrand from 0 to 2*pi.
       do i=1,NThetaIntegral
          Miller_QQ = Miller_QQ + QQIntegrand(2*pi*i/NThetaIntegral)
       end do
       Miller_QQ = Miller_kappa / (2*pi*Miller_A) * (Miller_QQ * 2*pi/NThetaIntegral)

    case (2)
       ! Circular concentric flux surfaces with Boozer poloidal angle

    case (3)
       ! EFIT interface
       print *,"Error! EFIT interface is not yet implemented"
       stop

       ! Here, put any steps that only need to be done once even in a parameter scan, such
       ! as reading in data from a file.

    case (4)
       !Read hdf5 file
       if (.not. len(geometryFilename)>=0) then
          print *,"If geometryToUse==4 then geometryFilename must be set."
          stop
       end if

    case default
       print *,"Error! Invalid geometry."
       stop
    end select

    initializedYet = .true.
  end subroutine initializeGeometry

  !-----------------------------------------------------------------------------------------
  subroutine computeMagneticQuantitiesOnGrids()
    ! This is the central subroutine of this module.  It is called at the beginning of
    ! each call to solveDKE, so it is called more than once in a parameter scan.
    ! Here we compute all magnetic quantities on the psi and theta grids.
    ! The arrays that must be initialized in this subroutine are:
    ! BHat(Npsi,Ntheta)
    ! dBHatdpsi(Npsi,Ntheta)
    ! dBHatdtheta(Npsi,Ntheta)
    ! JHat(Npsi,Ntheta)
    ! IHat(Npsi)
    ! dIHatdpsi(Npsi)

    implicit none

    PetscScalar, allocatable, dimension(:) :: bs_1D, dbdthetas_1D, oneOverqRbDotGradThetas_1D
    integer :: i

    allocate(BHat(Ntheta,Npsi))
    allocate(dBHatdpsi(Ntheta,Npsi))
    allocate(dBHatdtheta(Ntheta,Npsi))
    allocate(JHat(Ntheta,Npsi))
    allocate(IHat(Npsi))
    allocate(dIHatdpsi(Npsi))

    select case (geometryToUse)
    case (0,1,2)
       ! Simplistic profiles in which magnetic quantities have no radial variation

       allocate(bs_1D(Ntheta))  
       allocate(dbdthetas_1D(Ntheta))
       allocate(oneOverqRbDotGradThetas_1D(Ntheta))
       call computeBs_1D(theta, bs_1D)
       call computedBdthetas_1D(theta, dbdthetas_1D)
       call computeOneOverqRbDotGradThetas_1D(theta, oneOverqRbDotGradThetas_1D)

       do i=1,Npsi
          BHat(:,i) = bs_1D
          dBHatdtheta(:,i) = dbdthetas_1D
          JHat(:,i) = bs_1D / oneOverqRbDotGradThetas_1D / Miller_q
       end do
       dBHatdpsi = 0
       IHat = 1
       dIHatdpsi = 0

    case (3)
       ! EFIT interface would go here.
       ! The arrays that need to be filled are:
       ! BHat(Npsi,Ntheta)
       ! dBHatdpsi(Npsi,Ntheta)
       ! dBHatdtheta(Npsi,Ntheta)
       ! JHat(Npsi,Ntheta)
       ! IHat(Npsi)
       ! dIHatdpsi(Npsi)

    case (4)
       ! Read arrays from hdf5 file
       call h5fopen_f(trim(geometryFilename), H5F_ACC_RDONLY_F, HDF5FileID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error opening profiles input file"
          stop
       end if
   

       ! Get group for geometry data
       write (HDF5Groupname,"(A4,I0,A6,I0)") "Npsi", Npsi, "Ntheta",Ntheta
       call h5gopen_f(HDF5FileID, trim(HDF5Groupname), HDF5GroupID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error opening group for profiles with Npsi=",Npsi,", Ntheta=",Ntheta
          stop
       end if


       ! Specify the dimensions of the arrays to read from the hdf5 file
       arraydims2d(1) = Npsi
       arraydims2d(2) = Ntheta

       ! BHat
       call h5dopen_f(HDF5GroupID, "BHat", HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error opening BHat dataset"
          stop
       end if
       call h5dread_f(HDF5DatasetID, H5T_NATIVE_DOUBLE, BHat, arraydims2d, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error reading BHat"
          stop
       end if
       call h5dclose_f(HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error closing BHat dataset"
          stop
       end if 

       ! dBHatdpsi
       call h5dopen_f(HDF5GroupID, "dBHatdpsi", HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error opening dBHatdpsi dataset"
          stop
       end if
       call h5dread_f(HDF5DatasetID, H5T_NATIVE_DOUBLE, dBHatdpsi, arraydims2d, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error reading dBHatdpsi"
          stop
       end if
       call h5dclose_f(HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error closing dBHatdpsi dataset"
          stop
       end if

       ! dBHatdtheta
       call h5dopen_f(HDF5GroupID, "dBHatdtheta", HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error opening dBHatdtheta dataset"
          stop
       end if
       call h5dread_f(HDF5DatasetID, H5T_NATIVE_DOUBLE, dBHatdtheta, arraydims2d, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error reading dBHatdtheta"
          stop
       end if
       call h5dclose_f(HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error closing dBHatdtheta dataset"
          stop
       end if

       ! JHat
       call h5dopen_f(HDF5GroupID, "JHat", HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error opening JHat dataset"
          stop
       end if
       call h5dread_f(HDF5DatasetID, H5T_NATIVE_DOUBLE, JHat, arraydims2d, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error reading JHat"
          stop
       end if
       call h5dclose_f(HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error closing JHat dataset"
          stop
       end if

       !Read 1D IHat data
       arraydims1d(1) = Npsi
       
       ! IHat
       call h5dopen_f(HDF5GroupID, "IHat", HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error opening IHat dataset"
          stop
       end if
       call h5dread_f(HDF5DatasetID, H5T_NATIVE_DOUBLE, IHat, arraydims1d, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error reading IHat"
          stop
       end if
       call h5dclose_f(HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error closing IHat dataset"
          stop
       end if

       ! dIHatdpsi
       call h5dopen_f(HDF5GroupID, "dIHatdpsi", HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error opening dIHatdpsi dataset"
          stop
       end if
       call h5dread_f(HDF5DatasetID, H5T_NATIVE_DOUBLE, dIHatdpsi, arraydims1d, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error reading dIHatdpsi"
          stop
       end if
       call h5dclose_f(HDF5DatasetID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error closing dIHatdpsi dataset"
          stop
       end if
 
       call h5gclose_f(HDF5GroupID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error closing geometry group"
          stop
       end if

       call h5fclose_f(HDF5FileID, HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error closing geometry file"
          stop
       end if


    case default
       print *,"Error! Invalid setting for geometryToUse"
       stop
    end select

  end subroutine computeMagneticQuantitiesOnGrids

  ! Fill arrays that can be calculated from the magnetic geometry.
  subroutine computeDerivedMagneticQuantities()
    integer :: i

    allocate(VPrimeHat(Npsi))
    allocate(FSABHat2(Npsi))
    allocate(typicalB(Npsi))

    do i=1,Npsi
       VPrimeHat(i) = dot_product(thetaWeights, 1/JHat(:,i))
       FSABHat2(i) = dot_product(thetaWeights, BHat(:,i) * BHat(:,i) / JHat(:,i)) / VPrimeHat(i)
       typicalB(i) = sqrt(FSABHat2(i))
    end do

  end subroutine

  !-----------------------------------------------------------------------------------------
  ! Next are a set of functions needed only for simplistic profiles
  !-----------------------------------------------------------------------------------------

  subroutine computeBs_1D(thetas, bs)
    ! For simplistic profiles, compute the magnitude |B| on the theta grids, 
    ! assuming there is no interesting radial variation.

    implicit none

    PetscScalar, intent(in) :: thetas(:)
    PetscScalar, intent(out) :: bs(:)
    integer :: i

    if (.not. initializedYet) then
       print *,"Error!  No geometry has been initialized yet."
       stop
    end if

    select case (geometryToUse)
    case (0)
       ! Circular concentric flux surfaces
       bs = 1 / (1 + epsil * cos(thetas))
    case (1)
       ! Miller geometry
       do i=1,size(thetas)
          bs(i) = sqrt(BPoloidal(thetas(i))**2 + 1./((RHat(thetas(i)))**2))
       end do
    case (2)
       ! Circular concentric flux surfaces with Boozer poloidal angle
       bs = 1 + epsil * cos(thetas)
    case default
       print *,"Error! Invalid geometry."
       stop
    end select
  end subroutine computeBs_1D

  !-----------------------------------------------------------------------------------------

  subroutine computedBdthetas_1D(thetas, dBdthetas)
    ! For simplistic profiles, compute d |B| / d theta on the theta grid,
    ! assuming there is no interesting radial variation.

    implicit none

    PetscScalar, intent(in) :: thetas(:)
    PetscScalar, intent(out) :: dBdthetas(:)
    PetscScalar, dimension(:), allocatable :: temp
    PetscScalar, allocatable :: spectralDerivative(:,:), thetaFine(:), bs(:), dBdthetaFine(:)
    PetscScalar, allocatable :: d2dtheta2(:,:), weights(:)
    integer :: i, N, multipliedIndex
    integer, parameter :: dBdthetaResolutionMultiplier = 10

    if (.not. initializedYet) then
       print *,"Error!  No geometry has been initialized yet."
       stop
    end if

    select case (geometryToUse)
    case (0)
       ! Circular concentric flux surfaces
       allocate(temp(size(thetas)))
       temp = 1 + epsil * cos(thetas)
       dBdthetas = epsil * sin(thetas) / (temp*temp)
       deallocate(temp)
    case (1)
       ! Miller geometry

       ! It is not worth analytically differentiating the Miller formulae.
       ! Instead, just numerically differentiate b(theta).

       N=size(thetas)*dBdthetaResolutionMultiplier

       allocate(spectralDerivative(N,N))
       allocate(thetaFine(N))
       allocate(d2dtheta2(N,N))
       allocate(weights(N))
       allocate(bs(N))
       allocate(dBdthetaFine(N))
       call uniformDiffMatrices(N, zero, two*pi, 20, thetaFine, weights, spectralDerivative, d2dtheta2)

       call computeBs_1D(thetaFine, bs)
       dBdthetaFine = matmul(spectralDerivative, bs)
       do i=1,size(thetas)
          multipliedIndex = (i-1)*dBdthetaResolutionMultiplier+1
          dBdthetas(i) = dBdthetaFine(multipliedIndex)
          if (abs(thetas(i) - thetaFine(multipliedIndex)) > 1e-10) then
             print *,"Error! The input theta array to computedBdthetas was not of the expected form."
             print *,"thetas:",thetas
             print *,"thetaFine:",thetaFine
             stop
          end if
       end do

       deallocate(spectralDerivative)
       deallocate(thetaFine)
       deallocate(bs)
       deallocate(dBdthetaFine)
    case (2)
       ! Circular concentric flux surfaces with Boozer poloidal angle
       dBdthetas = - epsil * sin(thetas)
    case default
       print *,"Error! Invalid geometry."
       stop
    end select
  end subroutine computedBdthetas_1D

  !-----------------------------------------------------------------------------------------

  subroutine computeOneOverqRbDotGradThetas_1D(thetas, oneOverqRbDotGradThetas)
    ! For simplistic profiles, compute 1/(q R \vect{b} dot grad theta) on the theta grid,
    ! assuming there is no interesting radial variation.

    implicit none

    PetscScalar, intent(in) :: thetas(:)
    PetscScalar, intent(out) :: oneOverqRbDotGradThetas(:)
    integer :: i
    PetscScalar, allocatable :: bs(:)

    if (.not. initializedYet) then
       print *,"Error!  No geometry has been initialized yet."
       stop
    end if

    select case (geometryToUse)
    case (0)
       ! Circular concentric flux surfaces
       oneOverqRbDotGradThetas = [(one, i=1,size(thetas))]
    case (1)
       ! Miller geometry
       allocate(bs(size(thetas)))
       call computeBs_1D(thetas,bs)
       do i=1,size(thetas)
          oneOverqRbDotGradThetas(i) = bs(i) / (Miller_q*BDotGradTheta(thetas(i)));
       end do
       deallocate(bs)
    case (2)
       ! Circular concentric flux surfaces with Boozer poloidal angle
       oneOverqRbDotGradThetas = 1 / (1 + epsil*cos(thetas))
    case default
       print *,"Error! Invalid geometry."
       stop
    end select
  end subroutine computeOneOverqRbDotGradThetas_1D

  !-----------------------------------------------------------------------------------------
  ! Below are a set of functions needed only for Miller geometry
  !-----------------------------------------------------------------------------------------

  function RHat(theta)

    implicit none

    PetscScalar :: theta, RHat

    RHat = 1 + (1/Miller_A)*cos(theta + Miller_x*sin(theta))

  end function RHat

  !-----------------------------------------------------------------------------------------

  function ZHat(theta)

    implicit none

    PetscScalar :: theta, ZHat

    ZHat = (Miller_kappa/Miller_A)*sin(theta)

  end function ZHat

  !-----------------------------------------------------------------------------------------

  function QQIntegrand(theta)

    implicit none

    PetscScalar :: theta, QQIntegrand

    QQIntegrand = ((1+Miller_s_kappa)*sin(theta + Miller_x*sin(theta)) * (1+Miller_x*cos(theta)) * sin(theta) &
         + cos(theta) * (Miller_dRdr + cos(theta + Miller_x *sin(theta)) &
         - Miller_s_delta*sin(theta + Miller_x*sin(theta)) * sin(theta))) / RHat(theta)

  end function QQIntegrand

  !-----------------------------------------------------------------------------------------

  function BPoloidal(theta)

    implicit none

    PetscScalar :: theta, BPoloidal

    BPoloidal = Miller_QQ/(Miller_kappa*Miller_q)*sqrt((sin(theta+Miller_x*sin(theta)) &
         * (1+Miller_x*cos(theta)))**2 + (Miller_kappa*cos(theta))**2) &
         / (RHat(theta) * ( cos(Miller_x*sin(theta)) + Miller_dRdr*cos(theta) + (Miller_s_kappa-Miller_s_delta*cos(theta) &
         + (1+Miller_s_kappa)*Miller_x*cos(theta)) * sin(theta) * sin(theta + Miller_x*sin(theta))))

  end function BPoloidal

  !-----------------------------------------------------------------------------------------

  function BDotGradTheta(theta)

    implicit none

    PetscScalar :: theta, BDotGradTheta

    BDotGradTheta = - Miller_A*Miller_QQ/(Miller_kappa*Miller_q*RHat(theta) * &
         ((1+Miller_s_kappa)*sin(theta + Miller_x*sin(theta)) * (1+Miller_x*cos(theta)) * sin(theta) &
         + cos(theta) * (Miller_dRdr + cos(theta + Miller_x *sin(theta)) &
         - Miller_s_delta*sin(theta + Miller_x*sin(theta)) * sin(theta))))

  end function BDotGradTheta

end module geometry
