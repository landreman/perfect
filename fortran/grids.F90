module grids

  use globalVariables
  use petscksp
  use petscdmda
  use polynomialDiffMatrices
  use printToStdout
  use xGrid

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

  implicit none

  PetscScalar, dimension(:), allocatable :: psiWeights
  PetscScalar, dimension(:,:), allocatable :: ddpsi_accurate, ddpsiLeft, ddpsiRight, d2dpsi2
  PetscScalar, dimension(:), allocatable :: thetaWeights
  PetscScalar, dimension(:,:), allocatable :: ddthetaToUse
  PetscScalar, dimension(:), allocatable :: x, xWeights, xPotentials
  PetscScalar, dimension(:), allocatable :: x2
  PetscScalar, dimension(:,:), allocatable :: ddxPreconditioner, d2dx2Preconditioner, ddxToUse, d2dx2ToUse
  PetscScalar, dimension(:,:), allocatable :: ddx, d2dx2, ddxPotentials, d2dx2Potentials
  PetscScalar, dimension(:), allocatable :: expx2
  PetscScalar, dimension(:), allocatable :: LegendresOnXiUniform
  PetscScalar, dimension(:,:), allocatable :: ddpsiForPreconditioner
  PetscScalar, dimension(:,:), allocatable :: ddtheta, d2dtheta2
  PetscScalar, dimension(:,:), allocatable :: ddtheta_preconditioner, d2dtheta2_preconditioner
  PetscScalar, dimension(:), allocatable :: LegendresOnXiUniform_m1, LegendresOnXiUniform_m2
  PetscScalar, dimension(:,:), allocatable :: regridPolynomialToUniform
  PetscScalar, dimension(:,:), allocatable :: regridPolynomialToUniformForDiagnostics
  integer :: ipsiMin, ipsiMax, NxPotentials
  logical :: procThatHandlesLeftBoundary, procThatHandlesRightBoundary
  PetscScalar :: xMaxNotTooSmall

  contains

  subroutine createGrids(upwinding)
    PetscScalar, dimension(:,:), allocatable :: localddpsiLeftInterior, localddpsiRightInterior
    PetscScalar, dimension(:), allocatable :: theta_preconditioner, thetaWeights_preconditioner
    PetscScalar, dimension(:), allocatable :: xWeightsPotentials
    PetscScalar, dimension(:,:), allocatable :: dgdx, d2gdx2
    PetscScalar :: xMin, exp_px2, exp_mx2
    logical, intent(in) :: upwinding
    PetscErrorCode :: ierr
    integer :: i, ix, j
    integer :: localNpsi
    integer :: scheme, ipsiMinInterior, ipsiMaxInterior, localNpsiInterior
    PetscScalar :: temp, temp2
    DM :: myDM

        if (forceOddNtheta) then
       if (mod(Ntheta, 2) == 0) then
          Ntheta = Ntheta + 1
       end if
    end if

    call printInputs()

    localMatrixSize = Ntheta * Nxi * Nx * numSpecies
    matrixSize = Npsi * (localMatrixSize + 2*numSpecies)
    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] The matrix is ",matrixSize,"x",matrixSize," elements."
    end if

    didItConverge = integerToRepresentTrue

    ! Validate some input quantities:
    if (preconditioner_xi /= 0 .and. preconditioner_xi /= 1) then
       print *,"Error! preconditioner_xi must be 0 or 1."
       stop
    end if
    if (preconditioner_theta<0 .or. preconditioner_theta > 1) then
       print *,"Error! preconditioner_theta must be 0 or 1."
       stop
    end if
    if (preconditioner_x < 0 .or. preconditioner_x > 5) then
       print *,"Error! preconditioner_x must be in [0,5]."
       stop
    end if
    if (preconditioner_psi < 0 .or. preconditioner_psi > 4) then
       print *,"Error! preconditioner_psi must be between 0 and 4 (inclusive)."
       stop
    end if

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Create grids, integration weights, and differentiation matrices
    !
    ! *******************************************************************************
    ! *******************************************************************************

    ! Assign a range of psi indices to each processor.
    ! This is done by creating a PETSc DM that is not actually used for anything else.
    call DMDACreate1d(MPIComm, DM_BOUNDARY_NONE, Npsi, 1, 0, PETSC_NULL_INTEGER, myDM, ierr)
    call DMDAGetCorners(myDM, ipsiMin, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         localNpsi, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
    call DMDestroy(myDM, ierr) ! dubious
    ! Switch to 1-based indices:
    ipsiMin = ipsiMin + 1

    ipsiMax = ipsiMin+localNpsi-1
    procThatHandlesLeftBoundary = masterProcInSubComm
    if (numProcsInSubComm > Npsi) then
       procThatHandlesRightBoundary = (myRankInSubComm == Npsi)
    else
       procThatHandlesRightBoundary = (ipsiMax == Npsi)
    end if
    if (localNpsi<1 .and. procThatHandlesLeftBoundary) then
       print *,"Error! proc that handles left boundary should have localNpsi .ge. 1. myRank=",myRank
       stop
    end if
    if (localNpsi<1 .and. procThatHandlesRightBoundary) then
       print *,"Error! proc that handles right boundary should have localNpsi .ge. 1. myRank=",myRank
       stop
    end if

    CHKERRQ(ierr)
    print *,"[",myCommunicatorIndex,"] Processor ",myRank," owns psi indices ",ipsiMin," to ",ipsiMax

    ipsiMinInterior = max(2, ipsiMin)
    ipsiMaxInterior = min(Npsi-1, ipsiMax)
    localNpsiInterior = ipsiMaxInterior - ipsiMinInterior + 1

    ! Each processor is responsible for building the rows of the matrix and rhs corresponding
    ! to its ipsiMin:ipsiMax, and each processor is resposible for all columns of the matrix.

    ! *******************************************************************************
    ! Build psi grid, integration weights, and differentiation matrices:
    ! *******************************************************************************

    psiMin = psiMid - psiDiameter/two - widthExtender + leftBoundaryShift
    psiMax = psiMid + psiDiameter/two + widthExtender + rightBoundaryShift

    allocate(psi(Npsi))

    if (Npsi<5) then ! if Npsi<5 then we can do without psi-derivatives, but the simulation must be local
      
      !Sanity checks
      if (.not. makeLocalApproximation) then
        stop "Npsi is less than 5; this only makes sense when makeLocalApproximation is true, but it is false"
      end if
      if ( leftBoundaryScheme/=2 .or. rightBoundaryScheme/=2 ) then
        stop "Are you sure you want leftBoundaryScheme or rightBoundaryScheme other than 2 &
              when running with Npsi<5? It is likely to be inefficient because the local &
              solutions are calculated multiple times for the boundary points."
      end if
      
      ! Just build psi grid
      if (Npsi>1) then
        ! Include points at both psiMin and psiMax:
        psi = [( (psiMax-psiMin)*i/(Npsi-1)+psiMin, i=0,Npsi-1 )]
      else
        ! Unless Npsi=1 so that there is only one point, then:
        psi = psiMid
      end if
    else
      allocate(psiWeights(Npsi))
      allocate(ddpsiForPreconditioner(Npsi,Npsi))
      allocate(ddpsiLeft(Npsi,Npsi))
      allocate(ddpsiRight(Npsi,Npsi))
      allocate(d2dpsi2(Npsi,Npsi))
      allocate(localddpsiLeftInterior(localNpsiInterior,Npsi))
      allocate(localddpsiRightInterior(localNpsiInterior,Npsi))

      ! Build a less-accurate ddpsi which might be used for the preconditioner:
      if (preconditioner_psi == 1 .and. psiDerivativeScheme == 0) then
         print *,"Error! It is presently not compatible to set preconditioner_psi=1 and psiDerivativeScheme=0."
         stop
      end if
      ! centered finite differences, no upwinding, 3-point stencil
      scheme = 2
      call uniformDiffMatrices(Npsi, psiMin, psiMax, scheme, psi, psiWeights, ddpsiForPreconditioner, d2dpsi2)
      ! All of the returned arrays above will be over-written except for ddpsiForPreconditioner

      select case (psiDerivativeScheme)
      case (1)
         ! centered finite differences, 3-point stencil
         scheme = 2
         call uniformDiffMatrices(Npsi, psiMin, psiMax, scheme, psi, psiWeights, ddpsiLeft, d2dpsi2)
      case (2)
         ! centered finite differences, 5-point stencil
         scheme = 12
         call uniformDiffMatrices(Npsi, psiMin, psiMax, scheme, psi, psiWeights, ddpsiLeft, d2dpsi2)
      case default
         if (masterProcInSubComm) then
            print *,"[",myCommunicatorIndex,"] Error! Invalid setting for psiDerivativeScheme"
         end if
         stop
      end select

      allocate(ddpsi_accurate(Npsi,Npsi))

      if (.not. upwinding) then
         ! Copy ddpsiLeft to ddpsiRight
         do i=1,Npsi
            do j=1,Npsi
               ddpsiRight(i,j) = ddpsiLeft(i,j)
            end do
         end do
      end if
      localddpsiLeftInterior = ddpsiLeft(ipsiMinInterior:ipsiMaxInterior,:)
      localddpsiRightInterior = ddpsiRight(ipsiMinInterior:ipsiMaxInterior,:)
    end if


    ! *******************************************************************************
    ! Build theta grid, integration weights, and differentiation matrices:
    ! *******************************************************************************

    select case (thetaDerivativeScheme)
    case (0)
       scheme = 20
    case (1)
       scheme = 0
    case (2)
       scheme = 10
    case default
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] Error! Invalid setting for thetaDerivativeScheme"
       end if
       stop
    end select

    allocate(theta(Ntheta))
    allocate(thetaWeights(Ntheta))
    allocate(ddtheta(Ntheta,Ntheta))
    allocate(ddthetaToUse(Ntheta,Ntheta))
    allocate(d2dtheta2(Ntheta,Ntheta))
    call uniformDiffMatrices(Ntheta, 0d0, two*pi, scheme, theta, thetaWeights, ddtheta, d2dtheta2)

    ! Also make a sparser differentiation matrix for the preconditioner:
    allocate(theta_preconditioner(Ntheta))
    allocate(thetaWeights_preconditioner(Ntheta))
    allocate(ddtheta_preconditioner(Ntheta,Ntheta))
    allocate(d2dtheta2_preconditioner(Ntheta,Ntheta))
    scheme = 0
    call uniformDiffMatrices(Ntheta, 0d0, two*pi, scheme, theta_preconditioner, &
         thetaWeights_preconditioner, ddtheta_preconditioner, d2dtheta2_preconditioner)

    ! Find the theta grid point closest to 0 or 2*pi. We will call it the outboard side.
    temp = 9999
    do i=1,Ntheta
       temp2 = min(abs(theta(i)), abs(theta(i)-2*pi))
       if (temp2 < temp) then
          thetaIndexForOutboard = i
          temp = temp2
       end if
    end do

    ! *******************************************************************************
    ! Build x grids, integration weights, and differentiation matrices.
    ! Also build interpolation matrices to map functions from one x grid to the other.
    ! *******************************************************************************

    allocate(x(Nx))
    allocate(xWeights(Nx))
    allocate(x2(Nx))
    allocate(ddx(Nx,Nx))
    allocate(d2dx2(Nx,Nx))
    allocate(ddxPreconditioner(Nx,Nx))
    allocate(d2dx2Preconditioner(Nx,Nx))
    allocate(ddxToUse(Nx,Nx))
    allocate(d2dx2ToUse(Nx,Nx))

    select case (xDerivativeScheme)
    case (0)
       ! Polynomial spectral collocation
       call makeXGrid(Nx, x, xWeights)
       xWeights = xWeights / exp(-x*x)
       x = x * xScaleFactor
       xWeights = xWeights * xScaleFactor
       call makeXPolynomialDiffMatrices(x,ddx,d2dx2)
       
    case (1)
       ! centered finite differences, 5-point stencil

       ! We cannot have a point right at 0 due to singularities, so make the lowest point slightly positive:
       xMin = xMaxForDistribution/(Nx+one)/5

       ! Rather than differentiate f itself, it works much better to differentiate f*exp(x^2).
       allocate(dgdx(Nx,Nx))
       allocate(d2gdx2(Nx,Nx))

       scheme = 12
       call uniformDiffMatrices(Nx, xMin, xMaxForDistribution, scheme, x, xWeights, dgdx, d2gdx2)

       ! First, right-multiply the matrices by diag(exp(x^2)):
       do ix=1,Nx
          exp_px2 = exp(x(ix) * x(ix))
          dgdx(:,ix) = dgdx(:,ix) * exp_px2
          d2gdx2(:,ix) = d2gdx2(:,ix) * exp_px2
       end do

       
       do ix=1,Nx
          exp_mx2 = exp(-x(ix) * x(ix))

          ddx(ix,:) = exp_mx2 * dgdx(ix,:)
          ddx(ix,ix) = ddx(ix,ix) - 2*x(ix) ! Diagonal

          d2dx2(ix,:) = exp_mx2 * d2gdx2(ix,:) - 4*x(ix)*exp_mx2 * dgdx(ix,:)
          d2dx2(ix,ix) = d2dx2(ix,ix) + 4*x(ix)*x(ix) - 2 ! Diagonal
       end do

       deallocate(dgdx)
       deallocate(d2gdx2)
       
    case default
       print *,"Error! Invalid setting for xDerivativeScheme"
       stop
       
    end select

    xMaxNotTooSmall = max(x(Nx), xMax)
    x2=x*x
    NxPotentials = ceiling(xMaxNotTooSmall*NxPotentialsPerVth)

    allocate(xPotentials(NxPotentials))
    allocate(xWeightsPotentials(NxPotentials))
    allocate(ddxPotentials(NxPotentials, NxPotentials))
    allocate(d2dx2Potentials(NxPotentials, NxPotentials))
    call uniformDiffMatrices(NxPotentials, zero, xMaxNotTooSmall, 12, xPotentials, &
         xWeightsPotentials, ddxPotentials, d2dx2Potentials)

    allocate(regridPolynomialToUniform(NxPotentials, Nx))
    call polynomialInterpolationMatrix(Nx, NxPotentials, x, xPotentials, &
         exp(-x*x), exp(-xPotentials*xPotentials), regridPolynomialToUniform)
    !    allocate(regridUniformToPolynomial(Nx,NxPotentials))
    !    call interpolationMatrix(NxPotentials, Nx, xPotentials, x, regridUniformToPolynomial, -1, 0)

    allocate(expx2(Nx))
    expx2 = exp(-x*x)

  !!$    allocate(erfs(Nx))
  !!$    do i=1,Nx
  !!$#ifdef USE_GSL_ERF
  !!$       call erf(x(i), temp1)
  !!$#else
  !!$       temp1 = erf(x(i))
  !!$#endif
  !!$       erfs(i) = temp1
  !!$    end do

    ! For all preconditioner_x values except 5, we do not need to make
    ! d2dx2Preconditioner, because the xDot term only needs the first derivative.
    ! However, due to the way the collision operator is implemented later,
    ! we do need a d2dx2Preconditioner when preconditioner_x=5.
    ddxPreconditioner = 0
    select case (preconditioner_x)
    case (0)
       ! No simplification in x:
       ddxPreconditioner = ddx
    case (1)
       ! Keep only diagonal terms in x:
       do i=1,Nx
          ddxPreconditioner(i,i) = ddx(i,i)
       end do
    case (2)
       ! Keep only upper-triangular terms in x:
       do i=1,Nx
          do j=i,Nx
             ddxPreconditioner(i,j) = ddx(i,j)
          end do
       end do
    case (3)
       ! Keep only tridiagonal terms in x:
       do i=1,Nx
          do j=1,Nx
             if (abs(i-j) <= 1) then
                ddxPreconditioner(i,j) = ddx(i,j)
             end if
          end do
       end do
    case (4)
       ! Keep only diagonal and super-diagonal in x:
       do i=1,Nx
          ddxPreconditioner(i,i) = ddx(i,i)
       end do
       do i=1,(Nx-1)
          ddxPreconditioner(i,i+1) = ddx(i,i+1)
       end do

    case (5)
       ! Use finite differences with a 3-point stencil:

       print *,"Error! This code is not yet set up for preconditioner_x=5"
       stop
    case default
       print *,"Error! Invalid preconditioner_x"
       stop
    end select


    ! *******************************************************************************
    ! Create some uniform grids in x and xi used for diagnostics
    ! *******************************************************************************

    allocate(xUniform(NxUniform))
    allocate(xiUniform(NxiUniform))
    allocate(LegendresOnXiUniform(NxiUniform))
    allocate(LegendresOnXiUniform_m1(NxiUniform))
     allocate(LegendresOnXiUniform_m2(NxiUniform))

     do i=1,NxUniform
        xUniform(i) = (i-(1d+0))/(NxUniform-1)*xUniformMax
     end do

     do i=1,NxiUniform
        xiUniform(i) = (i-(1d+0))/(NxiUniform-1)*2-1
     end do

     allocate(regridPolynomialToUniformForDiagnostics(NxUniform, Nx))
     call polynomialInterpolationMatrix(Nx, NxUniform, x, xUniform, &
          exp(-x*x), exp(-xUniform*xUniform), regridPolynomialToUniformForDiagnostics)

  end subroutine createGrids

  subroutine deallocateInitializationGridArrays()

    if (Npsi>=5) then
      deallocate(ddpsi_accurate)
    end if
    deallocate(expx2)

  end subroutine deallocateInitializationGridArrays

end module grids
