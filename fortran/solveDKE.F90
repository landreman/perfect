! For compilers that do not include the error function erf(x), the line
! below should be un-commented, and you will need to link to GSL:
!#define USE_GSL_ERF

module solveDKE

  use DKEMatrices
  use DKERhs
  use geometry
  use globalVariables
  use grids
  use matlabOutput
  use moments
  use profiles
  use sparsify

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

  
  implicit none

  private

  Mat :: permutationMatrix
  KSP :: KSPInstance
  Vec, public :: soln

  public :: solveDKEMain

contains

  subroutine solveDKEMain()

    PetscErrorCode :: ierr
    ! integer :: scheme
    ! PetscScalar, dimension(:,:), allocatable :: ddpsiForKTheta
    logical :: upwinding
    ! PetscScalar, dimension(:), allocatable :: dnHatdpsi, detaHatdpsi
    !    PetscScalar, dimension(:), allocatable :: pPerpTermInKThetaFactors
    PetscLogDouble :: time1, time2, startTime
    !    PetscScalar :: scaleFactor, maxXForDistribution

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Start solveDKE()
    !
    ! *******************************************************************************
    ! *******************************************************************************

    call PetscTime(time1, ierr)
    startTime = time1

    if ((.not. isAParallelDirectSolverInstalled) .and. (numProcsInSubComm > 1)) then
       if (masterProcInSubComm) then
          print *,"Error! Neither mumps nor superlu_dist appears to be installed,"
          print *," yet you have asked for a matrix to be distributed across processsors."
       end if
       stop
    end if

    upwinding = .false.
    call createGrids(upwinding)

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Evaluate input physical quantities on the (psi, theta) grids:
    !
    ! *******************************************************************************
    ! *******************************************************************************

    ! This subroutine to initialize the magnetic geometry is in geometry.F90.
    ! It must fill the following arrays:
    ! BHat(Npsi,Ntheta)
    ! dBHatdpsi(Npsi,Ntheta)
    ! dBHatdtheta(Npsi,Ntheta)
    ! JHat(Npsi,Ntheta)
    ! IHat(Npsi)
    ! dIHatdpsi(Npsi)
    call computeMagneticQuantitiesOnGrids()
    
    ! Fill arrays that can be calculated from the magnetic geometry.
    call computeDerivedMagneticQuantities()

    ! This subroutine to initialize the radial profiles is in profiles.F90.
    ! It must fill the following arrays:
    ! PhiHat(Npsi)
    ! dPhiHatdpsi(Npsi)
    ! THats(numSpecies,Npsi)
    ! dTHatdpsis(numSpecies,Npsi)
    ! nHats(numSpecies,Npsi)
    ! dNHatdpsis(numSpecies,NPsi)
    ! etaHats(numSpecies,Npsi)
    ! detaHatdpsis(numSpecies,Npsi)
    call initializeProfiles()

    ! Fill some arrays that can be computed from the radial physics profiles.
    call computeDerivedProfileQuantities()

    ! *********************************************************
    ! *********************************************************
    !
    ! Now build the main matrix, as well as local matrices for 
    ! the left and right boundaries.
    !
    ! *********************************************************
    ! *********************************************************
    call DKECreateMainMatrix(upwinding,time1)
   
    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Create the right-hand side vector
    !
    ! *******************************************************************************
    ! *******************************************************************************
    call DKECreateRhsVector()

    call deallocateInitializationGridArrays()

    ! *********************************************************************************************
    ! If this process handles the left or right boundary, solve the local kinetic equation there:
    ! *********************************************************************************************
    call solveDKEBoundariesLocal(time1)
    ! Finally, assemble the global RHS:
    call VecAssemblyBegin(rhs, ierr)
    call VecAssemblyEnd(rhs, ierr)

    ! ***********************************************************************
    ! ***********************************************************************
    ! 
    !  Permute the rows and columns of the linear system, if desired:
    !
    ! ***********************************************************************
    ! ***********************************************************************
    ! call permuteRowsAndColumns()

    ! ***********************************************************************
    ! ***********************************************************************
    ! 
    !  Solve the main linear system:
    !
    ! ***********************************************************************
    ! ***********************************************************************
    call solveDKEGlobal(time1)

    !**************************************************************************
    !**************************************************************************
    ! 
    !  Calculate moments of the solution:
    !
    !**************************************************************************
    !**************************************************************************
    call calculateMoments(soln)


    ! *********************************************************
    ! Create a PETSc viewer to record output for Matlab
    ! *********************************************************
    if (saveMatlabOutput) then
      call writeMatlabOutput(soln,time1)
    end if

    ! Clean up PETSc objects
    call solveDKECleanUp()

    call PetscTime(time2, ierr)
    elapsedTime = time2 - startTime

    call printOutputs()

  end subroutine solveDKEMain

  ! *********************************************************************************************
  ! If this process handles the left or right boundary, solve the local kinetic equation there:
  ! *********************************************************************************************
  subroutine solveDKEBoundariesLocal(time1)

    PetscErrorCode :: ierr
    integer :: ix, itheta, ipsi, L, index
    integer :: ispecies
    KSP :: KSPBoundary
    PC :: PCBoundary
    KSPConvergedReason :: reason
    PetscScalar :: signOfPsiDot
    Vec :: solnLeft, solnRight
    PetscScalar, pointer :: solnArray(:)
    PetscLogDouble, intent(inout) :: time1
    PetscLogDouble :: time2

    if (procThatHandlesLeftBoundary) then
       ! This process handles the left boundary, so solve the local kinetic equation there.

       call VecAssemblyBegin(rhsLeft, ierr)
       call VecAssemblyEnd(rhsLeft, ierr)
       call VecDuplicate(rhsLeft, solnLeft, ierr)
       CHKERRQ(ierr)
       select case (leftBoundaryScheme)
       case (0)
          print *,"[",myCommunicatorIndex,"] Setting f_1 at left boundary to 0."
          call VecSet(solnLeft, 0d0, ierr)
       case (1)

          call KSPCreate(MPI_COMM_SELF, KSPBoundary, ierr)
       if (useIterativeBoundarySolver) then
             ! Use iterative solver
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
             ! Syntax for PETSc versions up through 3.4:
             call KSPSetOperators(KSPBoundary, leftMatrix, leftPreconditionerMatrix, SAME_PRECONDITIONER, ierr)
#else
             ! Syntax for PETSc version 3.5 and later
             call KSPSetOperators(KSPBoundary, leftMatrix, leftPreconditionerMatrix, ierr)
#endif
             call KSPGetPC(KSPBoundary, PCBoundary, ierr)
             call PCSetType(PCBoundary, PCLU, ierr)
             call KSPSetType(KSPBoundary, KSPBCGSL, ierr)
             call KSPSetTolerances(KSPBoundary, solverTolerance, 1.d-50, &
                  1.d10, PETSC_DEFAULT_INTEGER, ierr)
             call KSPSetFromOptions(KSPBoundary, ierr)
             call KSPMonitorSet(KSPBoundary, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
          else
             ! Direct solver:
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
             ! Syntax for PETSc versions up through 3.4:
             call KSPSetOperators(KSPBoundary, leftMatrix, leftMatrix, SAME_PRECONDITIONER, ierr)
#else
             ! Syntax for PETSc version 3.5 and later
             call KSPSetOperators(KSPBoundary, leftMatrix, leftMatrix, ierr)
#endif
             call KSPGetPC(KSPBoundary, PCBoundary, ierr)
             call PCSetType(PCBoundary, PCLU, ierr)
             call KSPSetType(KSPBoundary, KSPPREONLY, ierr)
             call KSPSetFromOptions(KSPBoundary, ierr)
          end if

          CHKERRQ(ierr)
          print *,"[",myCommunicatorIndex,"] Proc",myRank, &
               " is solving local kinetic equation at left boundary ..."
          if (solveSystem) then
             call KSPSolve(KSPBoundary, rhsLeft, solnLeft, ierr)
          end if
          CHKERRQ(ierr)

          call PetscTime(time2, ierr)
          print *,"[",myCommunicatorIndex,"] Done solving for left boundary.  Time to solve: ", time2-time1, " seconds."
          call PetscTime(time1, ierr)

          if (useIterativeBoundarySolver) then
             call KSPGetConvergedReason(KSPBoundary, reason, ierr)
             if (reason>0) then
                print *,"[",myCommunicatorIndex,"] Solution at left boundary converged!  KSPConvergedReason = ", reason
                !didItConverge = integerToRepresentTrue
             else
                print *,"[",myCommunicatorIndex,"] Solution at left boundary did not converge :(   KSPConvergedReason = ", reason
                didItConverge = integerToRepresentFalse
                if (solveSystem) then
                   stop
                end if
             end if
          else
             !didItConverge = integerToRepresentTrue
          end if

          ! Clean up PETSc objects
          call MatDestroy(leftPreconditionerMatrix, ierr) !dubious
          !call PCDestroy(PCBoundary, ierr) ! Don't need to destroy PC obtained from KSPGetPC?
          call MatDestroy(leftMatrix, ierr) !dubious
          call KSPDestroy(KSPBoundary, ierr)

          ! Before using the local solution, apply the constraints
          call applyBoundaryConstraints(solnLeft, "inner")
       case (2)
          ! Don't do anything
       case default
          print *,"Error! Invalid setting for leftBoundaryScheme"
          stop
       end select

       ! Where trajectories enter the domain, copy solnLeft to the global rhs:
       if (leftBoundaryScheme /= 2) then
         call VecGetArrayF90(solnLeft, solnArray, ierr)
         ipsi = 1
         do ispecies=1,numSpecies
            do itheta=1,Ntheta
               signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                    / (psiAHat*charges(ispecies))
               if (signOfPsiDot > -thresh) then
                  do L=0,(Nxi-1)
                     do ix=1,Nx
                        index = (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + L*Ntheta + itheta
                        call VecSetValue(rhs, index-1, solnArray(index), INSERT_VALUES, ierr)
                     end do
                  end do
               end if
            end do
         end do
         call VecRestoreArrayF90(solnLeft, solnArray, ierr)
       end if

       !       call VecView(solnLeft, PETSC_VIEWER_STDOUT_SELF, ierr)

       call VecDestroy(solnLeft, ierr)
       call VecDestroy(rhsLeft, ierr)
    end if

    if (procThatHandlesRightBoundary) then
       ! This process handles the right boundary, so solve the local kinetic equation there.

       call VecAssemblyBegin(rhsRight, ierr)
       call VecAssemblyEnd(rhsRight, ierr)
       call VecDuplicate(rhsRight, solnRight, ierr)
       CHKERRQ(ierr)
       select case (rightBoundaryScheme)
       case (0)
          print *,"[",myCommunicatorIndex,"] Setting f_1 at right boundary to 0."
          call VecSet(solnRight, 0d0, ierr)
       case (1)

          call KSPCreate(MPI_COMM_SELF, KSPBoundary, ierr)
       if (useIterativeBoundarySolver) then
             ! Use iterative solver
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
             ! Syntax for PETSc versions up through 3.4:
             call KSPSetOperators(KSPBoundary, rightMatrix, rightPreconditionerMatrix, SAME_PRECONDITIONER, ierr)
#else
             ! Syntax for PETSc version 3.5 and later
             call KSPSetOperators(KSPBoundary, rightMatrix, rightPreconditionerMatrix, ierr)
#endif
             call KSPGetPC(KSPBoundary, PCBoundary, ierr)
             call PCSetType(PCBoundary, PCLU, ierr)
             call KSPSetType(KSPBoundary, KSPBCGSL, ierr)
             call KSPSetTolerances(KSPBoundary, solverTolerance, 1.d-50, &
                  1.d10, PETSC_DEFAULT_INTEGER, ierr)
             call KSPSetFromOptions(KSPBoundary, ierr)
             call KSPMonitorSet(KSPBoundary, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
          else
             ! Direct solver:
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
             ! Syntax for PETSc versions up through 3.4:
             call KSPSetOperators(KSPBoundary, rightMatrix, rightMatrix, SAME_PRECONDITIONER, ierr)
#else
             ! Syntax for PETSc version 3.5 and later
             call KSPSetOperators(KSPBoundary, rightMatrix, rightMatrix, ierr)
#endif
             call KSPGetPC(KSPBoundary, PCBoundary, ierr)
             call PCSetType(PCBoundary, PCLU, ierr)
             call KSPSetType(KSPBoundary, KSPPREONLY, ierr)
             call KSPSetFromOptions(KSPBoundary, ierr)
          end if

          CHKERRQ(ierr)
          print *,"[",myCommunicatorIndex,"] Proc",myRank, &
               " is solving local kinetic equation at right boundary ..."
          if (solveSystem) then
             call KSPSolve(KSPBoundary, rhsRight, solnRight, ierr)
          end if
          CHKERRQ(ierr)

          call PetscTime(time2, ierr)
          print *,"[",myCommunicatorIndex,"] Done solving for right boundary.  Time to solve: ", time2-time1, " seconds."
          call PetscTime(time1, ierr)

          if (useIterativeBoundarySolver) then
             call KSPGetConvergedReason(KSPBoundary, reason, ierr)
             if (reason>0) then
                print *,"[",myCommunicatorIndex,"] Solution at right boundary converged!  KSPConvergedReason = ", reason
                !didItConverge = integerToRepresentTrue
             else
                print *,"[",myCommunicatorIndex,"] Solution at right boundary did not converge :(   KSPConvergedReason = ", reason
                didItConverge = integerToRepresentFalse
                if (solveSystem) then
                   stop
                end if
             end if
          else
             !didItConverge = integerToRepresentTrue
          end if

          ! Clean up PETSc objects
          call MatDestroy(rightPreconditionerMatrix, ierr) !dubious
          !call PCDestroy(PCBoundary, ierr) ! Don't need to destroy PC obtained from KSPGetPC?
          call MatDestroy(rightMatrix, ierr) !dubious
          call KSPDestroy(KSPBoundary, ierr)

          ! Before using the local solution, apply the constraints
          call applyBoundaryConstraints(solnRight, "outer")
       case (2)
          ! Don't do anything
       case default
          print *,"Error! Invalid setting for rightBoundaryScheme"
          stop
       end select

       ! Where trajectories enter the domain, copy solnRight to the global rhs:
       if (rightBoundaryScheme /= 2) then
         call VecGetArrayF90(solnRight, solnArray, ierr)
         ipsi = Npsi
         do ispecies = 1,numSpecies
            do itheta=1,Ntheta
               signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                    / (psiAHat*charges(ispecies))
               if (signOfPsiDot < thresh) then
                  do L=0,(Nxi-1)
                     do ix=1,Nx
                        index = (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + L*Ntheta + itheta
                        call VecSetValue(rhs, (ipsi-1)*localMatrixSize + index-1, solnArray(index), INSERT_VALUES, ierr)
                     end do
                  end do
               end if
            end do
         end do
         call VecRestoreArrayF90(solnRight, solnArray, ierr)
       end if

       call VecDestroy(solnRight, ierr)
       call VecDestroy(rhsRight, ierr)
    end if

  end subroutine solveDKEBoundariesLocal

  ! ***********************************************************************
  ! ***********************************************************************
  !
  !  Apply the constraints (flux surface averaged density and pressure of
  !  the perturbed distribution function should vanish) to the solution of
  !  the local DKE before using it as a boundary condition.
  !
  !  The constraints are applied by computing the appropriate moments of the
  !  distribution function found by the solver and subtracting them off: the
  !  resulting vector should still be a solution of the local DKE (apart
  !  from discretisation errors), since the moments on which constraints are
  !  needed are in the null space of the DKE differential operator (indeed
  !  that is why the constraints are needed).
  !
  ! ***********************************************************************
  ! ***********************************************************************
  subroutine applyBoundaryConstraints(soln, whichBoundary)

    Vec, intent(inout) :: soln
    PetscScalar, pointer :: solnArray(:)
    PetscErrorCode :: ierr
    character(len=*), intent(in) :: whichBoundary
    integer :: L, ispecies, ipsi, itheta, ix
    integer, dimension(:), allocatable :: indices
    !PetscScalar, dimension(:), allocatable :: densityIntegralWeights, pressureIntegralWeights
    PetscScalar :: FSALocalDensityPerturbation, FSALocalSecondMomentPerturbation
    PetscScalar, dimension(:), allocatable :: localDensityPerturbation, localSecondMomentPerturbation
    PetscScalar, dimension(:), allocatable :: x2
    PetscScalar :: densityFactors, pressureFactors

    allocate(indices(Nx))
    !allocate(densityIntegralWeights(Nx))
    !allocate(pressureIntegralWeights(Nx))
    allocate(localDensityPerturbation(Ntheta))
    allocate(localSecondMomentPerturbation(Ntheta))
    allocate(x2(Nx))

    call VecGetArrayF90(soln, solnArray, ierr)

    x2 = x*x

    L = 0
    if (whichBoundary == "inner") then
      ipsi = 1
    else if (whichBoundary == "outer") then
      ipsi = Npsi
    else
      stop "Invalid option given to whichBoundary"
    end if
    do ispecies=1,numSpecies
      ! Compute the density and 'second' moments (i.e. integrals of 1 and (x^2-3/2) times g) of the solution
      do itheta=1,Ntheta
        indices = (ispecies-1)*Nx*Nxi*Ntheta &
                  + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + itheta
        ! Removed densityFactors and pressureFactors that are present in moments.F90 (pretty sure they are just normalizations that are not needed here).
        localDensityPerturbation(itheta) = dot_product(xWeights, x2 * solnArray(indices))
        localSecondMomentPerturbation(itheta) = dot_product(xWeights, x2*(x2-1.5d0) * solnArray(indices))
      end do
      ! Take the flux surface averages
      FSALocalDensityPerturbation = dot_product(thetaWeights, localDensityPerturbation/JHat(:,ipsi)) / VPrimeHat(ipsi)
      FSALocalSecondMomentPerturbation = dot_product(thetaWeights, localSecondMomentPerturbation/JHat(:,ipsi)) / VPrimeHat(ipsi)
      !! print *,"ispecies=",ispecies,"ipsi=",ipsi,"density=",FSALocalDensityPerturbation, &
      !!        "secondmoment=",FSALocalSecondMomentPerturbation
      ! Subtract out the moments to set the flux surface averages to zero
      do itheta=1,Ntheta
        indices = (ispecies-1)*Nx*Nxi*Ntheta &
                  + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + itheta
        solnArray(indices) = solnArray(indices) - FSALocalDensityPerturbation*4d0/sqrt(pi)*exp(-x2)
        solnArray(indices) = solnArray(indices) - FSALocalSecondMomentPerturbation &
                                                  *8d0/3d0/sqrt(pi)*(x2-1.5d0)*exp(-x2)
      end do
    end do

    !! ! Re-calculate density and pressure perturbation flux surface averages for debugging
    !! do ispecies=1,numSpecies
    !!   densityFactors = Delta*4*pi*THats(ispecies,ipsi)*sqrtTHats(ispecies,ipsi) &
    !!            / (nHats(ispecies,ipsi)*masses(ispecies)*sqrt(masses(ispecies)))
    !!   pressureFactors = Delta*8*pi/(three*nHats(ispecies,ipsi)) * ((THats(ispecies,ipsi)/masses(ispecies)) ** (1.5d+0))
    !!   ! Compute the flux surface average of the density and pressure moments of the solution
    !!   do itheta=1,Ntheta
    !!     indices = (ispecies-1)*Nx*Nxi*Ntheta &
    !!               + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + itheta
    !!     ! Removed densityFactors and pressureFactors that are present in moments.F90 (pretty sure they are just normalizations that are not needed here).
    !!     localDensityPerturbation(itheta) = dot_product(xWeights, x2 * solnArray(indices))*densityFactors
    !!     localSecondMomentPerturbation(itheta) = dot_product(xWeights, x2*(x2) * solnArray(indices))*pressureFactors
    !!   end do
    !!   FSALocalDensityPerturbation = dot_product(thetaWeights, localDensityPerturbation/JHat(:,ipsi)) / VPrimeHat(ipsi)
    !!   FSALocalSecondMomentPerturbation = dot_product(thetaWeights, localSecondMomentPerturbation/JHat(:,ipsi)) / VPrimeHat(ipsi)
    !!   print *,"check flux surface averages are now zero:"
    !!   print *,"ispecies=",ispecies,"ipsi=",ipsi,"density=",FSALocalDensityPerturbation,"pressure=",FSALocalSecondMomentPerturbation
    !!   print *,"densityPerturbation=",localDensityPerturbation
    !! end do

    call VecRestoreArrayF90(soln, solnArray, ierr)

  end subroutine applyBoundaryConstraints

  ! ***********************************************************************
  ! ***********************************************************************
  ! 
  !  Solve the main linear system:
  !
  ! ***********************************************************************
  ! ***********************************************************************
  subroutine solveDKEGlobal(time1)

    PetscErrorCode :: ierr
    integer :: NNZMain, NNZPreconditioner
    PC :: preconditionerContext
    KSPConvergedReason :: reason
    Vec :: tempVec
    double precision :: myMatInfo(MAT_INFO_SIZE)
    PetscLogDouble, intent(inout) :: time1
    PetscLogDouble :: time2

    call KSPCreate(MPIComm, KSPInstance, ierr)
    CHKERRQ(ierr)

    if (useIterativeSolver) then
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
       ! Syntax for PETSc versions up through 3.4:
       call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, SAME_PRECONDITIONER, ierr)
#else
       ! Syntax for PETSc version 3.5 and later
       call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, ierr)
#endif
       CHKERRQ(ierr)
       call KSPGetPC(KSPInstance, preconditionerContext, ierr)
       CHKERRQ(ierr)
       call PCSetType(preconditionerContext, PCLU, ierr)
       CHKERRQ(ierr)
       call KSPSetType(KSPInstance, KSPBCGSL, ierr)
       CHKERRQ(ierr)
       call KSPSetTolerances(KSPInstance, solverTolerance, 1.d-50, &
           1.d10, PETSC_DEFAULT_INTEGER, ierr)
       CHKERRQ(ierr)
       call KSPSetFromOptions(KSPInstance, ierr)
       CHKERRQ(ierr)
       call KSPMonitorSet(KSPInstance, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
    else
       ! Direct solver:
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
       ! Syntax for PETSc versions up through 3.4:
       call KSPSetOperators(KSPInstance, matrix, matrix, SAME_PRECONDITIONER, ierr)
#else
       ! Syntax for PETSc version 3.5 and later
       call KSPSetOperators(KSPInstance, matrix, matrix, ierr)
#endif
       CHKERRQ(ierr)
       call KSPGetPC(KSPInstance, preconditionerContext, ierr)
       CHKERRQ(ierr)
       call PCSetType(preconditionerContext, PCLU, ierr)
       CHKERRQ(ierr)
       call KSPSetType(KSPInstance, KSPPREONLY, ierr)
       CHKERRQ(ierr)
       call KSPSetFromOptions(KSPInstance, ierr)
       !        call KSPSetTolerances(KSPInstance, 1d-5, PETSC_DEFAULT_REAL, &
       !          PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
       CHKERRQ(ierr)
    end if

    ! The line below is remarked out to force use of either superlu or mumps,
    ! since PETSc's built-in solver doesn't seem to work when there are zeros on the diagonal.
    ! I may eventually want to handle this problem more elegantly...
    !    if (numProcsInSubComm > 1) then
    if (.true.) then
       select case (whichParallelSolverToFactorPreconditioner)
       case (1)
          call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERMUMPS, ierr)
          if (masterProcInSubComm) then
             print *,"[",myCommunicatorIndex,"] Using mumps to factorize the preconditioner."
          end if
       case (2)
          call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERSUPERLU_DIST, ierr)
          if (masterProcInSubComm) then
             print *,"[",myCommunicatorIndex,"] Using superlu_dist to factorize the preconditioner."
          end if
       case default
          if (masterProcInSubComm) then
             print *,"Error: Invalid setting for whichParallelSolverToFactorPreconditioner"
             stop
          end if
       end select
    else
       if (masterProcInSubComm) then
          print *,"[",myCommunicatorIndex,"] Using PETSc's serial sparse direct solver to factorize the preconditioner."
       end if
    end if

    call MatGetInfo(matrix, MAT_GLOBAL_SUM, myMatInfo, ierr)
    NNZMain = nint(myMatInfo(MAT_INFO_NZ_USED))
    if (useIterativeSolver) then
       call MatGetInfo(preconditionerMatrix, MAT_GLOBAL_SUM, myMatInfo, ierr)
       NNZPreconditioner = nint(myMatInfo(MAT_INFO_NZ_USED))
    end if
    if (masterProcInSubComm) then
       if (useIterativeSolver) then
          print *,"[",myCommunicatorIndex,"] # of nonzeros in matrix:",NNZMain,", # of nonzeros in preconditioner:", &
               NNZPreconditioner,", ratio:",((one*NNZPreconditioner)/NNZMain)
       else
          print *,"[",myCommunicatorIndex,"] # of nonzeros in matrix:",NNZMain
       end if
    end if

    call VecDuplicate(rhs, soln, ierr)
    CHKERRQ(ierr)
    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] Beginning the main solve.  This could take a while ..."
    end if

    if (solveSystem) then
       call KSPSolve(KSPInstance, rhs, soln, ierr)
    end if
    CHKERRQ(ierr)

    call PetscTime(time2, ierr)
    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] Done with the main solve.  Time to solve: ", time2-time1, " seconds."
    end if
    call PetscTime(time1, ierr)

    if (useIterativeSolver) then
       call KSPGetConvergedReason(KSPInstance, reason, ierr)
       if (reason>0) then
          if (masterProcInSubComm) then
             print *,"[",myCommunicatorIndex,"] Converged!  KSPConvergedReason = ", reason
          end if
          !didItConverge = integerToRepresentTrue
       else
          if (masterProcInSubComm) then
             print *,"[",myCommunicatorIndex,"] Did not converge :(   KSPConvergedReason = ", reason
          end if
          didItConverge = integerToRepresentFalse
       end if
    else
       !didItConverge = integerToRepresentTrue
    end if

    if (layout /= 0) then
       call MatMult(permutationMatrix, soln, tempVec, ierr)
       call VecDestroy(soln, ierr)
       soln = tempVec
    end if

  end subroutine solveDKEGlobal

  ! ***********************************************************************
  ! ***********************************************************************
  ! 
  !  Permute the rows and columns of the linear system, if desired:
  !
  ! ***********************************************************************
  ! ***********************************************************************
  subroutine permuteRowsAndColumns()

    PetscErrorCode :: ierr
    integer :: i, ix, itheta, ipsi
    integer :: newIndex, oldIndex
    integer :: ixi
    Mat :: tempMat
    Vec :: tempVec

    if (layout /= 0) then
       call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
            1, PETSC_NULL_INTEGER, 1, PETSC_NULL_INTEGER, permutationMatrix, ierr)
       CHKERRQ(ierr)
 
       if (masterProcInSubComm) then
          select case (layout)
          case (1)
             ! Main section of the global matrix:
             do ipsi=1,Npsi
                do ix=1,Nx
                   do ixi=1,Nxi
                      do itheta=1,Ntheta
                         oldIndex = (ipsi-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
                         newIndex = (ix-1)*Npsi*Nxi*Ntheta + (ipsi-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
                         call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
                      end do
                   end do
                end do
             end do
 
             ! Extra rows/columns for sources and constraints:
             do i=1,2
                do ipsi=1,Npsi
                   oldIndex = Npsi*Nx*Nxi*Ntheta + (i-1)*Npsi + ipsi
                   newIndex = oldIndex
                   call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
                end do
             end do
 
          case (2)
             ! Main section of the global matrix:
             do ipsi=1,Npsi
                do ix=1,Nx
                   do ixi=1,Nxi
                      do itheta=1,Ntheta
                         oldIndex = (ipsi-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
                         newIndex = (ipsi-1)*(Nx*Nxi*Ntheta+2) + (ix-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
                         call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
                      end do
                   end do
                end do
             end do
 
             ! Extra rows/columns for sources and constraints:
             do i=1,2
                do ipsi=1,Npsi
                   oldIndex = Npsi*Nx*Nxi*Ntheta + (i-1)*Npsi + ipsi
                   newIndex = (ipsi-1)*(Nx*Nxi*Ntheta+2) + Nx*Nxi*Ntheta + i
                   call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
                end do
             end do
 
 
          case default
             print *,"Error! Invalid setting for layout"
             stop
          end select
       end if
 
       call MatAssemblyBegin(permutationMatrix, MAT_FINAL_ASSEMBLY, ierr)
       call MatAssemblyEnd(permutationMatrix, MAT_FINAL_ASSEMBLY, ierr)
 
       call MatPtAP(matrix, permutationMatrix, MAT_INITIAL_MATRIX, one, tempMat, ierr)
       call MatDestroy(matrix, ierr)
       matrix = tempMat
 
       CHKERRQ(ierr)
       call MatPtAP(preconditionerMatrix, permutationMatrix, MAT_INITIAL_MATRIX, one, tempMat, ierr)
       CHKERRQ(ierr)
       call MatDestroy(preconditionerMatrix, ierr)
       CHKERRQ(ierr)
       preconditionerMatrix = tempMat
       CHKERRQ(ierr)
 
       ! I get an error if I don't "create" tempVec using the statement below, even though
       ! it was not necessary to initialize tempMat in the same manner just above. This seems
       ! like an inconsistency in PETSc...
       call VecDuplicate(rhs, tempVec, ierr)
       CHKERRQ(ierr)
       call MatMultTranspose(permutationMatrix, rhs, tempVec, ierr)
       CHKERRQ(ierr)
       call VecDestroy(rhs, ierr)
       CHKERRQ(ierr)
       rhs = tempVec
       CHKERRQ(ierr)
    end if

  end subroutine permuteRowsAndColumns

  ! Clean up PETSc objects
  subroutine solveDKECleanUp()

    PetscErrorCode :: ierr

    call VecDestroy(rhs, ierr)
    call VecDestroy(soln, ierr)
    CHKERRQ(ierr)
    call MatDestroy(matrix, ierr)
    if (useIterativeSolver) then
       call MatDestroy(preconditionerMatrix, ierr)
    end if
    call KSPDestroy(KSPInstance, ierr)
    CHKERRQ(ierr)

  end subroutine solveDKECleanUp

end module solveDKE
