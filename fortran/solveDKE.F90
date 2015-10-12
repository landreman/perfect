! For compilers that do not include the error function erf(x), the line
! below should be un-commented, and you will need to link to GSL:
!#define USE_GSL_ERF

#include <finclude/petsckspdef.h>
#include <finclude/petscdmdadef.h>

#include "PETScVersions.F90"
  
subroutine solveDKE()

  use DKEMatrices
  use DKERhs
  use geometry
  use globalVariables
  use grids
  use profiles
  use sparsify

  implicit none

  PetscErrorCode :: ierr
  Vec :: soln, solnLeft, solnRight, solnOnProc0
  PetscViewer MatlabOutput
  integer :: i, ix, itheta, ipsi, L, index
  integer :: ispecies
  integer :: scheme
  integer, dimension(:), allocatable :: indices
  PetscScalar :: speciesFactor
  PetscScalar, dimension(:,:), allocatable :: ddpsiForKTheta
  PetscScalar :: signOfPsiDot, xPartOfSource
  logical :: upwinding
  integer :: constraintToAdd
  PetscScalar, dimension(:), allocatable :: dnHatdpsi, detaHatdpsi
  PetscScalar, dimension(:), allocatable :: densityFactors, densityIntegralWeights
  PetscScalar, dimension(:), allocatable :: flowFactors, flowIntegralWeights
  PetscScalar, dimension(:), allocatable :: pressureFactors, pressureIntegralWeights
  PetscScalar, dimension(:), allocatable :: particleFluxFactors, particleFluxIntegralWeights
  PetscScalar, dimension(:), allocatable :: momentumFluxFactors, momentumFluxIntegralWeights
  PetscScalar, dimension(:), allocatable :: heatFluxFactors, heatFluxIntegralWeights
  !    PetscScalar, dimension(:), allocatable :: pPerpTermInKThetaFactors
  PetscLogDouble :: time1, time2, startTime
  KSP :: KSPInstance, KSPBoundary
  PC :: preconditionerContext, PCBoundary
  KSPConvergedReason :: reason
  PetscScalar, pointer :: solnArray(:)
  VecScatter :: VecScatterContext
  integer :: ixi, newIndex, oldIndex
  PetscScalar :: stuffToAdd
  Mat :: permutationMatrix, tempMat
  Vec :: tempVec
  double precision :: myMatInfo(MAT_INFO_SIZE)
  integer :: NNZMain, NNZPreconditioner
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
   if (procThatHandlesLeftBoundary) then
      ! This process handles the left boundary, so solve the local kinetic equation there.

      call VecAssemblyBegin(rhsLeft, ierr)
      call VecAssemblyEnd(rhsLeft, ierr)

      call KSPCreate(MPI_COMM_SELF, KSPBoundary, ierr)
      if (useIterativeSolver) then
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
         call KSPSetTolerances(KSPBoundary, solverTolerance, PETSC_DEFAULT_REAL, &
              PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
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

      call VecDuplicate(rhsLeft, solnLeft, ierr)
      CHKERRQ(ierr)
      select case (leftBoundaryScheme)
      case (0)
         print *,"[",myCommunicatorIndex,"] Setting f_1 at left boundary to 0."
         call VecSet(solnLeft, 0, ierr)
      case (1)
         print *,"[",myCommunicatorIndex,"] Proc",myRank, &
              " is solving local kinetic equation at left boundary ..."
         if (solveSystem) then
            call KSPSolve(KSPBoundary, rhsLeft, solnLeft, ierr)
         end if
         CHKERRQ(ierr)

         call PetscTime(time2, ierr)
         print *,"[",myCommunicatorIndex,"] Done solving for left boundary.  Time to solve: ", time2-time1, " seconds."
         call PetscTime(time1, ierr)

         if (useIterativeSolver) then
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
      case default
         print *,"Error! Invalid setting for leftBoundaryScheme"
         stop
      end select

      ! Where trajectories enter the domain, copy solnLeft to the global rhs:
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

      !       call VecView(solnLeft, PETSC_VIEWER_STDOUT_SELF, ierr)

      call KSPDestroy(KSPBoundary, ierr)
      call VecDestroy(solnLeft, ierr)
      call VecDestroy(rhsLeft, ierr)
   end if

   if (procThatHandlesRightBoundary) then
      ! This process handles the right boundary, so solve the local kinetic equation there.

      call VecAssemblyBegin(rhsRight, ierr)
      call VecAssemblyEnd(rhsRight, ierr)

      call KSPCreate(MPI_COMM_SELF, KSPBoundary, ierr)
      if (useIterativeSolver) then
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
         call KSPSetTolerances(KSPBoundary, solverTolerance, PETSC_DEFAULT_REAL, &
              PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
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

      call VecDuplicate(rhsRight, solnRight, ierr)
      CHKERRQ(ierr)
      select case (rightBoundaryScheme)
      case (0)
         print *,"[",myCommunicatorIndex,"] Setting f_1 at right boundary to 0."
         call VecSet(solnRight, 0, ierr)
      case (1)
         print *,"[",myCommunicatorIndex,"] Proc",myRank, &
              " is solving local kinetic equation at right boundary ..."
         if (solveSystem) then
            call KSPSolve(KSPBoundary, rhsRight, solnRight, ierr)
         end if
         CHKERRQ(ierr)

         call PetscTime(time2, ierr)
         print *,"[",myCommunicatorIndex,"] Done solving for right boundary.  Time to solve: ", time2-time1, " seconds."
         call PetscTime(time1, ierr)

         if (useIterativeSolver) then
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
      case default
         print *,"Error! Invalid setting for rightBoundaryScheme"
         stop
      end select

      ! Where trajectories enter the domain, copy solnRight to the global rhs:
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

      !call PCDestroy(PCBoundary, ierr)
      call KSPDestroy(KSPBoundary, ierr)
      call VecDestroy(solnRight, ierr)
      call VecDestroy(rhsRight, ierr)
   end if

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

 !!$    if (layout /= 0) then
 !!$       call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
 !!$            1, PETSC_NULL_INTEGER, 1, PETSC_NULL_INTEGER, permutationMatrix, ierr)
 !!$       CHKERRQ(ierr)
 !!$
 !!$       if (masterProcInSubComm) then
 !!$          select case (layout)
 !!$          case (1)
 !!$             ! Main section of the global matrix:
 !!$             do ipsi=1,Npsi
 !!$                do ix=1,Nx
 !!$                   do ixi=1,Nxi
 !!$                      do itheta=1,Ntheta
 !!$                         oldIndex = (ipsi-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
 !!$                         newIndex = (ix-1)*Npsi*Nxi*Ntheta + (ipsi-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
 !!$                         call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
 !!$                      end do
 !!$                   end do
 !!$                end do
 !!$             end do
 !!$
 !!$             ! Extra rows/columns for sources and constraints:
 !!$             do i=1,2
 !!$                do ipsi=1,Npsi
 !!$                   oldIndex = Npsi*Nx*Nxi*Ntheta + (i-1)*Npsi + ipsi
 !!$                   newIndex = oldIndex
 !!$                   call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
 !!$                end do
 !!$             end do
 !!$
 !!$          case (2)
 !!$             ! Main section of the global matrix:
 !!$             do ipsi=1,Npsi
 !!$                do ix=1,Nx
 !!$                   do ixi=1,Nxi
 !!$                      do itheta=1,Ntheta
 !!$                         oldIndex = (ipsi-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
 !!$                         newIndex = (ipsi-1)*(Nx*Nxi*Ntheta+2) + (ix-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
 !!$                         call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
 !!$                      end do
 !!$                   end do
 !!$                end do
 !!$             end do
 !!$
 !!$             ! Extra rows/columns for sources and constraints:
 !!$             do i=1,2
 !!$                do ipsi=1,Npsi
 !!$                   oldIndex = Npsi*Nx*Nxi*Ntheta + (i-1)*Npsi + ipsi
 !!$                   newIndex = (ipsi-1)*(Nx*Nxi*Ntheta+2) + Nx*Nxi*Ntheta + i
 !!$                   call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
 !!$                end do
 !!$             end do
 !!$
 !!$
 !!$          case default
 !!$             print *,"Error! Invalid setting for layout"
 !!$             stop
 !!$          end select
 !!$       end if
 !!$
 !!$       call MatAssemblyBegin(permutationMatrix, MAT_FINAL_ASSEMBLY, ierr)
 !!$       call MatAssemblyEnd(permutationMatrix, MAT_FINAL_ASSEMBLY, ierr)
 !!$
 !!$       call MatPtAP(matrix, permutationMatrix, MAT_INITIAL_MATRIX, one, tempMat, ierr)
 !!$       call MatDestroy(matrix, ierr)
 !!$       matrix = tempMat
 !!$
 !!$       CHKERRQ(ierr)
 !!$       call MatPtAP(preconditionerMatrix, permutationMatrix, MAT_INITIAL_MATRIX, one, tempMat, ierr)
 !!$       CHKERRQ(ierr)
 !!$       call MatDestroy(preconditionerMatrix, ierr)
 !!$       CHKERRQ(ierr)
 !!$       preconditionerMatrix = tempMat
 !!$       CHKERRQ(ierr)
 !!$
 !!$       ! I get an error if I don't "create" tempVec using the statement below, even though
 !!$       ! it was not necessary to initialize tempMat in the same manner just above. This seems
 !!$       ! like an inconsistency in PETSc...
 !!$       call VecDuplicate(rhs, tempVec, ierr)
 !!$       CHKERRQ(ierr)
 !!$       call MatMultTranspose(permutationMatrix, rhs, tempVec, ierr)
 !!$       CHKERRQ(ierr)
 !!$       call VecDestroy(rhs, ierr)
 !!$       CHKERRQ(ierr)
 !!$       rhs = tempVec
 !!$       CHKERRQ(ierr)
 !!$    end if

   ! ***********************************************************************
   ! ***********************************************************************
   ! 
   !  Solve the main linear system:
   !
   ! ***********************************************************************
   ! ***********************************************************************

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
      call KSPSetTolerances(KSPInstance, solverTolerance, PETSC_DEFAULT_REAL, &
           PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
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

   !**************************************************************************
   !**************************************************************************
   ! 
   !  Calculate moments of the solution:
   !
   !**************************************************************************
   !**************************************************************************

   ! First, send the entire solution vector to the master process:
   call VecScatterCreateToZero(soln, VecScatterContext, solnOnProc0, ierr)
   CHKERRQ(ierr)
   call VecScatterBegin(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
   call VecScatterEnd(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
   call VecScatterDestroy(VecScatterContext, ierr)
   CHKERRQ(ierr)


   if (masterProcInSubComm) then
      ! All computation of moments of the distribution function is then done on the master process:

      allocate(particleSourceProfile(numSpecies,Npsi))
      allocate(heatSourceProfile(numSpecies,Npsi))

      allocate(densityPerturbation(numSpecies,Ntheta,Npsi))
      allocate(flow(numSpecies,Ntheta,Npsi))
      allocate(kPar(numSpecies,Ntheta,Npsi))
      allocate(pressurePerturbation(numSpecies,Ntheta,Npsi))
      allocate(particleFluxBeforeThetaIntegral(numSpecies,Ntheta,Npsi))
      allocate(momentumFluxBeforeThetaIntegral(numSpecies,Ntheta,Npsi))
      allocate(heatFluxBeforeThetaIntegral(numSpecies,Ntheta,Npsi))

      allocate(FSADensityPerturbation(numSpecies,Npsi))
      allocate(kParOutboard(numSpecies,Npsi))
      allocate(kParInboard(numSpecies,Npsi))
      allocate(FSAKPar(numSpecies,Npsi))
      allocate(flowOutboard(numSpecies,Npsi))
      allocate(flowInboard(numSpecies,Npsi))
      allocate(FSABFlow(numSpecies,Npsi))
      allocate(FSAPressurePerturbation(numSpecies,Npsi))
      allocate(particleFlux(numSpecies,Npsi))
      allocate(momentumFlux(numSpecies,Npsi))
      allocate(heatFlux(numSpecies,Npsi))

      allocate(deltaFOutboard(numSpecies,Npsi,NxUniform,NxiUniform))
      allocate(fullFOutboard(numSpecies,Npsi,NxUniform,NxiUniform))

      !       allocate(LHSOfKParEquation(Npsi))
      !       allocate(kThetaWith3PointStencil(Ntheta,Npsi))
      !       allocate(kThetaWith5PointStencil(Ntheta,Npsi))
      !       allocate(kThetaOutboardWith3PointStencil(Npsi))
      !       allocate(kThetaInboardWith3PointStencil(Npsi))
      !       allocate(kThetaOutboardWith5PointStencil(Npsi))
      !       allocate(kThetaInboardWith5PointStencil(Npsi))
      !       allocate(PhiTermInKTheta(Ntheta,Npsi))
      !       allocate(pPerpTermInKThetaWith3PointStencil(Ntheta,Npsi))
      !       allocate(pPerpTermInKThetaWith5PointStencil(Ntheta,Npsi))
      !       allocate(pPerpTermInKThetaBeforePsiDerivative(Ntheta,Npsi))
      !       allocate(ddpsiForKTheta(Npsi,Npsi))

      allocate(flowFactors(Npsi))
      allocate(densityFactors(Npsi))
      allocate(pressureFactors(Npsi))
      allocate(particleFluxFactors(Npsi))
      allocate(momentumFluxFactors(Npsi))
      allocate(heatFluxFactors(Npsi))
      !       allocate(pPerpTermInKThetaFactors(Npsi))

      allocate(flowIntegralWeights(Nx))
      allocate(densityIntegralWeights(Nx))
      allocate(pressureIntegralWeights(Nx))
      allocate(particleFluxIntegralWeights(Nx))
      allocate(momentumFluxIntegralWeights(Nx))
      allocate(heatFluxIntegralWeights(Nx))

      densityIntegralWeights = x*x
      flowIntegralWeights = x*x*x
      pressureIntegralWeights = x*x*x*x
      particleFluxIntegralWeights = x*x*x*x
      momentumFluxIntegralWeights = x*x*x*x*x
      heatFluxIntegralWeights = x*x*x*x*x*x

      allocate(indices(Nx))

      ! Convert the PETSc vector into a normal Fortran array:
      call VecGetArrayF90(solnOnProc0, solnArray, ierr)
      CHKERRQ(ierr)

      do ispecies = 1,numSpecies
         densityFactors = Delta*4*pi*THats(ispecies,:)*sqrtTHats(ispecies,:) &
              / (nHats(ispecies,:)*masses(ispecies)*sqrt(masses(ispecies)))
         flowFactors = 4*pi/(three*nHats(ispecies,:)) * ((THats(ispecies,:)/masses(ispecies)) ** 2)
         pressureFactors = Delta*8*pi/(three*nHats(ispecies,:)) * ((THats(ispecies,:)/masses(ispecies)) ** (1.5d+0))
         particleFluxFactors = -masses(ispecies) / charges(ispecies) * IHat * ((THats(ispecies,:)/masses(ispecies)) ** (5/two))
         momentumFluxFactors = -masses(ispecies) / charges(ispecies) * IHat*IHat * ((THats(ispecies,:)/masses(ispecies)) ** 3)
         heatFluxFactors = -masses(ispecies) / charges(ispecies) * THats(ispecies,:) &
              * IHat * ((THats(ispecies,:)/masses(ispecies)) ** (5/two))
         !       pPerpTermInKThetaFactors = THat ** (5/two)

         ! The final elements of the solution vector correspond to the source profiles:
         do ipsi=1,Npsi
            particleSourceProfile(ispecies,ipsi) = solnArray(localMatrixSize*Npsi + (ipsi-1)*numSpecies*2 + (ispecies-1)*2 + 1)
            heatSourceProfile(ispecies,ipsi) = solnArray(localMatrixSize*Npsi + (ipsi-1)*numSpecies*2 + (ispecies-1)*2 + 2)
         end do

         L = 0
         do ipsi=1,Npsi
            do itheta=1,Ntheta
               indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                    + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + itheta

               densityPerturbation(ispecies,itheta,ipsi) = dot_product(xWeights, densityIntegralWeights * solnArray(indices)) &
                    * densityFactors(ipsi)

               pressurePerturbation(ispecies,itheta,ipsi) = dot_product(xWeights, pressureIntegralWeights * solnArray(indices)) &
                    * pressureFactors(ipsi)

               particleFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = (8/three) * particleFluxFactors(ipsi) &
                    * dot_product(xWeights, particleFluxIntegralWeights * solnArray(indices))

               heatFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = (8/three) * heatFluxFactors(ipsi) &
                    * dot_product(xWeights, heatFluxIntegralWeights * solnArray(indices))

               !             pPerpTermInKThetaBeforePsiDerivative(itheta,ipsi) = &
               !                  (4/three) * pPerpTermInKThetaFactors(ipsi) &
               !                  * dot_product(xWeights, pressureIntegralWeights * solnArray(indices))

            end do
         end do

         L = 1
         do ipsi=1,Npsi
            do itheta=1,Ntheta
               indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                    + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + itheta

               flow(ispecies,itheta,ipsi) = dot_product(xWeights, flowIntegralWeights * solnArray(indices)) &
                    * flowFactors(ipsi)

               momentumFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = ((16d+0)/15) * momentumFluxFactors(ipsi) &
                    * dot_product(xWeights, momentumFluxIntegralWeights * solnArray(indices))

            end do
         end do

         L = 2
         do ipsi=1,Npsi
            do itheta=1,Ntheta
               indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                    + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + itheta

               particleFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = particleFluxBeforeThetaIntegral(ispecies,itheta,ipsi) &
                    + (four/15) * particleFluxFactors(ipsi) &
                    * dot_product(xWeights, particleFluxIntegralWeights * solnArray(indices))

               heatFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = heatFluxBeforeThetaIntegral(ispecies,itheta,ipsi) &
                    + (four/15) * heatFluxFactors(ipsi) &
                    * dot_product(xWeights, heatFluxIntegralWeights * solnArray(indices))

               !             pPerpTermInKThetaBeforePsiDerivative(itheta,ipsi) = &
               !                  pPerpTermInKThetaBeforePsiDerivative(itheta,ipsi) &
               !                  - ((4d+0)/15) * pPerpTermInKThetaFactors(ipsi) &
               !                  * dot_product(xWeights, pressureIntegralWeights * solnArray(indices))

            end do
         end do

         L = 3
         do ipsi=1,Npsi
            do itheta=1,Ntheta
               indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                    + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + itheta

               momentumFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = momentumFluxBeforeThetaIntegral(ispecies,itheta,ipsi) &
                    + (four/35) * momentumFluxFactors(ipsi) &
                    * dot_product(xWeights, momentumFluxIntegralWeights * solnArray(indices))

            end do
         end do

         do itheta=1,Ntheta
            kpar(ispecies,itheta,:) = FSABHat2/(BHat(itheta,:)*BHat(itheta,:)*dTHatdpsis(ispecies,:)) &
                 * (2*charges(ispecies)*psiAHat*BHat(itheta,:)/IHat*flow(ispecies,itheta,:) + dTHatdpsis(ispecies,:) &
                 + THats(ispecies,:)/nHats(ispecies,:)*dnHatdpsis(ispecies,:) + 2*charges(ispecies)*omega/Delta*dPhiHatdpsi)
         end do

         particleFluxBeforeThetaIntegral(ispecies,:,:) = particleFluxBeforeThetaIntegral(ispecies,:,:) &
              * dBHatdtheta / (BHat * BHat * BHat)

         momentumFluxBeforeThetaIntegral(ispecies,:,:) = momentumFluxBeforeThetaIntegral(ispecies,:,:) &
              * dBHatdtheta / (BHat * BHat * BHat * BHat)

         heatFluxBeforeThetaIntegral(ispecies,:,:) = heatFluxBeforeThetaIntegral(ispecies,:,:) &
              * dBHatdtheta / (BHat * BHat * BHat)

 !!$          if (psiDerivativeScheme == 0) then
 !!$             ddpsiForKTheta = ddpsiLeft
 !!$          else
 !!$             ! centered finite differences, no upwinding, 3-point stencil
 !!$             scheme = 2
 !!$             call uniformDiffMatrices(Npsi, psiMin, psiMax, scheme, psi, psiWeights, ddpsiForKTheta, d2dpsi2)
 !!$          end if
 !!$          !       pPerpTermInKThetaWith3PointStencil = transpose(matmul(ddpsiForKTheta, &
 !!$          !            transpose(pPerpTermInKThetaBeforePsiDerivative)))
 !!$
 !!$          if (psiDerivativeScheme == 0) then
 !!$             ddpsiForKTheta = ddpsiLeft
 !!$          else
 !!$             ! centered finite differences, no upwinding, 5-point stencil
 !!$             scheme = 12
 !!$             call uniformDiffMatrices(Npsi, psiMin, psiMax, scheme, psi, psiWeights, ddpsiForKTheta, d2dpsi2)
 !!$          end if
 !!$          !       pPerpTermInKThetaWith5PointStencil = transpose(matmul(ddpsiForKTheta, &
 !!$          !            transpose(pPerpTermInKThetaBeforePsiDerivative)))

         do ipsi = 1,Npsi
            particleFlux(ispecies,ipsi) = dot_product(thetaWeights, particleFluxBeforeThetaIntegral(ispecies,:,ipsi))
            momentumFlux(ispecies,ipsi) = dot_product(thetaWeights, momentumFluxBeforeThetaIntegral(ispecies,:,ipsi))
            heatFlux(ispecies,ipsi) = dot_product(thetaWeights, heatFluxBeforeThetaIntegral(ispecies,:,ipsi))

            !          pPerpTermInKThetaWith3PointStencil(:,ipsi) = pPerpTermInKThetaWith3PointStencil(:,ipsi) &
            !               * 2 * pi * delta * FSABHat2(ipsi) / (nHat(ipsi) * BHat(:,ipsi) * BHat(:,ipsi) * dTHatdpsi(ipsi))
            !
            !          pPerpTermInKThetaWith5PointStencil(:,ipsi) = pPerpTermInKThetaWith5PointStencil(:,ipsi) &
            !               * 2 * pi * delta * FSABHat2(ipsi) / (nHat(ipsi) * BHat(:,ipsi) * BHat(:,ipsi) * dTHatdpsi(ipsi))
            !
            !          PhiTermInKTheta(:,ipsi) = densityPerturbation(:,ipsi) * 2 * omega * dphidpsi(ipsi) * FSABHat2(ipsi) &
            !               / (delta * BHat(:,ipsi) * BHat(:,ipsi) * dTHatdpsi(ipsi))

            FSADensityPerturbation(ispecies,ipsi) = dot_product(thetaWeights, &
                 densityPerturbation(ispecies,:,ipsi)/JHat(:,ipsi)) / VPrimeHat(ipsi)

            FSABFlow(ispecies,ipsi) = dot_product(thetaWeights, flow(ispecies,:,ipsi)*BHat(:,ipsi)/JHat(:,ipsi)) / VPrimeHat(ipsi)

            FSAkPar(ispecies,ipsi) = dot_product(thetaWeights, kPar(ispecies,:,ipsi)/JHat(:,ipsi)) / VPrimeHat(ipsi)

            FSAPressurePerturbation(ispecies,ipsi) = dot_product(thetaWeights, &
                 pressurePerturbation(ispecies,:,ipsi)/JHat(:,ipsi)) / VPrimeHat(ipsi)
         end do

         !       LHSOfKParEquation = FSAKPar - omega*delta/(psiAHat*psiAHat) * IHat*IHat * dphidpsi / FSABHat2 &
         !            * matmul(ddpsiLeft,FSAKPar)

         !       kThetaWith3PointStencil = kpar + PhiTermInKTheta + pPerpTermInKThetaWith3PointStencil
         !       kThetaWith5PointStencil = kpar + PhiTermInKTheta + pPerpTermInKThetaWith5PointStencil

         kParOutboard(ispecies,:) = kPar(ispecies,1,:)
         flowOutboard(ispecies,:) = flow(ispecies,1,:)
         !       kThetaOutboardWith3PointStencil = kThetaWith3PointStencil(1,:)
         !       kThetaOutboardWith5PointStencil = kThetaWith5PointStencil(1,:)
         if (mod(Ntheta,2)==0) then
            kParInboard(ispecies,:) = kPar(ispecies,Ntheta/2+1,:)
            flowInboard(ispecies,:) = flow(ispecies,Ntheta/2+1,:)
            !          kThetaInboardWith3PointStencil = kThetaWith3PointStencil(Ntheta/2+1,:)
            !          kThetaInboardWith5PointStencil = kThetaWith5PointStencil(Ntheta/2+1,:)
         else
            index = (Ntheta+1)/2
            kParInboard(ispecies,:) = (kPar(ispecies,index,:) + kPar(ispecies,index+1,:))*oneHalf
            flowInboard(ispecies,:) = (flow(ispecies,index,:) + flow(ispecies,index+1,:))*oneHalf
            !          kThetaInboardWith3PointStencil = (kThetaWith3PointStencil(index,:) + kThetaWith3PointStencil(index+1,:))*oneHalf
            !          kThetaInboardWith5PointStencil = (kThetaWith5PointStencil(index,:) + kThetaWith5PointStencil(index+1,:))*oneHalf
         end if

      end do

      LegendresOnXiUniform_m1 = 1
      deltaFOutboard = 0
      do L = 0,(Nxi-1)
         ! Recursively evaluate Legendre polynomials on a uniform grid in xi.
         ! The results will be used to map the distribution function from the modal discretization
         ! to a uniform grid.
         if (L == 0) then
            LegendresOnXiUniform = 1
         elseif (L == 1) then
            LegendresOnXiUniform = xiUniform
         else
            LegendresOnXiUniform_m2 = LegendresOnXiUniform_m1
            LegendresOnXiUniform_m1 = LegendresOnXiUniform
            LegendresOnXiUniform = ((2*L-1)*xiUniform * LegendresOnXiUniform_m1 - (L-1)*LegendresOnXiUniform_m2)/L
        end if
        
        do ispecies=1,numSpecies
           do ipsi = 1,Npsi
              indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                   + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + thetaIndexForOutboard

              do ixi = 1,NxiUniform
                 deltaFOutboard(ispecies,ipsi,:,ixi) = deltaFOutboard(ispecies,ipsi,:,ixi) + &
                      LegendresOnXiUniform(ixi) * matmul(regridPolynomialToUniformForDiagnostics, solnArray(indices))
              end do
           end do
        end do
     end do

     ! Now build the full f from delta f:
     do ispecies = 1,numSpecies
        do ipsi = 1,Npsi
           speciesFactor = nHats(ispecies,ipsi)/(pi*sqrtpi*((THats(ispecies,ipsi)/masses(ispecies)) ** (1.5d+0)))
           do ixi = 1,NxiUniform
              do ix = 1,NxUniform
                 fullFOutboard(ispecies,ipsi,ix,ixi) = deltaFOutboard(ispecies,ipsi,ix,ixi) * Delta &
                      + speciesFactor * exp(-xUniform(ix)*xUniform(ix)) ! This line adds the Maxwellian F_M
              end do
           end do
        end do
     end do

     deallocate(indices)

     call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)
     CHKERRQ(ierr)
  end if


  ! *********************************************************
  ! Create a PETSc viewer to record output
  ! *********************************************************

  if (saveMatlabOutput) then
     call PetscViewerASCIIOpen(MPIComm, &
          & MatlabOutputFilename,&
          & MatlabOutput, ierr)
     CHKERRQ(ierr)
     call PetscViewerSetFormat(MatlabOutput, PETSC_VIEWER_ASCII_MATLAB, ierr)
     CHKERRQ(ierr)

     call PetscObjectSetName(rhs, "rhs", ierr)
     CHKERRQ(ierr)
     call VecView(rhs, MatlabOutput, ierr)
     CHKERRQ(ierr)
     call PetscObjectSetName(soln, "soln", ierr)
     CHKERRQ(ierr)
     call VecView(soln, MatlabOutput, ierr)
     CHKERRQ(ierr)

     if (useIterativeSolver) then
        call PetscObjectSetName(preconditionerMatrix, "preconditionerMatrix", ierr)
        CHKERRQ(ierr)
        call MatView(preconditionerMatrix, MatlabOutput, ierr)
        CHKERRQ(ierr)
     end if
     call PetscObjectSetName(matrix, "matrix", ierr)
     CHKERRQ(ierr)
     call MatView(matrix, MatlabOutput, ierr)
     CHKERRQ(ierr)

     call PetscTime(time2, ierr)
     if (masterProcInSubComm) then
        print *,"[",myCommunicatorIndex,"] Time to write output: ", time2-time1, " seconds."
     end if
     call PetscTime(time1, ierr)

     call PetscViewerDestroy(MatlabOutput, ierr)
     CHKERRQ(ierr)
  end if

  !    call VecDestroy(rhs, ierr)
  call VecDestroy(soln, ierr)
  CHKERRQ(ierr)
  call MatDestroy(matrix, ierr)
  if (useIterativeSolver) then
     call MatDestroy(preconditionerMatrix, ierr)
  end if
  call KSPDestroy(KSPInstance, ierr)
  CHKERRQ(ierr)


  call PetscTime(time2, ierr)
  elapsedTime = time2 - startTime

  call printOutputs()

end subroutine solveDKE

