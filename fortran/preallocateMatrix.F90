#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

subroutine preallocateMatrix(matrix, whichMatrix, finalMatrix)
  ! whichMatrix = 0 for global, 1 for local.
  ! finalMatrix = 0 for preconditioner matrix, 1 for final matrix.

  use petscmat
  use globalVariables, only: Nx, Nxi, Ntheta, Npsi, Nspecies, Nsources,matrixSize, localMatrixSize, &
       MPIComm, masterProcInSubComm, numProcsInSubComm, PETSCPreallocationStrategy, &
       psiDerivativeScheme, thetaDerivativeScheme, xDerivativeScheme, &
       preconditioner_species, preconditioner_x, preconditioner_x_min_L, &
       preconditioner_psi, preconditioner_theta, preconditioner_xi, &
       lowestEnforcedIpsi, highestEnforcedIpsi, NEnforcedPsi
  use indices

  implicit none

  integer, intent(in) :: whichMatrix, finalMatrix
  integer :: thisThetaDerivativeScheme, thisPsiDerivativeScheme, thisxDerivativeScheme
  Mat :: matrix
  integer :: predictedNNZForEachRowOfTotalMatrix, tempInt
  integer, dimension(:), allocatable :: predictedNNZsForEachRow, predictedNNZsForEachRowDiagonal
  PetscErrorCode :: ierr
  integer :: predictedNNZPerRow_DKE, i, itheta, ipsi, ispecies, ix, index, isources, iextraSources
  integer :: firstRowThisProcOwns, lastRowThisProcOwns, numLocalRows, thisMatrixSize

  MPI_Comm :: MPICommToUse

  if (masterProcInSubComm) then
     print *,"Beginning preallocation for whichMatrix = ",whichMatrix
  end if

  select case (whichMatrix)
  case (0)
     thisMatrixSize = matrixSize
     MPICommToUse = MPIComm
  case (1)
     thisMatrixSize = localMatrixSize
     MPICommToUse = PETSC_COMM_SELF
  case default
     stop "Invalid whichMatrix!"
  end select

  select case (finalMatrix)
  case (0)
     
     select case (preconditioner_theta)
     case (0)
       thisThetaDerivativeScheme = thetaDerivativeScheme
     case (1)
       ! Use 3-point finite difference stencil for preconditioner
       thisThetaDerivativeScheme = 1
     case default
       stop "Unrecognized preconditioner_theta"
     end select

     select case (preconditioner_psi)
     case (0)
       ! Use full coupling
       thisPsiDerivativeScheme = psiDerivativeScheme
     case (1)
       ! Use 'less accurate derivative', i.e. 3-point stencil
       thisPsiDerivativeScheme = 1
     case (2)
       ! Drop ddpsi term
       thisPsiDerivativeScheme = 8
     case (3)
       ! Keep only diagonal of ddpsi. Diagonal is already going to be allocated, so do nothing
       thisPsiDerivativeScheme = 8
     case (4)
       ! Drop all global terms (no ddpsi stencil again?)
       thisPsiDerivativeScheme = 8
     case default
       stop "Unrecognized preconditioner_psi"
     end select

     if (preconditioner_x_min_L == 0) then
       select case (preconditioner_x)
       case (0)
         ! Keep full coupling
         thisXDerivativeScheme = xDerivativeScheme
       case (1)
         ! Drop everything off-diagonal
         thisXDerivativeScheme = 8
       case (2)
         ! Keep only upper-triangular part, just preallocate for full derivative (for now)
         thisXDerivativeScheme = xDerivativeScheme
       case (3)
         ! Keep only tridiagonal terms
         thisXDerivativeScheme = 2
       case (4)
         ! Keep only the diagonal and superdiagonal, preallocate for tridiagonal (for now)
         thisXDerivativeScheme = 2
       case (5)
         ! Use 3-point finite differences
         thisXDerivativeScheme = 2
       case default
         stop "Unrecognized preconditioner_x"
       end select
     else
       ! Only simplify x derivative for some values of L, just preallocate for the full derivative for all L
       thisXDerivativeScheme = xDerivativeScheme
     end if
     
  case (1)
     thisThetaDerivativeScheme = thetaDerivativeScheme
     thisPsiDerivativeScheme = psiDerivativeScheme
     thisXDerivativeScheme = xDerivativeScheme
  case default
     stop "Invalid finalMatrix value"
  end select

  allocate(predictedNNZsForEachRow(thisMatrixSize))
  allocate(predictedNNZsForEachRowDiagonal(thisMatrixSize))

  ! Set predictedNNZPerRow_DKE to the expected number of nonzeros in a row of the kinetic equation block:

  if (finalMatrix==0 .and. preconditioner_xi==1) then
    ! Assumes each source only gives a non-zero for one L? 
    predictedNNZPerRow_DKE = 3 & ! d/dxi term is tridiagonal for the preconditioner
                            + Nsources + NextraSources + 2*NspeciesIndepSources 
  else
    predictedNNZPerRow_DKE = 5 & ! d/dxi term is pentadiagonal
                            + Nsources + NextraSources + 2*NspeciesIndepSources 
  end if

  select case (thisThetaDerivativeScheme)
  case (0)
     ! Spectral collocation
     predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + Ntheta*5-1  ! d/dtheta terms (dense in theta, pentadiagonal in L, -1 since we already counted the diagonal)
  case (1)
     ! 3 point stencil
     predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 3*5-1       ! d/dtheta terms (tridiagonal in theta, pentadiagonal in L, -1 since we already counted the diagonal)
  case (2)
     ! 5 point stencil
     predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 5*5-1       ! d/dtheta terms (pentadiagonal in theta, pentadiagonal in L, -1 since we already counted the diagonal)
  case default
     stop "Invalid thetaDerivativeScheme"
  end select

  if (whichMatrix == 0) then
     select case (thisPsiDerivativeScheme)
     case (1)
        ! 3 point stencil
        predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 3*3-1      ! d/dpsi terms (tridiagonal in psi, tridiagonal in L, -1 since we already counted the diagonal)
     case (2)
        ! 5 point stencil
        predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 5*3-1      ! d/dpsi terms (tridiagonal in psi, tridiagonal in L, -1 since we already counted the diagonal)
     case (8)
        ! no ddpsi term, so do nothing
     case default
        stop "Invalid psiDerivativeScheme"
     end select
  end if

  select case (xDerivativeScheme)
  case (0,2)
     ! Spectral collocation
     predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + Nx*3-1          ! xdot*d/dx terms (dense in x, tridiagonal in L, -1 since we already counted the diagonal)
     if (.not. (finalMatrix==0 .and. preconditioner_species==1) ) then ! if we are building the preconditioner and preconditioner_species=1 then no collisional coupling in the matrix
       predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + Nx*(Nspecies-1) ! collision operator (dense in x, dense in species, -Nx since we already counted the terms diagonal in both x and species.)
     end if
  case (1)
     ! 5 point stencil
     predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 5*3-1      ! xdot*d/dx terms (pentadiagonal in psi, tridiagonal in L, -1 since we already counted the diagonal)
     if (.not. (finalMatrix==0 .and. preconditioner_species==1) ) then ! if we are building the preconditioner and preconditioner_species=1 then no collisional coupling in the matrix
       predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 5*(Nspecies-1) ! collision operator (pentadiagonal in x, dense in species, -5 since we already counted the terms diagonal in both x and species.)
     end if
!!$  case (2)
!!$     ! 3 point stencil, only used for preconditioner at the moment
!!$     predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 3*3-1      ! xdot*d/dx terms (tridiagonal in psi, tridiagonal in L, -1 since we already counted the diagonal)
!!$     if (.not. (finalMatrix==0 .and. preconditioner_species==1) ) then ! if we are building the preconditioner and preconditioner_species=1 then no collisional coupling in the matrix
!!$       predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 3*(Nspecies-1) ! collision operator (tridiagonal in x, dense in species, -3 since we already counted the terms diagonal in both x and species.)
!!$     end if
  case (8)
     ! Drop everything off-diagonal in x for the preconditioner, so do nothing
  case default
     stop "Invalid xDerivativeScheme"
  end select

  ! PETSc gets angry if you request more nonzeros than the matrix size:
  if (predictedNNZPerRow_DKE > thisMatrixSize) then
     predictedNNZPerRow_DKE = thisMatrixSize
  end if

  predictedNNZsForEachRow = predictedNNZPerRow_DKE

  ! The rows for the constraints have more nonzeros:
  if (whichMatrix==0) then
     do isources = 1,Nsources
        do ispecies = 1,Nspecies
           do ipsi = lowestEnforcedIpsi, highestEnforcedIpsi
           
              index = Npsi*localMatrixSize + (ipsi-lowestEnforcedIpsi)*Nspecies*Nsources + (ispecies-1)*Nsources + isources
              predictedNNZsForEachRow(index) = Ntheta*Nx + 1  !+1 for diagonal
           end do
        end do
     end do

     do iextraSources = 1, NextraSources
        do ipsi = lowestEnforcedIpsi, highestEnforcedIpsi
           index = getIndexExtraSources(iextraSources,ipsi)
           ! Would be more accurate to use non-zero species structure element
           predictedNNZsForEachRow(index) = Nspecies*Ntheta*Nx + 1 
        end do
     end do

     do iextraSources = 1, NspeciesIndepSources
        do ipsi = 1, Npsi
           index = getIndexSpeciesIndepSources(iextraSources,ipsi)
           ! Would be more accurate to use non-zero species structure element
           ! factor "2" is from NL
           predictedNNZsForEachRow(index) = 2*Nspecies*Ntheta*Nx + 1
        end do
     end do
  end if

  predictedNNZsForEachRowDiagonal = predictedNNZsForEachRow
  
  select case (PETSCPreallocationStrategy)
  case (0)
     ! 0 = Old method with high estimated number-of-nonzeros.
     ! This method is simpler but does not work consistently.

     predictedNNZForEachRowOfTotalMatrix = predictedNNZPerRow_DKE

     call MatCreateAIJ(MPICommToUse, PETSC_DECIDE, PETSC_DECIDE, thisMatrixSize, thisMatrixSize, &
          predictedNNZsForEachRow, PETSC_NULL_INTEGER, &
          predictedNNZsForEachRow, PETSC_NULL_INTEGER, &
          matrix, ierr)
     
  case (1)
     ! 1 = New method with more precise estimated number-of-nonzeros.
     ! This method is more complicated, but it should use much less memory.
     
     call MatCreate(MPICommToUse, matrix, ierr)
     !call MatSetType(matrix, MATMPIAIJ, ierr)
     call MatSetType(matrix, MATAIJ, ierr)
     
     numLocalRows = PETSC_DECIDE
     call PetscSplitOwnership(MPICommToUse, numLocalRows, thisMatrixSize, ierr)
     
     call MatSetSizes(matrix, numLocalRows, numLocalRows, PETSC_DETERMINE, PETSC_DETERMINE, ierr)
     
     ! We first pre-allocate assuming number-of-nonzeros = 0, because due to a quirk in PETSc,
     ! MatGetOwnershipRange only works after MatXXXSetPreallocation is called:
     if (numProcsInSubComm == 1 .or. whichMatrix==1) then
        call MatSeqAIJSetPreallocation(matrix, 0, PETSC_NULL_INTEGER, ierr)
     else
        call MatMPIAIJSetPreallocation(matrix, 0, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER, ierr)
     end if
     
     call MatGetOwnershipRange(matrix, firstRowThisProcOwns, lastRowThisProcOwns, ierr)
     !print *,"I am proc ",myRank," and I own rows ",firstRowThisProcOwns," to ",lastRowThisProcOwns-1
     
     ! To avoid a PETSc error message, the predicted nnz for each row of the diagonal blocks must be no greater than the # of columns this proc owns:
     ! But we must not lower the predicted nnz for the off-diagonal blocks, because then the total predicted nnz for the row
     ! would be too low.
     tempInt = lastRowThisProcOwns - firstRowThisProcOwns
     do i=firstRowThisProcOwns+1,lastRowThisProcOwns
        if (predictedNNZsForEachRowDiagonal(i) > tempInt) then
           predictedNNZsForEachRowDiagonal(i) = tempInt
        end if
     end do
     
     ! Now, set the real estimated number-of-nonzeros:
     if (numProcsInSubComm == 1 .or. whichMatrix==1) then
        call MatSeqAIJSetPreallocation(matrix, 0, predictedNNZsForEachRow(firstRowThisProcOwns+1:lastRowThisProcOwns), ierr)
     else
        call MatMPIAIJSetPreallocation(matrix, &
             0, predictedNNZsForEachRowDiagonal(firstRowThisProcOwns+1:lastRowThisProcOwns), &
             0, predictedNNZsForEachRow(firstRowThisProcOwns+1:lastRowThisProcOwns), ierr)
     end if
     
  case default
     print *,"Error! Invalid setting for PETSCPreallocationStrategy."
     stop
  end select
  

  ! If any mallocs are required during matrix assembly, do not generate an error:
  call MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
  
  !if (masterProcInSubComm) then
  !   print *,"Done with preallocation for whichMatrix = ",whichMatrix
  !end if

end subroutine preallocateMatrix
