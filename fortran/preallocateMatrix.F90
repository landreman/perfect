#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

subroutine preallocateMatrix(matrix, whichMatrix)
  ! whichMatrix = 0 for global, 1 for local.

  use petscmat
  use globalVariables, only: Nx, Nxi, Ntheta, Npsi, numSpecies, matrixSize, localMatrixSize, &
       MPIComm, masterProcInSubComm, numProcsInSubComm, PETSCPreallocationStrategy, &
       psiDerivativeScheme, thetaDerivativeScheme, xDerivativeScheme

  implicit none

  integer, intent(in) :: whichMatrix
  Mat :: matrix
  integer :: predictedNNZForEachRowOfTotalMatrix, tempInt
  integer, dimension(:), allocatable :: predictedNNZsForEachRow, predictedNNZsForEachRowDiagonal
  PetscErrorCode :: ierr
  integer :: predictedNNZPerRow_DKE, i, itheta, ipsi, ispecies, ix, index
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

  allocate(predictedNNZsForEachRow(thisMatrixSize))
  allocate(predictedNNZsForEachRowDiagonal(thisMatrixSize))

  ! Set predictedNNZPerRow_DKE to the expected number of nonzeros in a row of the kinetic equation block:

  predictedNNZPerRow_DKE = 5 & ! d/dxi term is pentadiagonal
                          + 2  ! particle and heat sources

  select case (thetaDerivativeScheme)
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
     select case (psiDerivativeScheme)
     case (1)
        ! 3 point stencil
        predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 3*3-1      ! d/dpsi terms (tridiagonal in psi, tridiagonal in L, -1 since we already counted the diagonal)
     case (2)
        ! 5 point stencil
        predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 5*3-1      ! d/dpsi terms (tridiagonal in psi, tridiagonal in L, -1 since we already counted the diagonal)
     case default
        stop "Invalid psiDerivativeScheme"
     end select
  end if

  select case (xDerivativeScheme)
  case (0)
     ! Spectral collocation
     predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + Nx*3-1          ! xdot*d/dx terms (dense in x, tridiagonal in L, -1 since we already counted the diagonal)
     predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + Nx*(numSpecies-1) ! collision operator (dense in x, dense in species, -Nx since we already counted the terms diagonal in both x and species.)
  case (1)
     ! 5 point stencil
     predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 5*3-1      ! xdot*d/dx terms (pentadiagonal in psi, tridiagonal in L, -1 since we already counted the diagonal)
     predictedNNZPerRow_DKE = predictedNNZPerRow_DKE + 5*(numSpecies-1) ! collision operator (pentadiagonal in x, dense in species, -5 since we already counted the terms diagonal in both x and species.)
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
     do ispecies = 1,numSpecies
        do ipsi = 1, Npsi
           index = Npsi*localMatrixSize + (ipsi-1)*numSpecies*2 + (ispecies-1)*2 + 1
           predictedNNZsForEachRow(index) = Ntheta*Nx + 1  !+1 for diagonal
           predictedNNZsForEachRow(index+1) = Ntheta*Nx + 1  !+1 for diagonal
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
  !call MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
  
  !if (masterProcInSubComm) then
  !   print *,"Done with preallocation for whichMatrix = ",whichMatrix
  !end if

end subroutine preallocateMatrix