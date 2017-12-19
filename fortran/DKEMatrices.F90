module DKEMatrices

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#include <finclude/petscdmdadef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#include <petsc/finclude/petscdmdadef.h>
#endif
  
  use globalVariables
  use grids
  use indices
  use petscksp
  use petscdmda
  use sparsify
  use sourcesConstraints

  implicit none

  private

  public :: DKECreateMainMatrix

  Mat, public :: matrix, preconditionerMatrix
  Mat, public :: leftMatrix, leftPreconditionerMatrix
  Mat, public :: rightMatrix, rightPreconditionerMatrix

contains
  
  subroutine DKECreateMainMatrix(upwinding,time1)

    ! *********************************************************
    ! *********************************************************
    !
    ! Now build the main matrix, as well as local matrices for 
    ! the left and right boundaries.
    !
    ! *********************************************************
    ! *********************************************************

    ! PetscErrorCode :: ierr
    ! PetscScalar, dimension(:), allocatable :: xAndThetaPartOfConstraint
    PetscScalar :: sqrtMass
    PetscScalar :: signOfPsiDot
    logical :: makeLocalApproximationOriginal
    logical, intent(in) :: upwinding
    !integer :: i, j, ix, itheta, ipsi, L, index
    integer :: scheme
    integer :: whichMatrix, whichMatrixMin
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscLogDouble, intent(inout) :: time1
    PetscLogDouble :: time2

    !predictedNNZForEachRowOfTotalMatrix = 5*(Nx + 5*3 + 5*5 + 3*Nx + 5 +2)
    !predictedNNZForEachRowOfTotalMatrixBoundary = (2*Ntheta + 2 + Nx)

    makeLocalApproximationOriginal = makeLocalApproximation

    if (useIterativeSolver) then
       whichMatrixMin = 0
    else
       whichMatrixMin = 1
    end if
    do whichMatrix = whichMatrixMin,1
       ! When whichMatrix = 0, build the preconditioner.
       ! When whichMatrix = 1, build the final matrix.

       makeLocalApproximation = makeLocalApproximationOriginal
       if (whichMatrix==0 .and. preconditioner_psi==4) then
          makeLocalApproximation = .true.
       end if

       call preallocateMatrices(whichMatrix)

       call setDiagonalsToZero()
       
       ! *********************************************************
       ! Add the streaming and mirror terms which persist even 
       ! in the (delta,omega)=(0,0) limit:
       ! *********************************************************
       call streamingAndMirrorTerms0(whichMatrix)
       
       ! *********************************************************
       ! Add the streaming and mirror terms which vanish
       ! in the (delta,omega)=(0,0) limit:
       ! *********************************************************
       call streamingAndMirrorTerms1(whichMatrix)
       
       ! *********************************************************
       ! Add the collisionless d/dx term:
       ! *********************************************************
       call collisionlessDdx(whichMatrix)
       
       ! *********************************************************
       ! Add the collisionless d/dpsi term:
       ! *********************************************************
       call collisionlessDdpsi(whichMatrix,upwinding)
       
       ! *********************************************************
       ! Add the collision operator
       ! *********************************************************
       call collisionOperator(whichMatrix)
       
       ! *******************************************************************************
       ! Put a 1 on the matrix diagonal where appropriate to enforce the radial boundary condition
       ! *******************************************************************************
       call radialBoundaryConditionDiagonal()
       
       ! *******************************************************************************
       ! Enforce boundary condition at x=0 if needed (if there is a point there)
       ! *******************************************************************************
       call xBoundaryCondition(whichMatrix)

       ! *******************************************************************************
       ! Add sources:
       ! *******************************************************************************
       call sources()
       
       ! *******************************************************************************
       ! Add constraints:
       ! *******************************************************************************
       call constraints()
       
       ! *******************************************************************************
       ! Done inserting values into the matrices.
       ! Now finalize the matrices:
       ! *******************************************************************************
       call finalizeMatrices(whichMatrix,time1)

    end do

    makeLocalApproximation = makeLocalApproximationOriginal

  ! *******************************************************************************
  ! *******************************************************************************
  ! Delete this next section eventually
  ! *******************************************************************************
  ! *******************************************************************************
 !!$       call PetscViewerASCIIOpen(MPIComm, &
 !!$            & "boundaries.m",&
 !!$            & MatlabOutput, ierr)
 !!$       CHKERRQ(ierr)
 !!$       call PetscViewerSetFormat(MatlabOutput, PETSC_VIEWER_ASCII_MATLAB, ierr)
 !!$       CHKERRQ(ierr)
 !!$
 !!$       !call PetscObjectSetName(leftPreconditionerMatrix, "leftPreconditionerMatrix", ierr)
 !!$       !call MatView(leftPreconditionerMatrix, MatlabOutput, ierr)
 !!$       call PetscObjectSetName(leftMatrix, "leftMatrix", ierr)
 !!$       call MatView(leftMatrix, MatlabOutput, ierr)
 !!$
 !!$       !call PetscObjectSetName(rightPreconditionerMatrix, "rightPreconditionerMatrix", ierr)
 !!$       !call MatView(rightPreconditionerMatrix, MatlabOutput, ierr)
 !!$       call PetscObjectSetName(rightMatrix, "rightMatrix", ierr)
 !!$       call MatView(rightMatrix, MatlabOutput, ierr)
 !!$
 !!$       call PetscViewerDestroy(MatlabOutput, ierr)

  end subroutine DKECreateMainMatrix

  subroutine preallocateMatrices(finalMatrix)

    ! finalMatrix==0 for the preconditioner, finalMatrix==1 for the final matrix
    integer, intent(in) :: finalMatrix

    ! In  preallocateMatrix, the last parameter is 0 for the global matrix, or 1 for the local matrices.
    call preallocateMatrix(matrix, 0, finalMatrix)
    if (procThatHandlesLeftBoundary .and. boundaryScheme /= 2  .and. boundaryScheme /= 3) then
       call preallocateMatrix(leftMatrix, 1, finalMatrix)
    end if
    if (procThatHandlesRightBoundary .and. boundaryScheme /= 1  .and. boundaryScheme /= 3) then
       call preallocateMatrix(rightMatrix, 1, finalMatrix)
    end if

  end subroutine preallocateMatrices

  subroutine setDiagonalsToZero()

    ! Sometimes PETSc complains if any of the diagonal elements are not set.
    ! Therefore, set the entire diagonal to 0 to be safe.

    PetscErrorCode :: ierr
    integer :: i

    if (masterProcInSubComm) then
       do i=1,matrixSize
          call MatSetValue(matrix, i-1, i-1, zero, ADD_VALUES, ierr)
       end do
    end if
    if (procThatHandlesLeftBoundary .and. boundaryScheme /= 2  .and. boundaryScheme /= 3) then
       do i=1,localMatrixSize
          call MatSetValue(leftMatrix, i-1, i-1, zero, ADD_VALUES, ierr)
       end do
    end if
    if (procThatHandlesRightBoundary .and. boundaryScheme /= 1  .and. boundaryScheme /= 3) then
       do i=1,localMatrixSize
          call MatSetValue(rightMatrix, i-1, i-1, zero, ADD_VALUES, ierr)
       end do
    end if

  end subroutine setDiagonalsToZero

  subroutine streamingAndMirrorTerms0(whichMatrix)

    ! *********************************************************
    ! Add the streaming and mirror terms which persist even 
    ! in the (delta,omega)=(0,0) limit:
    ! *********************************************************

    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    PetscScalar, dimension(:), allocatable :: thetaPartOfMirrorTerm
    PetscScalar, dimension(:,:), allocatable :: thetaPartMatrix
    PetscScalar, dimension(:,:), allocatable :: thetaPartOfStreamingTerm
    integer, dimension(:), allocatable :: localRowIndices, localColIndices, globalRowIndices, globalColIndices
    integer :: i, ix, itheta, ipsi, L, ixMin
    integer :: ispecies
    integer :: rowIndexArray(1)
    PetscScalar :: signOfPsiDot

    if (whichMatrix==1 .or. preconditioner_theta==0) then
       ddthetaToUse = ddtheta
    else
       ddthetaToUse = ddtheta_preconditioner
    end if

    if (pointAtX0) then
       ixMin = 2
    else
       ixMin = 1
    end if

    allocate(thetaPartOfMirrorTerm(Ntheta))
    allocate(thetaPartMatrix(Ntheta,Ntheta))
    allocate(thetaPartOfStreamingTerm(Ntheta,Ntheta))

    allocate(localRowIndices(Ntheta))
    allocate(localColIndices(Ntheta))
    allocate(globalRowIndices(Ntheta))
    allocate(globalColIndices(Ntheta))
    do ispecies = 1,Nspecies
       do ipsi = ipsiMin, ipsiMax
          thetaPartOfMirrorTerm = -oneHalf * sqrtTHats(ispecies,ipsi) &
               * JHat(:,ipsi) * dBHatdtheta(:,ipsi) &
               / (BHat(:,ipsi) * BHat(:,ipsi))
          do itheta=1,Ntheta
             thetaPartOfStreamingTerm(itheta,:) = JHat(itheta,ipsi)*sqrtTHats(ispecies,ipsi) &
                  / BHat(itheta,ipsi)*ddthetaToUse(itheta,:)
          end do
          do ix=ixMin,Nx
             do L=0,Nxi_for_x(ix)-1
                localRowIndices = [(getIndex(ispecies,ix,L,itheta,1), itheta=1,Ntheta)]
                globalRowIndices = [(getIndex(ispecies,ix,L,itheta,ipsi), itheta=1,Ntheta)]
                if (L < Nxi_for_x(ix)-1) then
                   ! Super-diagonal in L:
                   localColIndices = [(getIndex(ispecies,ix,L+1,itheta,1), itheta=1,Ntheta)]
                   globalColIndices = [(getIndex(ispecies,ix,L+1,itheta,ipsi), itheta=1,Ntheta)]

                   ! Add streaming term:
                   thetaPartMatrix = x(ix)*(L+1)/(two*L+3)*thetaPartOfStreamingTerm

                   ! Add mirror term:
                   ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                   do itheta=1,Ntheta
                      thetaPartMatrix(itheta,itheta) = x(ix)*(L+1)*(L+2)/(two*L+3)*thetaPartOfMirrorTerm(itheta)
                   end do

                   ! Petsc uses a transposed format relative to Fortran:
                   thetaPartMatrix = transpose(thetaPartMatrix)
                   ! Put values in matrix, 
                   if (ipsi==1 .and. boundaryScheme /= 3 .and. leftBoundaryScheme /= 2 .and. boundaryScheme /= 2) then
                      call MatSetValuesSparse(leftMatrix, Ntheta, localRowIndices, Ntheta, localColIndices, &
                           thetaPartMatrix, ADD_VALUES, ierr)
                      do itheta=1,Ntheta
                         signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                              / (psiAHat*charges(ispecies))
                         if (signOfPsiDot < -thresh .and. boundaryScheme == 0) then
                            ! do not enforce boundary conditions
                            ! add default matrix
                            call MatSetValuesSparse(matrix, 1, globalRowIndices(itheta), Ntheta, globalColIndices, &
                                 thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                         end if
                      end do
                   elseif (ipsi==Npsi .and. boundaryScheme /= 3 .and. rightBoundaryScheme /= 2 .and. boundaryScheme /= 1) then
                      call MatSetValuesSparse(rightMatrix, Ntheta, localRowIndices, Ntheta, localColIndices, &
                           thetaPartMatrix, ADD_VALUES, ierr)
                      do itheta=1,Ntheta
                         signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                              / (psiAHat*charges(ispecies))
                         if (signOfPsiDot > thresh .and. boundaryScheme == 0) then
                            rowIndexArray = globalRowIndices(itheta)
                            call MatSetValuesSparse(matrix, 1, rowIndexArray, &
                                 Ntheta, globalColIndices, thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                         end if
                      end do
                   else
                      call MatSetValuesSparse(matrix, Ntheta, globalRowIndices, Ntheta, &
                           globalColIndices, &
                           thetaPartMatrix, ADD_VALUES, ierr)
                   end if
                end if
                if (L>0) then
                   ! Sub-diagonal in L:
                   localColIndices = [(getIndex(ispecies,ix,L-1,itheta,1), itheta=1,Ntheta)]
                   globalColIndices = [(getIndex(ispecies,ix,L-1,itheta,ipsi), itheta=1,Ntheta)]

                   ! Streaming term
                   thetaPartMatrix = x(ix)*L/(two*L-1)*thetaPartOfStreamingTerm

                   ! Mirror term
                   ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                   do itheta=1,Ntheta
                      thetaPartMatrix(itheta,itheta) = -x(ix)*L*(L-1)/(two*L-1)*thetaPartOfMirrorTerm(itheta)
                   end do

                   thetaPartMatrix = transpose(thetaPartMatrix)

                   ! Put values in matrix, noting that Petsc uses a transposed format relative to Fortran
                   if (ipsi==1 .and. boundaryScheme /= 3 .and. leftBoundaryScheme /= 2 .and. boundaryScheme /= 2) then
                      call MatSetValuesSparse(leftMatrix, Ntheta, localRowIndices, Ntheta, localColIndices, &
                           thetaPartMatrix, ADD_VALUES, ierr)
                      do itheta=1,Ntheta
                         signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi)&
                              / (psiAHat*charges(ispecies))
                         if (signOfPsiDot < -thresh .and. boundaryScheme == 0 ) then
                            call MatSetValuesSparse(matrix, 1, globalRowIndices(itheta), Ntheta, globalColIndices, &
                                 thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                         end if
                      end do
                   elseif (ipsi==Npsi .and. boundaryScheme /= 3 .and. rightBoundaryScheme /= 2 .and. boundaryScheme /= 1) then
                      call MatSetValuesSparse(rightMatrix, Ntheta, localRowIndices, Ntheta, localColIndices, &
                           thetaPartMatrix, ADD_VALUES, ierr)
                      do itheta=1,Ntheta
                         signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                              / (psiAHat*charges(ispecies))
                         if (signOfPsiDot > thresh .and. boundaryScheme == 0) then
                            rowIndexArray = globalRowIndices(itheta)
                            call MatSetValuesSparse(matrix, 1, rowIndexArray, &
                                 Ntheta, globalColIndices, thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                         end if
                      end do
                   else
                      call MatSetValuesSparse(matrix, Ntheta, globalRowIndices, &
                           Ntheta, globalColIndices, thetaPartMatrix, ADD_VALUES, ierr)
                   end if
                end if
             end do
          end do
       end do
    end do
    deallocate(localRowIndices)
    deallocate(localColIndices)
    deallocate(globalRowIndices)
    deallocate(globalColIndices)
    deallocate(thetaPartOfStreamingTerm)
    deallocate(thetaPartMatrix)
    deallocate(thetaPartOfMirrorTerm)

  end subroutine streamingAndMirrorTerms0

  subroutine streamingAndMirrorTerms1(whichMatrix)

    ! *********************************************************
    ! Add the streaming and mirror terms which vanish
    ! in the (delta,omega)=(0,0) limit:
    ! *********************************************************

    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    PetscScalar, dimension(:), allocatable :: thetaPartOfMirrorTerm
    PetscScalar, dimension(:,:), allocatable :: thetaPartMatrix
    PetscScalar, dimension(:,:), allocatable :: spatialPartOfStreamingTermDiagonal1
    PetscScalar, dimension(:,:), allocatable :: spatialPartOfStreamingTermDiagonal2
    PetscScalar, dimension(:,:), allocatable :: spatialPartOfStreamingTermDiagonal3
    PetscScalar, dimension(:,:), allocatable :: spatialPartOfStreamingTermOffDiagonal
    integer, dimension(:), allocatable :: rowIndices, colIndices
    integer :: i, ix, itheta, ipsi, L, ixMin
    integer :: ispecies
    PetscScalar :: signOfPsiDot
    PetscScalar :: sqrtMass

    if (pointAtX0) then
       ixMin = 2
    else
       ixMin = 1
    end if

    if (.not. makeLocalApproximation) then

      allocate(thetaPartOfMirrorTerm(Ntheta))
      allocate(thetaPartMatrix(Ntheta,Ntheta))
      allocate(spatialPartOfStreamingTermDiagonal1(Ntheta,Ntheta))
      allocate(spatialPartOfStreamingTermDiagonal2(Ntheta,Ntheta))
      allocate(spatialPartOfStreamingTermDiagonal3(Ntheta,Ntheta))
      allocate(spatialPartOfStreamingTermOffDiagonal(Ntheta,Ntheta))

      allocate(rowIndices(Ntheta))
      allocate(colIndices(Ntheta))

      do ispecies = 1,Nspecies
         sqrtMass = sqrt(masses(ispecies))
         do ipsi = ipsiMin, ipsiMax
            do itheta=1,Ntheta
               spatialPartOfStreamingTermDiagonal1(itheta,:) = globalTermMultiplier(ipsi) &
                    * sqrtMass * omega*JHat(itheta,ipsi)*IHat(ipsi)*dPhiHatdpsi(ipsi) &
                    / (psiAHatArray(ipsi)*BHat(itheta,ipsi)*BHat(itheta,ipsi)) &
                    * ddthetaToUse(itheta,:)

               spatialPartOfStreamingTermDiagonal2(itheta,:) =  globalTermMultiplier(ipsi) &
                    * sqrtMass/charges(ispecies)*delta*THats(ispecies,ipsi)*JHat(itheta,ipsi) &
                    /(psiAHatArray(ipsi)*BHat(itheta,ipsi)*BHat(itheta,ipsi)) &
                    * IHat(ipsi)*dBHatdpsi(itheta,ipsi)/BHat(itheta,ipsi) &
                    * ddthetaToUse(itheta,:)

               spatialPartOfStreamingTermDiagonal3(itheta,:) =  globalTermMultiplier(ipsi) &
                    * sqrtMass/charges(ispecies)*delta*THats(ispecies,ipsi)*JHat(itheta,ipsi) &
                    /(psiAHatArray(ipsi)*BHat(itheta,ipsi)*BHat(itheta,ipsi)) &
                    * dIHatdpsi(ipsi) &
                    * ddthetaToUse(itheta,:)

               spatialPartOfStreamingTermOffDiagonal(itheta,:) =  globalTermMultiplier(ipsi) &
                   * sqrtMass/charges(ispecies)*delta*THats(ispecies,ipsi)*JHat(itheta,ipsi) &
                    / (BHat(itheta,ipsi)*BHat(itheta,ipsi)*psiAHatArray(ipsi)) &
                    * (IHat(ipsi)/(two*BHat(itheta,ipsi))*dBHatdpsi(itheta,ipsi) - dIHatdpsi(ipsi)) &
                    * ddthetaToUse(itheta,:)
            end do
            do ix=ixMin,Nx
               thetaPartOfMirrorTerm = globalTermMultiplier(ipsi) &
               * sqrt(masses(ispecies))*(omega*dPhiHatdpsi(ipsi)*IHat(ipsi) &
                    + delta*x2(ix)*THats(ispecies,ipsi)/charges(ispecies)*dIHatdpsi(ipsi)) &
                    * JHat(:,ipsi)*dBHatdtheta(:,ipsi) / (two*psiAHatArray(ipsi)*(BHat(:,ipsi) ** 3))

               do L=0,Nxi_for_x(ix)-1
                  rowIndices = [(getIndex(ispecies,ix,L,itheta,ipsi), itheta=1,Ntheta)]

                  ! Term that is diagonal in L:
                  colIndices = rowIndices

                  ! Add streaming term:
                  thetaPartMatrix = spatialPartOfStreamingTermDiagonal1 &
                       + x2(ix) * ((3*L*L+3*L-2)*spatialPartOfStreamingTermDiagonal2 &
                       - (2*L*L+2*L-1)*spatialPartOfStreamingTermDiagonal3) / ((two*L+3)*(two*L-1))

                  ! Add mirror term:
                  ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                  do itheta=1,Ntheta
                     thetaPartMatrix(itheta,itheta) = L*(L+1)/((two*L-1)*(two*L+3))&
                          * thetaPartOfMirrorTerm(itheta)
                  end do

                  thetaPartMatrix = transpose(thetaPartMatrix)

                  ! Put values in matrix
                  if ((ipsi>1 .and. ipsi<Npsi) .or. boundaryScheme == 3) then
                     call MatSetValuesSparse(matrix, Ntheta, rowIndices, &
                          Ntheta, colIndices, thetaPartMatrix, ADD_VALUES, ierr)
                  else
                     do itheta=1,Ntheta
                        signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                             / (psiAHat*charges(ispecies))
                        if ( &
                             (ipsi==1 .and. ((boundaryScheme == 0 .and. signOfPsiDot<-thresh) .or. boundaryScheme == 2)) &
                           .or. (ipsi==Npsi .and. ((boundaryScheme == 0 .and. signOfPsiDot>thresh) .or. boundaryScheme == 1))) then
                           call MatSetValuesSparse(matrix, 1, rowIndices(itheta), &
                                Ntheta, colIndices, thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                        end if
                     end do
                  end if

                  ! End of term that is diagonal in L.

                  if ((L < Nxi_for_x(ix)-2) .and. (whichMatrix==1 .or. preconditioner_xi==0)) then
                     ! Super-super-diagonal in L:
                     !colIndices = rowIndices + Ntheta*2
                     colIndices = [(getIndex(ispecies,ix,L+2,itheta,ipsi), itheta=1,Ntheta)]


                     ! Add streaming term:
                     thetaPartMatrix = x2(ix)*(L+2)*(L+1)/((two*L+5)*(two*L+3)) &
                          * spatialPartOfStreamingTermOffDiagonal

                     ! Add mirror term:
                     ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                     do itheta=1,Ntheta
                        thetaPartMatrix(itheta,itheta) = (L+3)*(L+2)*(L+1)/((two*L+5)*(two*L+3)) &
                             * thetaPartOfMirrorTerm(itheta)
                     end do

                     thetaPartMatrix = transpose(thetaPartMatrix)

                     ! Put values in matrix
                     if ((ipsi>1 .and. ipsi<Npsi) .or. boundaryScheme == 3) then
                        call MatSetValuesSparse(matrix, Ntheta, rowIndices, &
                             Ntheta, colIndices, thetaPartMatrix, ADD_VALUES, ierr)
                     else
                        do itheta=1,Ntheta
                           signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                                / (psiAHat*charges(ispecies))
                           if ( &
                             (ipsi==1 .and. ((boundaryScheme == 0 .and. signOfPsiDot<-thresh) .or. boundaryScheme == 2)) &
                           .or. (ipsi==Npsi .and. ((boundaryScheme == 0 .and. signOfPsiDot>thresh) .or. boundaryScheme == 1))) then
                              call MatSetValuesSparse(matrix, 1, rowIndices(itheta), &
                                   Ntheta, colIndices, thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                           end if
                        end do
                     end if
 !!$                 ! Put values in matrix, noting that Petsc uses a transposed format relative to Fortran
 !!$                 call MatSetValuesSparse(preconditionerMatrix, Ntheta, rowIndices, &
 !!$                      Ntheta, colIndices, transpose(thetaPartMatrix), ADD_VALUES, ierr)
                  end if

                  if ((L>1) .and. (whichMatrix==1 .or. preconditioner_xi==0)) then
                     ! Sub-sub-diagonal in L:
                     !colIndices = rowIndices - Ntheta*2
                     colIndices = [(getIndex(ispecies,ix,L-2,itheta,ipsi), itheta=1,Ntheta)]

                     ! Streaming term
                     thetaPartMatrix = x2(ix)*L*(L-1)/((two*L-3)*(two*L-1)) &
                          * spatialPartOfStreamingTermOffDiagonal

                     ! Mirror term
                     ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                     do itheta=1,Ntheta
                        thetaPartMatrix(itheta,itheta) = -L*(L-1)*(L-2)/((two*L-3)*(two*L-1)) &
                             * thetaPartOfMirrorTerm(itheta)
                     end do

                     thetaPartMatrix = transpose(thetaPartMatrix)

                     ! Put values in matrix
                     if ((ipsi>1 .and. ipsi<Npsi) .or. boundaryScheme == 3) then
                        call MatSetValuesSparse(matrix, Ntheta, rowIndices, &
                             Ntheta, colIndices, thetaPartMatrix, ADD_VALUES, ierr)
                     else
                        do itheta=1,Ntheta
                           signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                                / (psiAHat*charges(ispecies))
                           if ( &
                                (ipsi==1 .and. ((boundaryScheme == 0 .and. signOfPsiDot<-thresh) .or. boundaryScheme == 2)) &
                                .or. (ipsi==Npsi .and. ((boundaryScheme == 0 .and. signOfPsiDot>thresh) .or. boundaryScheme == 1)) &
                                ) then
                              call MatSetValuesSparse(matrix, 1, rowIndices(itheta), &
                                   Ntheta, colIndices, thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                           end if
                        end do
                     end if
 !!$                 ! Put values in matrix, noting that Petsc uses a transposed format relative to Fortran
 !!$                 call MatSetValuesSparse(preconditionerMatrix, Ntheta, rowIndices, &
 !!$                      Ntheta, colIndices, transpose(thetaPartMatrix), ADD_VALUES, ierr)
                  end if
               end do
            end do
         end do
      end do
      deallocate(rowIndices)
      deallocate(colIndices)
      deallocate(spatialPartOfStreamingTermDiagonal1)
      deallocate(spatialPartOfStreamingTermDiagonal2)
      deallocate(spatialPartOfStreamingTermDiagonal3)
      deallocate(spatialPartOfStreamingTermOffDiagonal)
      deallocate(thetaPartOfMirrorTerm)
      deallocate(thetaPartMatrix)
   end if

  end subroutine streamingAndMirrorTerms1

  subroutine collisionlessDdx(whichMatrix)

    ! *********************************************************
    ! Add the collisionless d/dx term:
    ! *********************************************************

    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    integer :: rowIndex, colIndex
    PetscScalar, dimension(:,:), allocatable :: xPartOfXDot
    PetscScalar, dimension(:), allocatable :: diagonalOfXDot
    integer :: ix, itheta, ipsi, L, ell, ix_row, ix_col, ixMin
    integer :: ispecies
    integer :: keepXCoupling
    PetscScalar :: xDotFactor, LFactor
    PetscScalar :: signOfPsiDot

    keepXCoupling = 1
    if (whichMatrix==0 .and. preconditioner_x==1) then
       keepXCoupling = 0
    end if
    ! When keepXCoupling==1, keep the full x coupling.
    ! When keepXCoupling==0, drop everything off-diagonal in x.

    if (pointAtX0) then
       ixMin = 2
    else
       ixMin = 1
    end if

    if (.not. makeLocalApproximation) then

       allocate(xPartOfXDot(Nx,Nx))
       allocate(diagonalOfXDot(Nx))
       do ipsi=ipsiMin,ipsiMax
          do L=0,(Nxi-1)
             if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
                ddxToUse = ddxPreconditioner
             else
                ddxToUse = ddx
             end if
             do ispecies = 1,Nspecies
                do ix=1,Nx
                   xPartOfXDot(ix,:) = x(ix)*(delta*dTHatdpsis(ispecies,ipsi)*x2(ix)/(two*charges(ispecies))&
                        + omega*dPhiHatdpsi(ipsi)) * sqrt(masses(ispecies)) * ddxToUse(ix,:)
                end do
                ! This next line with the transpose was commented out since I switched from MatSetValuesSparse to MatSetValueSparse:
                !xPartOfXDot = transpose(xPartOfXDot)  ! PETSc uses the opposite convention of Fortran
                do itheta=1,Ntheta
                   signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                        / (psiAHat*charges(ispecies))
                   if ((ipsi > 1 .and. ipsi < Npsi) .or. (boundaryScheme == 3) &
                        .or. ((ipsi == 1) .and. (boundaryScheme == 2 .or. (signOfPsiDot<-thresh .and. boundaryScheme == 0))) &
                        .or. ((ipsi==Npsi) .and. (boundaryScheme == 1 .or. (signOfPsiDot>thresh .and. boundaryScheme == 0)))) then
                      ! We're either in the interior, or on a boundary point at which trajectories leave the domain,
                      ! so impose the kinetic equation here.

                      xDotFactor =  globalTermMultiplier(ipsi) &
                           * JHat(itheta,ipsi) * IHat(ipsi) * dBHatdtheta(itheta,ipsi) &
                           / (two*psiAHatArray(ipsi)*(BHat(itheta,ipsi) ** 3))

                      ! Term that is diagonal in L:
                      LFactor = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))*xDotFactor
                      ell = L
                      do ix_col=max(ixMin,min_x_for_L(ell)),Nx
                         colIndex = getIndex(ispecies,ix_col,ell,itheta,ipsi)
                         do ix_row=max(ixMin,min_x_for_L(L)),Nx
                            rowIndex = getIndex(ispecies,ix_row,L,itheta,ipsi)
                            call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                 LFactor*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                         end do
                      end do

                      if (whichMatrix==1 .or. preconditioner_xi==0) then
                         ! Term that is super-super-diagonal in L:
                         if (L<(Nxi-2)) then
                            ell = L+2
                            LFactor = (L+1)*(L+2)/((two*L+5)*(2*L+3))*xDotFactor
                            do ix_col=max(ixMin,min_x_for_L(ell)),Nx
                               colIndex = getIndex(ispecies,ix_col,ell,itheta,ipsi)
                               do ix_row=max(ixMin,min_x_for_L(L)),Nx
                                  rowIndex = getIndex(ispecies,ix_row,L,itheta,ipsi)
                                  call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                       LFactor*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                               end do
                            end do
                         end if

                         ! Term that is sub-sub-diagonal in L:
                         if (L>1) then
                            ell = L-2
                            LFactor = L*(L-1)/((two*L-3)*(2*L-1))*xDotFactor
                            do ix_col=max(ixMin,min_x_for_L(ell)),Nx
                               colIndex = getIndex(ispecies,ix_col,ell,itheta,ipsi)
                               do ix_row=max(ixMin,min_x_for_L(L)),Nx                            
                                  rowIndex = getIndex(ispecies,ix_row,L,itheta,ipsi)
                                  call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                       LFactor*xPartOfXDot(ix_row,ix_col), ADD_VALUES, ierr)
                               end do
                            end do
                         end if
                      end if

                   end if
                end do
             end do
          end do
       end do
       deallocate(xPartOfXDot)
       deallocate(diagonalOfXDot)
    end if

  end subroutine collisionlessDdx

  subroutine collisionlessDdpsi(whichMatrix,upwinding)

    ! *********************************************************
    ! Add the collisionless d/dpsi term:
    ! *********************************************************

    integer, intent(in) :: whichMatrix
    logical, intent(in) :: upwinding
    PetscErrorCode :: ierr
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar, dimension(:), allocatable :: thetaPartOfPsiDot
    PetscScalar, dimension(:,:), allocatable :: ddpsiToUse
    PetscScalar, dimension(:,:), allocatable :: localddpsiToUse, everythingButLInPsiDot
    integer :: ipsiMinInterior, ipsiMaxInterior, ell
    integer :: i, ix, itheta, ipsi, L, ixMin
    integer :: ispecies
    logical :: includeddpsiTermThisTime
    integer :: ipsiMinForThisTheta, ipsiMaxForThisTheta, NpsiToUse
    PetscScalar :: LFactor

    if (pointAtX0) then
       ixMin = 2
    else
       ixMin = 1
    end if

    includeddpsiTermThisTime = includeddpsiTerm
    if (makeLocalApproximation) then
       includeddpsiTermThisTime = .false.
    end if
    if (whichMatrix==0) then
       if (preconditioner_psi==2 .or. preconditioner_psi == 4) then
          includeddpsiTermThisTime = .false.
       end if
    end if
    if (includeddpsiTermThisTime) then
       allocate(colIndices(Npsi))
       allocate(thetaPartOfPsiDot(Npsi))
       allocate(ddpsiToUse(Npsi,Npsi))

       if (whichMatrix==1) then
          ddpsiToUse = ddpsiLeft
       else
          select case (preconditioner_psi)
          case (0)
             ddpsiToUse = ddpsiLeft
          case (1)
             ddpsiToUse = ddpsiForPreconditioner
          case (2,4)
             print *,"Error! Program should not get here."
             stop
          case (3)
             ! Keep only the diagonal of ddpsi, which is typically nonzero only near the boundaries.
             ddpsiToUse = 0
             do i=1,Npsi
                ddpsiToUse(i,i) = ddpsiLeft(i,i)
             end do
          case default
             print *,"Error! Invalid setting for preconditioner_psi."
             stop
          end select
       end if
       
       !allocate(localddpsiToUse(localNpsiInterior,Npsi))
       do ispecies = 1,Nspecies
          do itheta = 1, Ntheta
             thetaPartOfPsiDot = -oneHalf * globalTermMultiplier(:) &
             * sqrt(masses(ispecies)) * delta * JHat(itheta,:) &
                  * IHat(:) * THats(ispecies,:) * dBHatdtheta(itheta,:) &
                  / (charges(ispecies) * psiAHatArray(:) * (BHat(itheta,:) ** 3))
             if (upwinding .and. (maxval(thetaPartOfPsiDot)>0) .and. (minval(thetaPartOfPsiDot)<0)) then
                print *,"Warning: psiDot at itheta =",itheta,&
                     " changes sign with psi, so upwinding is not well-defined."
             end if
             ipsiMinForThisTheta = ipsiMin
             ipsiMaxForThisTheta = ipsiMax
             if ((procThatHandlesLeftBoundary .and. (thetaPartOfPsiDot(1) .ge. 0)) &
                  .and. boundaryScheme /= 3 .and. boundaryScheme /= 2) then
                ipsiMinForThisTheta = 2
             end if
             if ((procThatHandlesRightBoundary .and. (thetaPartOfPsiDot(Npsi) .le. 0)) &
                  .and. boundaryScheme /= 3 .and. boundaryScheme /= 1) then
                ipsiMaxForThisTheta = Npsi-1
             end if

             NpsiToUse = ipsiMaxForThisTheta - ipsiMinForThisTheta + 1
             allocate(localddpsiToUse(NpsiToUse,Npsi))
             allocate(rowIndices(NpsiToUse))
             allocate(everythingButLInPsiDot(Npsi,NpsiToUse))
             ! The line below must be modified if we want to use upwinded derivatives:
             localddpsiToUse = ddpsiToUse(ipsiMinForThisTheta:ipsiMaxForThisTheta,:)

             !signOfPsiDot = thetaPartOfPsiDot(1)
             !if (signOfPsiDot > 0) then
             !   localddpsiToUse = localddpsiLeftInterior
             !else
             !   localddpsiToUse = localddpsiRightInterior
             !end if
             do ix = ixMin, Nx
                do ipsi=1,NpsiToUse
                   ! Build everythingButLInPsiDot as the transpose of what it would normally be
                   ! because of the difference between Fortran and Petsc convention.
                   everythingButLInPsiDot(:,ipsi) = x2(ix) * thetaPartOfPsiDot(ipsi-1+ipsiMinForThisTheta) &
                        * localddpsiToUse(ipsi,:)
                end do
                do L=0,(Nxi_for_x(ix)-1)
                   !rowIndices = [(ipsi-1, ipsi=ipsiMinInterior,ipsiMaxInterior)]*localMatrixSize &
                   rowIndices = [(getIndex(ispecies,ix,L,itheta,ipsi), ipsi=ipsiMinForThisTheta,ipsiMaxForThisTheta)]

                   ! Term that is diagonal in L:
                   ell = L
                   colIndices = [(getIndex(ispecies,ix,ell,itheta,ipsi), ipsi=1,Npsi)]
                   LFactor = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))
                   call MatSetValuesSparse(matrix, NpsiToUse, rowIndices, Npsi, colIndices, &
                        LFactor*everythingButLInPsiDot, ADD_VALUES, ierr)

                   if (whichMatrix==1 .or. preconditioner_xi==0) then
                      ! Term that is super-super-diagonal in L:
                      if (L<(Nxi_for_x(ix)-2)) then
                         ell = L+2
                         colIndices = [(getIndex(ispecies,ix,ell,itheta,ipsi), ipsi=1,Npsi)]
                         LFactor = (L+1)*(L+2)/((two*L+5)*(2*L+3))
                         call MatSetValuesSparse(matrix, NpsiToUse, rowIndices, Npsi, colIndices, &
                              LFactor*everythingButLInPsiDot, ADD_VALUES, ierr)
                      end if

                      ! Term that is sub-sub-diagonal in L:
                      if (L>1) then
                         ell = L-2
                         colIndices = [(getIndex(ispecies,ix,ell,itheta,ipsi), ipsi=1,Npsi)]
                         LFactor = L*(L-1)/((two*L-3)*(2*L-1))
                         call MatSetValuesSparse(matrix, NpsiToUse, rowIndices, Npsi, colIndices, &
                              LFactor*everythingButLInPsiDot, ADD_VALUES, ierr)
                      end if
                   end if

                end do
             end do
             deallocate(localddpsiToUse)
             deallocate(rowIndices)
             deallocate(everythingButLInPsiDot)
          end do
       end do
       deallocate(colIndices)
       deallocate(thetaPartOfPsiDot)
       deallocate(ddpsiToUse)
    end if

  end subroutine collisionlessDdpsi

  subroutine collisionOperator(whichMatrix)

    ! *********************************************************
    ! Add the collision operator
    ! *********************************************************

    integer, intent(in) :: whichMatrix
    PetscErrorCode :: ierr
    integer, dimension(:), allocatable :: rowIndices, colIndices
    PetscScalar, dimension(:), allocatable :: erfs, xb, expxb2, Psi_Chandra
    PetscScalar, dimension(:,:), allocatable :: tempMatrix, tempMatrix2
    PetscScalar, dimension(:,:), allocatable :: M11, M21, M32, LaplacianTimesX2WithoutL
    PetscScalar, dimension(:,:), allocatable :: M22BackslashM21, M33BackslashM32
    PetscScalar, dimension(:,:,:), allocatable :: M22BackslashM21s, M33BackslashM32s
    PetscScalar, dimension(:,:), allocatable :: KWithoutThetaPart, M22, M33, M12, M13
    PetscScalar, dimension(:,:), allocatable :: nuDHat
    PetscScalar, dimension(:,:), allocatable :: fToFInterpolationMatrix
    PetscScalar, dimension(:,:), allocatable :: potentialsToFInterpolationMatrix
    PetscScalar, dimension(:,:,:,:), allocatable :: CECD
    PetscScalar :: temp, temp1, temp2, speciesFactor, speciesFactor2
    PetscScalar :: signOfPsiDot
    PetscScalar :: T32, sqrt_m
    integer :: i, j, ix, ix_row, ix_col, itheta, ipsi, L, ixMin, ixMinCol
    integer :: rowIndex, colIndex
    integer :: iSpeciesA, iSpeciesB
    integer :: scheme
    integer :: LAPACKInfo
    integer, dimension(:), allocatable :: IPIV  ! Needed by LAPACK

    allocate(xb(Nx))
    allocate(expxb2(Nx))
    allocate(erfs(Nx))
    allocate(Psi_Chandra(Nx))
    allocate(nuDHat(Nspecies, Nx))

    if (xDerivativeScheme == 0 .or. xDerivativeScheme == 1) then
      allocate(M21(NxPotentials, Nx))
      allocate(M32(NxPotentials, NxPotentials))
      allocate(M22BackslashM21(NxPotentials, Nx))
      allocate(M33BackslashM32(NxPotentials, NxPotentials))
      allocate(M22BackslashM21s(NL,NxPotentials, Nx))
      allocate(M33BackslashM32s(NL,NxPotentials, NxPotentials))
      allocate(LaplacianTimesX2WithoutL(NxPotentials, NxPotentials))

      allocate(M12(Nx,NxPotentials))
      allocate(M13(Nx,NxPotentials))
      allocate(M22(NxPotentials,NxPotentials))
      allocate(M33(NxPotentials,NxPotentials))

      allocate(potentialsToFInterpolationMatrix(Nx, NxPotentials))
    endif

    allocate(M11(Nx,Nx))

    allocate(KWithoutThetaPart(Nx,Nx))
    allocate(fToFInterpolationMatrix(Nx,Nx))
    allocate(CECD(Nspecies, Nspecies, Nx, Nx))


    allocate(IPIV(NxPotentials))

    ! *********************************************************
    ! In preparation for adding the collision operator,
    ! create several matrices which will be needed.
    ! *********************************************************

    allocate(rowIndices(Nx))
    allocate(colIndices(Nx))
    if (xDerivativeScheme == 0 .or. xDerivativeScheme == 1) then
      allocate(tempMatrix(Nx, NxPotentials))
      allocate(tempMatrix2(NxPotentials, NxPotentials))
    endif

    if (pointAtX0) then
       ixMin = 2
    else
       ixMin = 1
    end if

    if (whichMatrix==0 .and. preconditioner_x==5) then
       ! For other preconditioner_x values, the simplified x derivative
       ! is implemented later instead of here.
       ddxToUse = ddxPreconditioner
       d2dx2ToUse = d2dx2Preconditioner
       ! Since I set ddxToUse here, the same ddxToUse is used for all values of L.
       ! This means that preconditioner_x_min_L does not behave as you might expect when
       ! preconditioner_x=5.
    else
       ddxToUse = ddx
       d2dx2ToUse = d2dx2
    end if

    if (xDerivativeScheme == 0 .or. xDerivativeScheme == 1) then
      ! Matrices that are only needed if using the old Rosenbluth potential scheme

      ! First assemble rows 2 and 3 of the block linear system, since they
      ! are independent of psi and independent of species.

      M32 = zero
      M21 = 4*pi*regridPolynomialToUniform
      do i=2,NxPotentials-1
         M21(i,:) = M21(i,:)*xPotentials(i)*xPotentials(i)
         M32(i,i) = -2*xPotentials(i)*xPotentials(i)
      end do
      M21(1,:)=zero
      M21(NxPotentials,:)=zero
      M32(1,:)=zero
      M32(NxPotentials,:)=zero
      do i=1,NxPotentials
         LaplacianTimesX2WithoutL(i,:) = xPotentials(i)*xPotentials(i)*d2dx2Potentials(i,:) &
              + 2 * xPotentials(i) * ddxPotentials(i,:)
      end do

      do L=0,(NL-1)
         M22 = LaplacianTimesX2WithoutL
         do i=1,NxPotentials
            M22(i,i) = M22(i,i) - L*(L+1)
         end do

         ! Add Dirichlet or Neumann boundary condition for potentials at x=0:
         if (L==0) then
            M22(1,:)=ddxPotentials(1,:)
         else
            M22(1,:) = 0
            M22(1,1) = 1
         end if
         M33 = M22;

         ! Add Robin boundary condition for potentials at x=xMax:
         M22(NxPotentials,:) = xMaxNotTooSmall*ddxPotentials(NxPotentials,:)
         M22(NxPotentials,NxPotentials) = M22(NxPotentials,NxPotentials) + L+1

         ! Boundary condition for G:
         M33(NxPotentials,:) = xMaxNotTooSmall*xMaxNotTooSmall*d2dx2Potentials(NxPotentials,:) &
              + (2*L+1)*xMaxNotTooSmall*ddxPotentials(NxPotentials,:)
         M33(NxPotentials,NxPotentials) = M33(NxPotentials,NxPotentials) + (L*L-1)

         if (L /= 0) then
            M22(NxPotentials,1)=0
            M33(NxPotentials,1)=0
         end if

         ! Call LAPACK subroutine DGESV to solve a linear system
         ! Note: this subroutine changes M22 and M33!
         M22BackslashM21 = M21  ! This will be overwritten by LAPACK.
#if defined(PETSC_USE_REAL_SINGLE)
         call SGESV(NxPotentials, Nx, M22, NxPotentials, IPIV, M22BackslashM21, NxPotentials, LAPACKInfo)
#else
         call DGESV(NxPotentials, Nx, M22, NxPotentials, IPIV, M22BackslashM21, NxPotentials, LAPACKInfo)
#endif
         if (LAPACKInfo /= 0) then
            print *, "Error in LAPACK call: info = ", LAPACKInfo
            stop
         end if
         M33BackslashM32 = M32  ! This will be overwritten by LAPACK.
#if defined(PETSC_USE_REAL_SINGLE)
         call SGESV(NxPotentials, NxPotentials, M33, NxPotentials, IPIV, M33BackslashM32, NxPotentials, LAPACKInfo)
#else
         call DGESV(NxPotentials, NxPotentials, M33, NxPotentials, IPIV, M33BackslashM32, NxPotentials, LAPACKInfo)
#endif
         if (LAPACKInfo /= 0) then
            print *, "Error in LAPACK call: info = ", LAPACKInfo
            stop
         end if

         M33BackslashM32s(L+1,:,:) = M33BackslashM32
         M22BackslashM21s(L+1,:,:) = M22BackslashM21
      end do
    end if

    do ipsi = ipsiMin, ipsiMax

       nuDHat = zero
       CECD = zero
       ! Before adding the collision operator, we must loop over both species
       ! to build several terms in the operator.
       ! row is species a, column is species b
       do iSpeciesA = 1,Nspecies
          do iSpeciesB = 1,Nspecies
             speciesFactor = sqrt(THats(iSpeciesA,ipsi)*masses(iSpeciesB) &
                  / (THats(iSpeciesB,ipsi) * masses(iSpeciesA)))
             xb =  x * speciesFactor
             expxb2 = exp(-xb*xb)
             do ix=1,Nx
                ! erf is vectorized in gfortran but not pathscale
                temp1 = xb(ix)
#ifdef USE_GSL_ERF
                call erf(temp1, temp2)
#else
                temp2 = erf(temp1)
#endif
                erfs(ix) = temp2
             end do
             Psi_Chandra = (erfs - 2/sqrtpi * xb * expxb2) / (2*xb*xb)

             T32 = THats(iSpeciesA,ipsi) * sqrt(THats(iSpeciesA,ipsi))

             ! Build the pitch-angle scattering frequency:
             nuDHat(iSpeciesA, :) =  nuDHat(iSpeciesA, :) &
                  + (three*sqrtpi/four) / T32 &
                  * charges(iSpeciesA)*charges(iSpeciesA)*charges(iSpeciesB)*charges(iSpeciesB) &
                  * nHats(iSpeciesB,ipsi)*(erfs - Psi_Chandra)/(x*x*x)

             ! Given a vector of function values on the species-B grid, multiply the vector
             ! by this regridding matrix to obtain its values on the species-A grid:
             if (iSpeciesA /= iSpeciesB) then
                call polynomialInterpolationMatrix(Nx, Nx, x, xb, expx2, &
                     expxb2, fToFInterpolationMatrix)
             else
                fToFInterpolationMatrix = zero
                do ix=1,Nx
                   fToFInterpolationMatrix(ix, ix) = one
                end do
             end if

             ! Using the resulting interpolation matrix,
             ! add CD (the part of the field term independent of Rosenbluth potentials.
             ! CD is dense in the species indices.

             speciesFactor = 3 * nHats(iSpeciesA,ipsi)  * masses(iSpeciesA)/masses(iSpeciesB) &
                  * charges(iSpeciesA)*charges(iSpeciesA)*charges(iSpeciesB)*charges(iSpeciesB) / T32

             do ix=1,Nx
                CECD(iSpeciesA, iSpeciesB, ix, :) = CECD(iSpeciesA, iSpeciesB, ix, :) &
                     + speciesFactor * expx2(ix) * fToFInterpolationMatrix(ix, :)
             end do

             ! Done adding CD. Now add energy scattering (CE).
             ! Unlike CD, CE is diagonal in the species index.

             speciesFactor = 3*sqrtpi/four * nHats(iSpeciesB,ipsi)  &
                  * charges(iSpeciesA)*charges(iSpeciesA)*charges(iSpeciesB)*charges(iSpeciesB) / T32

             do ix=1,Nx
                !Now add the d2dx2 and ddx terms in CE:
                !CE is diagonal in the species indices, so use iSpeciesA for both indices in CECD:
                CECD(iSpeciesA, iSpeciesA, ix, :) = CECD(iSpeciesA, iSpeciesA, ix, :) &
                     + speciesFactor * (Psi_Chandra(ix)/x(ix)*d2dx2ToUse(ix,:) &
                     + (-2*THats(iSpeciesA,ipsi)*masses(iSpeciesB)/(THats(iSpeciesB,ipsi)*masses(iSpeciesA)) &
                     * Psi_Chandra(ix)*(1-masses(iSpeciesA)/masses(iSpeciesB)) &
                     + (erfs(ix)-Psi_Chandra(ix))/x2(ix)) * ddxToUse(ix,:))

                ! Lastly, add the part of CE for which f is not differentiated:
                ! CE is diagonal in the species indices, so use iSpeciesA for both indices in CECD:
                CECD(iSpeciesA, iSpeciesA, ix, ix) = CECD(iSpeciesA, iSpeciesA, ix, ix) &
                     + speciesFactor *4/sqrtpi*THats(iSpeciesA,ipsi)/THats(iSpeciesB,ipsi) &
                     *sqrt(THats(iSpeciesA,ipsi)*masses(iSpeciesB)/(THats(iSpeciesB,ipsi)*masses(iSpeciesA))) &
                     * expxb2(ix)

             end do

          end do
       end do


       ! *****************************************************************
       ! Now we are ready to add the collision operator to the main matrix.
       ! *****************************************************************

       do L=0, Nxi-1
          if (L>0 .and. pointAtX0) then
             ixMinCol = 2
          else
             ixMinCol = 1
          end if

          !print *,"Adding L=",L
          do iSpeciesA = 1,Nspecies
             sqrt_m = sqrt(masses(iSpeciesA))
             do iSpeciesB = 1,Nspecies
                if (includeCollisionOperator .and. ((iSpeciesA == iSpeciesB .or. preconditioner_species==0) &
                     .or. whichMatrix==1)) then
                   ! Build M11
                   M11 = -nu_r * CECD(iSpeciesA, iSpeciesB,:,:)
                   if (iSpeciesA == iSpeciesB) then
                      do ix=1,Nx
                         M11(ix,ix) = M11(ix,ix) - nu_r * (-oneHalf*nuDHat(iSpeciesA,ix)*L*(L+1))
                      end do
                   end if

                   if (L < NL) then
                   !   if (.false.) then
                      ! Add Rosenbluth potential terms.

                      if (xDerivativeScheme==2 .or. xDerivativeScheme==3) then
                         ! New scheme for the Rosenbluth potential terms.
                         do ix=1,Nx
                            ! The DKE normalization in perfect has an extra sqrt(m) compared to the normalization in SFINCS. Add the factor here:
                            M11(ix, :) = M11(ix,:) &
                                 - nu_r * sqrt_m * RosenbluthPotentialTerms(iSpeciesA,iSpeciesB,L+1,ix,:,ipsi-ipsiMin+1) 
                         end do

                         KWithoutThetaPart = M11

                      else
                         ! Old scheme for the Rosenbluth potential terms.

                         speciesFactor2 = sqrt(THats(iSpeciesA,ipsi)*masses(iSpeciesB) &
                              / (THats(iSpeciesB,ipsi) * masses(iSpeciesA)))

                         ! Build M13:
                         scheme = 2
                         call interpolationMatrix(NxPotentials, Nx, xPotentials, x*speciesFactor2, &
                              potentialsToFInterpolationMatrix, scheme, L)
                         
                         speciesFactor = -nu_r * 3/(2*pi)*nHats(iSpeciesA,ipsi) &
                              * charges(iSpeciesA)*charges(iSpeciesA)*charges(iSpeciesB)*charges(iSpeciesB) &
                              / (THats(iSpeciesA,ipsi) * sqrt(THats(iSpeciesA,ipsi))) &
                              * THats(iSpeciesB,ipsi)*masses(iSpeciesA)/(THats(iSpeciesA,ipsi)*masses(iSpeciesB))
                         
                         tempMatrix = matmul(potentialsToFInterpolationMatrix, d2dx2Potentials)
                         do ix=1,Nx
                            M13(ix, :) = speciesFactor*expx2(ix)*x2(ix)*tempMatrix(ix,:)
                         end do
                         
                         ! Build M12:
                         scheme = 1
                         call interpolationMatrix(NxPotentials, Nx, xPotentials, x*speciesFactor2, &
                              potentialsToFInterpolationMatrix, scheme, L)
                         
                         temp = 1-masses(iSpeciesA)/masses(iSpeciesB)
                         do i=1,NxPotentials
                            tempMatrix2(i,:) = temp*xPotentials(i)*ddxPotentials(i,:)
                            tempMatrix2(i,i) = tempMatrix2(i,i) + one
                         end do
                         tempMatrix = matmul(potentialsToFInterpolationMatrix, tempMatrix2)
                         do ix=1,Nx
                            M12(ix,:) = -speciesFactor*expx2(ix)*tempMatrix(ix,:)
                         end do
                         
                         ! Possibly add Dirichlet boundary condition for potentials at x=0:
                         if (L /= 0) then
                            M12(:,1) = 0
                            M13(:,1) = 0
                         end if
                         
                         !KWithoutThetaPart = M11 -  (M12 - M13 * (M33 \ M32)) * (M22 \ M21);
                         KWithoutThetaPart = M11 - matmul(M12 - matmul(M13, M33BackslashM32s(L+1,:,:)),&
                              M22BackslashM21s(L+1,:,:))
                      end if
                   else
                      KWithoutThetaPart = M11;
                   end if

                   if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
                      ! We're making the preconditioner, so simplify the x part of the matrix if desired.
                      select case (preconditioner_x)
                      case (0)
                         ! Do nothing.
                      case (1)
                         ! Keep only diagonal in x:
                         do i=1,Nx
                            do j=1,Nx
                               if (i /= j) then
                                  KWithoutThetaPart(i,j) = zero
                               end if
                            end do
                         end do
                      case (2)
                         ! Keep only upper-triangular part:
                         do i=2,Nx
                            do j=1,(i-1)
                               KWithoutThetaPart(i,j) = zero
                            end do
                         end do
                      case (3,5)
                         ! Keep only tridiagonal part:
                         do i=1,Nx
                            do j=1,Nx
                               if (abs(i-j)>1) then
                                  KWithoutThetaPart(i,j) = zero
                               end if
                            end do
                         end do
                      case (4)
                         ! Keep only the diagonal and super-diagonal:
                         do i=1,Nx
                            do j=1,Nx
                               if (i /= j .and. j /= (i+1)) then
                                  KWithoutThetaPart(i,j) = zero
                               end if
                            end do
                         end do
                      case default
                         print *,"Error! Invalid preconditioner_x"
                         stop
                      end select

                   end if

                   !! Don't need any more, now using MatSetValueSparse below
                   !! ! PETSc and Fortran use row-major vs column-major:
                   !! KWithoutThetaPart = transpose(KWithoutThetaPart)

                   do itheta=1,Ntheta
                      signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                           / (psiAHat*charges(iSpeciesA))
                      if ((ipsi > 1 .and. ipsi < Npsi) .or. (boundaryScheme == 3) &
                           .or. (ipsi == 1 .and. ((signOfPsiDot < -thresh .and. boundaryScheme == 0) &
                           .or. boundaryScheme == 2 .or. leftBoundaryScheme == 2)) &
                           .or. (ipsi==Npsi .and. ((signOfPsiDot > thresh .and. boundaryScheme == 0) &
                           .or. boundaryScheme == 1 .or. rightBoundaryScheme == 2))) then
                         ! We're either in the interior, or on a boundary point at which trajectories leave the domain,
                         ! so impose the kinetic equation here.

                         do ix_row=max(ixMin,min_x_for_L(L)),Nx
                           rowIndex = getIndex(iSpeciesA,ix_row,L,itheta,ipsi)
                           do ix_col=max(ixMinCol,min_x_for_L(L)),Nx
                             colIndex = getIndex(iSpeciesB,ix_col,L,itheta,ipsi)
                             call MatSetValueSparse(matrix, rowIndex, colIndex, &
                                  KWithoutThetaPart(ix_row,ix_col), ADD_VALUES, ierr)
                           end do
                         end do
                      end if
                   end do

                   if (procThatHandlesLeftBoundary .and. ipsi==1 .and. boundaryScheme /= 2  .and. boundaryScheme /= 3) then
                      do itheta=1,Ntheta
                         do ix_row=min_x_for_L(L),Nx
                           rowIndex = getIndex(iSpeciesA,ix_row,L,itheta,1)
                           do ix_col=min_x_for_L(L),Nx
                             colIndex = getIndex(iSpeciesB,ix_col,L,itheta,1)
                             call MatSetValueSparse(leftMatrix, rowIndex, colIndex, &
                                  KWithoutThetaPart(ix_row,ix_col), ADD_VALUES, ierr)
                           end do
                         end do
                      end do
                   end if

                   if (procThatHandlesRightBoundary .and. ipsi==Npsi .and. boundaryScheme /= 1  .and. boundaryScheme /= 3) then
                      do itheta=1,Ntheta
                         do ix_row=min_x_for_L(L),Nx
                           rowIndex = getIndex(iSpeciesA,ix_row,L,itheta,1)
                           do ix_col=min_x_for_L(L),Nx
                             colIndex = getIndex(iSpeciesB,ix_col,L,itheta,1)
                             call MatSetValueSparse(rightMatrix, rowIndex, colIndex, &
                                  KWithoutThetaPart(ix_row,ix_col), ADD_VALUES, ierr)
                           end do
                         end do
                      end do
                   end if

                end if
             end do
          end do
       end do
    end do

    deallocate(rowIndices)
    deallocate(colIndices)
    if (xDerivativeScheme == 0 .or. xDerivativeScheme == 1) then
      deallocate(tempMatrix)
      deallocate(tempMatrix2)
    endif

    deallocate(xb)
    deallocate(expxb2)
    deallocate(erfs)
    deallocate(Psi_Chandra)
    deallocate(nuDHat)

    if (xDerivativeScheme == 0 .or. xDerivativeScheme == 1) then
      deallocate(M21)
      deallocate(M32)
      deallocate(M22BackslashM21)
      deallocate(M33BackslashM32)
      deallocate(M22BackslashM21s)
      deallocate(M33BackslashM32s)
      deallocate(LaplacianTimesX2WithoutL)

      deallocate(M12)
      deallocate(M13)
      deallocate(M22)
      deallocate(M33)

      deallocate(potentialsToFInterpolationMatrix)
    endif

    deallocate(M11)

    deallocate(KWithoutThetaPart)
    deallocate(fToFInterpolationMatrix)
    deallocate(CECD)

    deallocate(IPIV)

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Done adding collision operator.
    !
    ! *******************************************************************************
    ! *******************************************************************************

  end subroutine collisionOperator

  subroutine radialBoundaryConditionDiagonal()

    ! *******************************************************************************
    ! Put a 1 on the matrix diagonal where appropriate to enforce the radial boundary condition
    ! *******************************************************************************

    PetscErrorCode :: ierr
    integer :: ix, itheta, ipsi, L, index
    integer :: ispecies
    PetscScalar :: signOfPsiDot

    if (procThatHandlesLeftBoundary .and. leftBoundaryScheme /= 2 .and. boundaryScheme /= 3 .and. boundaryScheme /= 2) then
       ipsi=1
       do ispecies = 1,Nspecies
          do itheta=1,Ntheta
             signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                  /(psiAHat*charges(ispecies))
             if ((signOfPsiDot > -thresh .and. boundaryScheme == 0) .or. boundaryScheme == 1) then
                do ix=1,Nx
                   do L=0,(Nxi_for_x(ix)-1)
                      index = getIndex(ispecies,ix,L,itheta,1)
                      call MatSetValueSparse(matrix, index, index, one, ADD_VALUES, ierr)
                   end do
                end do
             end if
          end do
       end do
    end if
    if (procThatHandlesRightBoundary .and. rightBoundaryScheme /= 2 .and. boundaryScheme /= 3 .and. boundaryScheme /= 1) then
       ipsi = Npsi
       do ispecies = 1,Nspecies
          do itheta=1,Ntheta
             signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                  / (psiAHat*charges(ispecies))
             if ((signOfPsiDot < thresh .and. boundaryScheme == 0) .or. boundaryScheme == 2) then
                do ix=1,Nx
                   do L=0,(Nxi_for_x(ix)-1)
                      index = getIndex(ispecies,ix,L,itheta,ipsi)
                      call MatSetValueSparse(matrix, index, index, one, ADD_VALUES, ierr)
                   end do
                end do
             end if
          end do
       end do

    end if

  end subroutine radialBoundaryConditionDiagonal

  subroutine xBoundaryCondition(whichMatrix)

    ! *******************************************************************************
    ! Enforce boundary condition at x=0 if needed (if there is a point there)
    ! *******************************************************************************

    PetscErrorCode :: ierr
    integer :: ix, L, itheta, ipsi, index, ix_row, ix_col, rowIndex, colIndex, ispecies
    integer, intent(in) :: whichMatrix

    if (pointAtX0) then
       ! For L > 0 modes, impose f=0 at x=0:
       ix = 1
       do ipsi = ipsiMin,ipsiMax
          do L = 1,(Nxi_for_x(ix)-1)
             do itheta = 1,Ntheta
                do ispecies = 1,Nspecies
                   index = getIndex(ispecies,ix,L,itheta,ipsi)
                   call MatSetValue(matrix, index, index, one, ADD_VALUES, ierr)
                end do
             end do
          end do
       end do

       ! For L=0 mode, impose regularity (df/dx=0) at x=0:
       L=0
       if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
          ddxToUse = ddxPreconditioner
       else
          ddxToUse = ddx
       end if
       ix_row = 1
       do ipsi = ipsiMin,ipsiMax
          do itheta = 1,Ntheta
             do ispecies = 1,Nspecies
                rowIndex = getIndex(ispecies,ix_row,L,itheta,ipsi)
                do ix_col = 1,Nx
                   colIndex = getIndex(ispecies,ix_col,L,itheta,ipsi)
                   call MatSetValueSparse(matrix, rowIndex, colIndex, ddxToUse(ix_row,ix_col), ADD_VALUES, ierr)
                end do
             end do
          end do
       end do
    end if

  end subroutine xBoundaryCondition

  subroutine sources()

    ! *******************************************************************************
    ! Add sources:
    ! *******************************************************************************

    PetscErrorCode :: ierr
    PetscScalar, dimension(:), allocatable :: sourceXPart
  
    PetscScalar :: signOfPsiDot, xPartOfSource
    integer :: i, ix, itheta, ipsi, L
    integer :: ispecies, isources, iextraSources, ispeciesIndepSources
    integer :: rowIndex, colIndex
    integer :: this_ipsiMin,this_ipsiMax

    if (ipsiMax < lowestEnforcedIpsi) then
       ! this processor do not own any indices with sources
       ! UNTESTED
       return
    end if
    
    if (ipsiMin > highestEnforcedIpsi) then
       ! this processor do not own any indices with sources
       ! UNTESTED
       return 
    end if

    this_ipsiMin = max(ipsiMin,lowestEnforcedIpsi)
    this_ipsiMax = min(ipsiMax,highestEnforcedIpsi)

    allocate(sourceXPart(Nx))
    
    do isources = 1,Nsources
       !call initializeSourcesThetaPart(isources,sourceThetaPart(isource,:),sourceThetaPartFSA)

       !call initializeSourcesVPart(isources,sourceXPart,L)
       call initializeSources(isources,sourceXPart,L,sourceThetaPart(isources,:),sourceThetaPartFSA(isources,:))
          
       do ix=1,Nx             
          do ipsi=this_ipsiMin,this_ipsiMax
             do ispecies = 1,Nspecies
                do itheta=1,Ntheta
                   signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                        / (psiAHat*charges(ispecies))
                   if ((ipsi > 1 .and. ipsi < Npsi) .or. (boundaryScheme == 3) &
                           .or. (ipsi == 1 .and. ((signOfPsiDot < -thresh .and. boundaryScheme == 0) &
                           .or. boundaryScheme == 2 .or. leftBoundaryScheme == 2)) &
                           .or. (ipsi==Npsi .and. ((signOfPsiDot > thresh .and. boundaryScheme == 0) &
                           .or. boundaryScheme == 1 .or. rightBoundaryScheme == 2))) then
                      ! We're either in the interior, or on a boundary point at which trajectories leave the domain,
                      ! so impose the kinetic equation here.
                      rowIndex = getIndex(ispecies,ix,L,itheta,ipsi)
                      colIndex = getIndexSources(isources,ispecies,ipsi)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           sourceThetaPart(isources,itheta) * sourceXPart(ix), ADD_VALUES, ierr)
                   end if
                end do
             end do
          end do
       end do
    end do

    do ispeciesIndepSources = 1,NspeciesIndepSources
       call initializeSpeciesIndepSources(ispeciesIndepSources,sourceXPart,L,&
            speciesIndepSourceThetaPart(ispeciesIndepSources,:),speciesIndepSourceThetaPartFSA(ispeciesIndepSources,:),&
            speciesIndepSourceSpeciesPart(ispeciesIndepSources,:))
       do ix=1,Nx         
          do ipsi=this_ipsiMin,this_ipsiMax
             do ispecies = 1,Nspecies
                do itheta=1,Ntheta
                   signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                        / (psiAHat*charges(ispecies))
                   if ((ipsi > 1 .and. ipsi < Npsi) .or. (boundaryScheme == 3) &
                        .or. (ipsi == 1 .and. ((signOfPsiDot < -thresh .and. boundaryScheme == 0) &
                        .or. boundaryScheme == 2 .or. leftBoundaryScheme == 2)) &
                        .or. (ipsi==Npsi .and. ((signOfPsiDot > thresh .and. boundaryScheme == 0) &
                        .or. boundaryScheme == 1 .or. rightBoundaryScheme == 2))) then
                      ! We're either in the interior, or on a boundary point at which trajectories leave the domain,
                      ! so impose the kinetic equation here.
                      rowIndex = getIndex(ispecies,ix,L,itheta,ipsi)
                      colIndex = getIndexSpeciesIndepSources(ispeciesIndepSources,ipsi)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           speciesIndepSourceSpeciesPart(ispeciesIndepSources,ispecies) &
                           * speciesIndepSourceThetaPart(ispeciesIndepSources,itheta) &
                           * sourceXPart(ix), ADD_VALUES, ierr)
                   end if
                end do
             end do
          end do
       end do
    end do

    do iextraSources = 1,NextraSources
       call initializeExtraSources(iextraSources,sourceXPart,L,&
            extraSourceThetaPart(iextraSources,:),extraSourceThetaPartFSA(iextraSources,:),extraSourceSpeciesPart(iextraSources,:))
          
       do ix=1,Nx         
          do ipsi=this_ipsiMin,this_ipsiMax
             do ispecies = 1,Nspecies
                do itheta=1,Ntheta
                   signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                        / (psiAHat*charges(ispecies))
                   if ((ipsi > 1 .and. ipsi < Npsi) .or. (boundaryScheme == 3) &
                        .or. (ipsi == 1 .and. ((signOfPsiDot < -thresh .and. boundaryScheme == 0) &
                        .or. boundaryScheme == 2 .or. leftBoundaryScheme == 2)) &
                        .or. (ipsi==Npsi .and. ((signOfPsiDot > thresh .and. boundaryScheme == 0) &
                        .or. boundaryScheme == 1 .or. rightBoundaryScheme == 2))) then
                      ! We're either in the interior, or on a boundary point at which trajectories leave the domain,
                      ! so impose the kinetic equation here.
                      rowIndex = getIndex(ispecies,ix,L,itheta,ipsi)
                      colIndex = getIndexExtraSources(iextraSources,ipsi)
                      call MatSetValueSparse(matrix, rowIndex, colIndex, &
                           extraSourceSpeciesPart(iextraSources,ispecies)*extraSourceThetaPart(iextraSources,itheta)* &
                           sourceXPart(ix), ADD_VALUES, ierr)
                   end if
                end do
             end do
          end do
       end do
    end do

    deallocate(sourceXPart)

  end subroutine sources

  subroutine constraints()

    ! *******************************************************************************
    ! Add constraints:
    ! *******************************************************************************

    PetscErrorCode :: ierr
    PetscScalar, dimension(:), allocatable :: constraintXAndThetaPart
    PetscScalar, dimension(:), allocatable :: constraintXPart
    PetscScalar, dimension(:,:), allocatable :: constraintPsiAndSpeciesPart, constraintXandLPart
    PetscScalar, dimension(:,:), allocatable :: constraintPsiAndThetaPart
    
    integer, dimension(:), allocatable :: rowIndices, colIndices
    integer :: ix, itheta, ipsi, L, iL
    integer, dimension(:), allocatable :: Larray
    integer :: ispecies, isources,iextraSources, ispeciesIndepSources
    integer :: rowIndexArray(1)
    integer :: rowIndex, colIndex

    if (procThatHandlesRightBoundary) then
       ! The processor that owns the right-most psi point handles the constraints.

       allocate(colIndices(Ntheta))
       allocate(constraintXAndThetaPart(Ntheta))
       allocate(constraintXPart(Nx))
       allocate(constraintPsiAndThetaPart(Ntheta,Npsi))

       do isources = 1,Nsources
          call initializeDeltaFConstraints(isources,constraintXPart,L)

          do ispecies = 1,Nspecies
             do ipsi = lowestEnforcedIpsi, highestEnforcedIpsi
                rowIndexArray = getIndexSources(isources,ispecies,ipsi)
                do ix = 1, Nx
                   constraintXAndThetaPart = constraintXPart(ix) * thetaWeights / JHat(:,ipsi)
                   colIndices = [(getIndex(ispecies,ix,L,itheta,ipsi), itheta=1,Ntheta)]
                   call MatSetValuesSparse(matrix, 1, rowIndexArray, Ntheta, colIndices, constraintXAndThetaPart,&
                           ADD_VALUES, ierr)
                   CHKERRQ(ierr)
                end do
             end do
          end do
       end do

       do ispeciesIndepSources = 1,NspeciesIndepSources
          ! allocate(constraintPsiAndSpeciesPart(Nspecies,Npsi))
          ! allocate(constraintXandLPart(NL,Nx))
          ! allocate(constraintPsiAndThetaPart(Ntheta,Npsi))
          ! allocate(L(NL))
          call initializeSpeciesIndepDeltaFConstraints(ispeciesIndepSources,constraintXandLPart,Larray,&
               constraintPsiAndSpeciesPart,constraintPsiAndThetaPart)
          do ipsi = 1, Npsi
             rowIndex = getIndexSpeciesIndepSources(ispeciesIndepSources,ipsi)
             do ix = 1, Nx
                do ispecies = 1,Nspecies
                   do iL = lbound(Larray,1), ubound(Larray,1)
                      do itheta = 1,Ntheta
                         colIndex = getIndex(ispecies,ix,Larray(iL),itheta,ipsi)
                         !print *,colIndex
                         call MatSetValueSparse(matrix, rowIndex, colIndex, &
                              constraintPsiAndSpeciesPart(ispecies,ipsi)*constraintXandLPart(iL,ix) &
                              * constraintPsiAndThetaPart(itheta,ipsi),ADD_VALUES, ierr)
                         CHKERRQ(ierr)
                      end do
                   end do
                end do
             end do
          end do    
       end do

       do iextraSources = 1,NextraSources
          allocate(constraintPsiAndSpeciesPart(Nspecies,Npsi))
          call initializeSourceConstraints(iextraSources,constraintPsiAndSpeciesPart,isources)
          do ispecies = 1,Nspecies
             do ipsi = lowestEnforcedIpsi, highestEnforcedIpsi
                rowIndex = getIndexExtraSources(iextraSources,ipsi)
                colIndex = getIndexSources(isources,ispecies,ipsi)
                call MatSetValueSparse(matrix, rowIndex, colIndex, constraintPsiAndSpeciesPart(ispecies,ipsi),&
                     ADD_VALUES, ierr)
                CHKERRQ(ierr)
             end do
          end do
       end do

       ! THIS CODE HAD CONSTRAINTS ON SOURCE CONSTRAINTS SOURCES
       !
       !if (noChargeSource == 3) then
       !   iextraSources = 1
       !   allocate(constraintPsiAndSpeciesPart(Nspecies,Npsi))
       !   call initializeSourceConstraints(iextraSources,constraintPsiAndSpeciesPart)
          ! momentum sources only enters into this constraint
          ! this is reflected in the col-indices used, which should match
          ! the row-indices in the source subroutine
       !   do ispecies = 1,Nspecies
       !      do ipsi = lowestEnforcedIpsi, highestEnforcedIpsi
       !         rowIndex = getIndexExtraSources(iextraSources,ipsi)
       !         colIndex = getIndexExtraSources(iextraSources,ipsi)
       !         call MatSetValueSparse(matrix, rowIndex, colIndex, constraintPsiAndSpeciesPart(ispecies,ipsi),&
       !              ADD_VALUES, ierr)
       !         CHKERRQ(ierr)
       !      end do
       !   end do
       !end if
       
       deallocate(colIndices)
       deallocate(constraintXAndThetaPart)
       deallocate(constraintXPart) 
    end if
  end subroutine constraints

  subroutine finalizeMatrices(whichMatrix,time1)

    ! *******************************************************************************
    ! Done inserting values into the matrices.
    ! Now finalize the matrices:
    ! *******************************************************************************

    integer, intent(in) :: whichMatrix
    PetscLogDouble, intent(inout) :: time1
    PetscErrorCode :: ierr
    PetscLogDouble :: time2

    call PetscTime(time2, ierr)
    if (masterProcInSubComm) then
       if (whichMatrix==0) then
          print *,"[",myCommunicatorIndex,"] Time to pre-assemble preconditioner matrix: ", time2-time1, " seconds."
       else
          print *,"[",myCommunicatorIndex,"] Time to pre-assemble matrix: ", time2-time1, " seconds."
       end if
    end if
    call PetscTime(time1, ierr)

    call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
    if (procThatHandlesLeftBoundary .and. boundaryScheme /= 2  .and. boundaryScheme /= 3) then
       call MatAssemblyBegin(leftMatrix, MAT_FINAL_ASSEMBLY, ierr)
    end if
    if (procThatHandlesRightBoundary .and. boundaryScheme /= 1  .and. boundaryScheme /= 3) then
       call MatAssemblyBegin(rightMatrix, MAT_FINAL_ASSEMBLY, ierr)
    end if
    CHKERRQ(ierr)
    call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
    if (procThatHandlesLeftBoundary .and. boundaryScheme /= 2  .and. boundaryScheme /= 3) then
       call MatAssemblyEnd(leftMatrix, MAT_FINAL_ASSEMBLY, ierr)
    end if
    if (procThatHandlesRightBoundary .and. boundaryScheme /= 1  .and. boundaryScheme /= 3) then
       call MatAssemblyEnd(rightMatrix, MAT_FINAL_ASSEMBLY, ierr)
    end if
    CHKERRQ(ierr)
    
    if (whichMatrix==0) then
       preconditionerMatrix = matrix
       if (procThatHandlesLeftBoundary .and. boundaryScheme /= 2  .and. boundaryScheme /= 3) then
          leftPreconditionerMatrix = leftMatrix
       end if
       if (procThatHandlesRightBoundary .and. boundaryScheme /= 1  .and. boundaryScheme /= 3) then
          rightPreconditionerMatrix = rightMatrix
       end if
    end if

    call PetscTime(time2, ierr)
    if (masterProcInSubComm) then
       if (whichMatrix==0) then
          print *,"[",myCommunicatorIndex,"] Time to assemble preconditioner matrices: ", time2-time1, " seconds."
       else
          print *,"[",myCommunicatorIndex,"] Time to assemble matrices: ", time2-time1, " seconds."
       end if
    end if
    call PetscTime(time1, ierr)

    ! Does not actually print matrix.
    ! Use PETSc's -mat_view instead
    !
    ! if (printMatrixDebug == 1) then
    !    print *,"Printing matrix ", whichMatrix, "..."
    !    print *,matrix
    ! end if

  end subroutine finalizeMatrices

end module DKEMatrices
