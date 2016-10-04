module DKERhs

  use globalVariables
  use grids
  !use petscksp

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

  implicit none

  Vec :: rhs, rhsLeft, rhsRight

contains

  ! *******************************************************************************
  ! *******************************************************************************
  !
  ! Create the right-hand side vector
  !
  ! *******************************************************************************
  ! *******************************************************************************
  subroutine DKECreateRhsVector()

    PetscErrorCode :: ierr
    integer :: ix, itheta, ipsi, L, index
    integer :: ispecies
    PetscScalar :: LFactor
    PetscScalar :: stuffToAdd
    print *,matrixSize
    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, rhs, ierr)
    CHKERRQ(ierr)
    call VecSet(rhs, zero,ierr)
    CHKERRQ(ierr)
    if (procThatHandlesLeftBoundary) then
       ! This process handles the left boundary, so solve for the local solution there.
       call VecCreateSeq(MPI_COMM_SELF, localMatrixSize, rhsLeft, ierr)
       CHKERRQ(ierr)
       call VecSet(rhsLeft, zero,ierr)
    end if
    if (procThatHandlesRightBoundary) then
       ! This process handles the right boundary, so solve for the local solution there.
       call VecCreateSeq(MPI_COMM_SELF, localMatrixSize, rhsRight, ierr)
       CHKERRQ(ierr)
       call VecSet(rhsRight, zero, ierr)
    end if
    CHKERRQ(ierr)

    do ispecies = 1, numSpecies
       do ipsi = ipsiMin, ipsiMax
          do itheta = 1, Ntheta
             do ix = 1, Nx

                stuffToAdd = masses(ispecies)*masses(ispecies) * nHats(ispecies,ipsi) * IHat(ipsi) &
                     * JHat(itheta,ipsi) * dBHatdtheta(itheta,ipsi) * x2(ix) * expx2(ix) &
                     /(2*pi*sqrtpi*charges(ispecies)*sqrtTHats(ispecies,ipsi)*psiAHatArray(ipsi) &
                     * (BHat(itheta,ipsi) ** 3)) &
                     * (dnHatdpsis(ispecies,ipsi)/nHats(ispecies, ipsi) &
                     + 2*charges(ispecies)/THats(ispecies,ipsi)*omega/Delta*dPhiHatdpsi(ipsi) &
                     + (x2(ix)-3/two)/THats(ispecies,ipsi)*dTHatdpsis(ispecies,ipsi))

                ! It's okay to assign the rhs even at ipsi=1 and ipsi=Npsi because we will
                ! over-write these values later.

                L = 0
                LFactor = 4/three
                index = getIndex(ispecies,ix,L,itheta,ipsi)
                !call VecSetValues(rhs, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                call VecSetValue(rhs, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                index = getIndex(ispecies,ix,L,itheta,1)
                if (ipsi==1) then
                   ! This is the left boundary
                   !call VecSetValues(rhsLeft, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                   call VecSetValue(rhsLeft, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                elseif (ipsi==Npsi) then
                   ! This is the right boundary
                   !call VecSetValues(rhsRight, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                   call VecSetValue(rhsRight, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                end if

                L = 2
                LFactor = 2/three
                index = getIndex(ispecies,ix,L,itheta,ipsi)
                !call VecSetValues(rhs, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                call VecSetValue(rhs, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                index = getIndex(ispecies,ix,L,itheta,1)
                if (ipsi==1) then
                   ! This is the left boundary
                   !call VecSetValues(rhsLeft, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                   call VecSetValue(rhsLeft, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                elseif (ipsi==Npsi) then
                   ! This is the right boundary
                   !call VecSetValues(rhsRight, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                   call VecSetValue(rhsRight, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                end if
                CHKERRQ(ierr)

             end do
          end do
       end do
    end do

    if ((noChargeSource == 1) .or. (noChargeSource == 2)) then
       ! Add the RHS of the relation between particle sources
       ! We do not care about boundaries except possibly through the enforced ipsi parameters.
       do ipsi =lowestEnforcedIpsi, highestEnforcedIpsi
          index = Npsi * localMatrixSize + NEnforcedPsi * Nsources * numSpecies + (ipsi - lowestEnforcedIpsi)
          ! prefactors needed to translate to source that gives psiN derivative of particleFlux output
          call VecSetValue(rhs, index, chargeSource(ipsi) * Delta/(sqrtpi* abs(VPrimeHat) * psiAHat), INSERT_VALUES, ierr)
          ! 
       end do
    end if
    ! 
       
  end subroutine DKECreateRhsVector

end module DKERhs
