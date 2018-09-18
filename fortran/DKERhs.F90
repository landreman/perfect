module DKERhs

  use globalVariables
  use grids
  use sourcesConstraints
  !use petscksp

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#elif (PETSC_VERSION_MAJOR < 3 && PETSC_VERSION_MAJOR<=7)
#include <petsc/finclude/petsckspdef.h>
#else
#include <petsc/finclude/petscksp.h>
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
    integer :: ix, itheta, ipsi, L, index, iextraSources
    integer :: ispecies
    PetscScalar :: LFactor
    PetscScalar :: stuffToAdd
    PetscScalar, dimension(:), allocatable :: this_sourceConstraintsRHS, &
         constantSourceXPart, constantSourceThetaPart

    call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, rhs, ierr)
    CHKERRQ(ierr)
    call VecSet(rhs, zero,ierr)
    CHKERRQ(ierr)
    if (procThatHandlesLeftBoundary .and. boundaryScheme /= 2  .and. boundaryScheme /= 3) then
       ! This process handles the left boundary, so solve for the local solution there.
       call VecCreateSeq(MPI_COMM_SELF, localMatrixSize, rhsLeft, ierr)
       CHKERRQ(ierr)
       call VecSet(rhsLeft, zero,ierr)
       CHKERRQ(ierr)
    end if
    if (procThatHandlesRightBoundary .and. boundaryScheme /= 1  .and. boundaryScheme /= 3) then
       ! This process handles the right boundary, so solve for the local solution there.
       call VecCreateSeq(MPI_COMM_SELF, localMatrixSize, rhsRight, ierr)
       CHKERRQ(ierr)
       call VecSet(rhsRight, zero, ierr)
       CHKERRQ(ierr)
    end if

    do ispecies = 1, Nspecies
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
                !call VecSetValues(rhs, 1, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                call VecSetValue(rhs, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                index = getIndex(ispecies,ix,L,itheta,1)
                if (ipsi==1 .and. boundaryScheme /= 2  .and. boundaryScheme /= 3) then
                   ! This is the left boundary
                   !call VecSetValues(rhsLeft, 1, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                   call VecSetValue(rhsLeft, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                elseif (ipsi==Npsi  .and. boundaryScheme /= 1  .and. boundaryScheme /= 3) then
                   ! This is the right boundary
                   !call VecSetValues(rhsRight, 1, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                   call VecSetValue(rhsRight, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                end if

                L = 2
                LFactor = 2/three
                index = getIndex(ispecies,ix,L,itheta,ipsi)
                !call VecSetValues(rhs, 1, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                call VecSetValue(rhs, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                index = getIndex(ispecies,ix,L,itheta,1)
                if (ipsi==1  .and. boundaryScheme /= 2  .and. boundaryScheme /= 3) then
                   ! This is the left boundary
                   !call VecSetValues(rhsLeft, 1, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                   call VecSetValue(rhsLeft, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                elseif (ipsi==Npsi  .and. boundaryScheme /= 1  .and. boundaryScheme /= 3) then
                   ! This is the right boundary
                   !call VecSetValues(rhsRight, 1, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                   call VecSetValue(rhsRight, index, LFactor*stuffToAdd, ADD_VALUES, ierr)
                end if
                CHKERRQ(ierr)

             end do
          end do
       end do
    end do

    if (NextraSources > 0 .or. NspeciesIndepSources > 0) then
       allocate(this_sourceConstraintsRHS(NEnforcedPsi))
    end if
    
    do iextraSources = 1,NextraSources
       ! Add the RHS of the relation between particle sources
       ! We do not care about boundaries except possibly through the enforced ipsi parameters.
       select case(sourceConstraints(iextraSources))
       case(0)
          this_sourceConstraintsRHS = sourceConstraintsRHS(iextraSources,:) * Delta/(sqrtpi* abs(VPrimeHat) * psiAHat)
       case(1)
          this_sourceConstraintsRHS = sourceConstraintsRHS(iextraSources,:) * two * Delta/(sqrtpi* abs(VPrimeHat) * psiAHat)
       case default
          print *,"Error! Invalid sourceConstraints. Currently supported values are: 0,1."
          stop
       end select
       do ipsi =lowestEnforcedIpsi, highestEnforcedIpsi
          index = getIndexExtraSources(iextraSources,ipsi)
          call VecSetValue(rhs, index, this_sourceConstraintsRHS(ipsi), ADD_VALUES, ierr) 
       end do
    end do

    
    do iextraSources = 1,NspeciesIndepSources
       ! Add the RHS of the constraint on species indep g moments
       ! We do not care about boundaries except possibly through the enforced ipsi parameters.
       select case(speciesIndepGConstraints(iextraSources))
       case(1)
          this_sourceConstraintsRHS = -speciesIndepRHS(iextraSources,:)
       case(2)
          this_sourceConstraintsRHS = -speciesIndepRHS(iextraSources,:)
       case default
          print *,"Error! Invalid sourceConstraints. Currently supported values are: 1,2."
          stop
       end select
       do ipsi =1, Npsi
          index = getIndexSpeciesIndepSources(iextraSources,ipsi)
          call VecSetValue(rhs, index, this_sourceConstraintsRHS(ipsi), ADD_VALUES, ierr) 
       end do
    end do

    
    
    
    if (NconstantSources > 0) then
       allocate(constantSourceXPart(Nx))
       allocate(constantSourceThetaPart(Ntheta))
    end if
    
    do iextraSources = 1,NconstantSources
       ! Add constant sources to the RHS
       ! We do not care about boundaries except possibly through the enforced ipsi parameters. Boundaries will be overwritten later to enforce boundary conditions.
       call initializeConstantSources(iextraSources,constantSourceXPart,L,constantSourceThetaPart)
           
       do ispecies = 1, Nspecies
          do ipsi = ipsiMin, ipsiMax
             do itheta = 1, Ntheta
                do ix = 1, Nx  
                   index = getIndex(ispecies,ix,L,itheta,ipsi)
                   !minus sign since we have moved the sources to the RHS
                   call VecSetValue(rhs, index, -constantSourceThetaPart(itheta) * constantSourceXPart(ix) &
                        * constantSourceProfile(iextraSources,ispecies,ipsi), ADD_VALUES, ierr)
                end do
             end do
          end do
       end do
    end do

    

    ! Assemble the rhs vector here, since we will insert values later
    ! (ADD_VALUE and INSERT_VALUE cannot be used without assembling first)
    call VecAssemblyBegin(rhs, ierr)
    CHKERRQ(ierr)
    call VecAssemblyEnd(rhs, ierr)
    CHKERRQ(ierr)
    
    
  end subroutine DKECreateRhsVector

end module DKERhs
