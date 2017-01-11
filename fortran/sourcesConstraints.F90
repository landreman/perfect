module sourcesConstraints

  use globalVariables
  use grids

  implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

contains

  subroutine initializeSources(isources,sourceXPart,L,sourceThetaPart,sourceThetaPartFSA)

    ! Here we must calculate the following variables:
    ! sourceXPart(Nx)
    ! sourceThetaPart(Ntheta)
    ! sourceSpeciesPart(Nspecies)
    ! sourceThetaPartFSA(Npsi)
    ! L

    integer, intent(in) :: isources
    PetscScalar, dimension (:), intent(out) :: sourceXPart, sourceThetaPart, sourceThetaPartFSA
    integer, intent(out) :: L

    call initializeSourcesVPart(isources,sourceXPart,L)
    
    call initializeSourcesThetaPart(isources,sourceThetaPart,sourceThetaPartFSA)
    
  end subroutine initializeSources

  subroutine initializeSourcesVPart(isources,sourceXPart,L)
    integer, intent(in) :: isources
    PetscScalar, intent(out) :: sourceXPart(:)
    integer, intent(out) :: L
    select case(sourcesVStructure(isources))
    case(1)
       iparticleSource = isources
       L = 0
       sourceXPart = (x2-5/two)*exp(-x2)
       

    case(2)
       iheatSource = isources
       L = 0
       sourceXPart = (x2-3/two)*exp(-x2)

    case(3)
       imomentumSource = isources
       L = 1
       sourceXPart = x*(x2-7/two)*exp(-x2)
       
    case default
       print *,"Error! Invalid sourcesVStructure. 1,2,3 implemented"
       stop
    end select

    
  end subroutine initializeSourcesVPart
  
  subroutine initializeSourcesThetaPart(isources,sourceThetaPart,sourceThetaPartFSA)
    integer, intent(in) :: isources
    PetscScalar, intent(out) :: sourceThetaPartFSA(:), sourceThetaPart(:)
    integer i

    select case(sourcesThetaStructure(isources))
    case (0)
       do i = 1,Ntheta
          sourceThetaPart(i) = one
       end do
    case (1)
       do i=1,Ntheta
          sourceThetaPart(i) = one + sourcePoloidalVariationStrength * cos(theta(i) + sourcePoloidalVariationPhase)           
       end do
    case default
       print *,"Error! Invalid sourcesThetaStructure. 0,1 implemented"
       stop
    end select

    ! FSA of theta structure
    do i=lowestEnforcedIpsi,highestEnforcedIpsi
       sourceThetaPartFSA(i - lowestEnforcedIpsi + 1) = dot_product(thetaWeights, sourceThetaPart/JHat(:,i)) / VPrimeHat(i)
    end do


  end subroutine initializeSourcesThetaPart

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  deltaF constraints
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine initializeDeltaFConstraints(isources,constraintXPart,L)
    ! Here we must calculate the following variables:
    ! constraintXPart(Nx)
    ! L

    integer, intent(in) :: isources
    PetscScalar, dimension(:), intent(out), allocatable :: constraintXPart
    integer, intent(out) :: L
    allocate(constraintXPart(Nx))

    select case(gConstraints(isources))
    case(1)
       ! Enforce <n_1> = 0 at psi between lowestEnforcedIpsi and highestEnforcedIpsi
       L=0
       constraintXPart = xWeights*x2
    case(2)
       ! Enforce <p_1> = 0 at psi between lowestEnforcedIpsi and highestEnforcedIpsi
       L=0
       constraintXPart = xWeights*x2*x2
    case default
       print *,"Error! Invalid gConstraints. 1,2 implemented "
       stop
    end select
    
  end subroutine initializeDeltaFConstraints

  subroutine initializeSpeciesIndepDeltaFConstraints(isources,constraintXandLPart,L,constraintPsiAndSpeciesPart,&
       constraintPsiAndThetaPart)
    ! Here we must calculate the following variables:
    ! constraintXPart(Nx)
    ! L
    ! constraintPsiPart

    integer, intent(in) :: isources
    PetscScalar, dimension(:,:), intent(out), allocatable :: constraintPsiAndSpeciesPart, constraintXandLPart
    PetscScalar, dimension(:,:), intent(out), allocatable :: constraintPsiAndThetaPart
    integer, dimension(:), intent(out), allocatable :: L
    integer :: NL, iL, ispecies, ipsi

    allocate(constraintPsiAndSpeciesPart(Nspecies,Npsi))
    allocate(constraintPsiAndThetaPart(Ntheta,Npsi))
    
    select case(speciesIndepGConstraints(isources))
    case(1)
       ! Enforce j_r = 0
       NL = 2
       allocate(L(NL))
       allocate(constraintXandLPart(NL,Nx))
       
       
       L= (/ 0, 2 /)
       constraintXandLPart(1,:) = (8/three) * xWeights*x2*x2 !L=0
       constraintXandLPart(2,:) = (four/15) * xWeights*x2*x2 !L=2
       
       do ispecies = 1,Nspecies
          ! extra charges since we want current to add up to zero
          constraintPsiAndSpeciesPart(ispecies,:) = &
               -masses(ispecies)  * IHat * ((THats(ispecies,:)/masses(ispecies)) ** (5/two))
       end do

       do ipsi = 1,Npsi
          constraintPsiAndThetaPart(:,ipsi) = dBHatdtheta(:,ipsi) / (BHat(:,ipsi)**3)*thetaWeights
       end do
    
    case default
       print *,"Error! Invalid speciesIndepGConstraints. 1 implemented. "
       stop
    end select
    
  end subroutine initializeSpeciesIndepDeltaFConstraints

  subroutine initializeSourceConstraints(isources,constraintPsiAndSpeciesPart, whichSource)
    integer, intent(in) :: isources
    PetscScalar, dimension(:,:), intent(out) :: constraintPsiAndSpeciesPart
    integer, intent(out) :: whichSource
    PetscScalar :: sourceThetaPartIHatOverBHatFSA, constraintSpeciesPart
    integer :: ispecies, ipsi
    ! Here we must calculate the following variables:
    ! constraintPsiAndSpeciesPart(Nspecies,Npsi)

    select case(sourceConstraints(isources))
    case(0)
       if (iparticleSource /= sourcesNotInitialized) then
          whichSource = iparticleSource
       else
          print *,"Error! Constraints on particle sources enforced, but no particle source are initialized."
          stop
       end if
       
       do ispecies = 1,Nspecies
          constraintSpeciesPart = charges(ispecies)/(masses(ispecies)**2)
          do ipsi = lowestEnforcedIpsi, highestEnforcedIpsi
             constraintPsiAndSpeciesPart(ispecies,ipsi) = constraintSpeciesPart * sourceThetaPartFSA(isources,ipsi) &
                  * THats(ispecies,ipsi)**(3.0/2.0)
          end do
       end do
       
       
    case(1)
       if (imomentumSource /= sourcesNotInitialized) then
          whichSource = imomentumSource
       else
          print *,"Error! Constraints on momentum sources enforced, but no momentum source are initialized."
          stop
       end if
       
       do ispecies = 1,Nspecies
          constraintSpeciesPart = extraSourceSpeciesPart(isources,ispecies)/(masses(ispecies)**(3.0/2.0))
          do ipsi = lowestEnforcedIpsi, highestEnforcedIpsi
             sourceThetaPartIHatOverBHatFSA = &
                  dot_product(thetaWeights, sourceThetaPart(isources,:)*IHat(ipsi)/(BHat(:,ipsi)*JHat(:,ipsi))) / VPrimeHat(ipsi)
             constraintPsiAndSpeciesPart(ispecies,ipsi) = constraintSpeciesPart * sourceThetaPartIHatOverBHatFSA &
                  * THats(ispecies,ipsi)**(2.0)
          end do
       end do
    case default
       print *,"Error! Invalid sourceConstraints. 0,1 implemented"
       stop
    end select
    
  end subroutine initializeSourceConstraints


  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  Extra sources
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine initializeExtraSources(isources,&
       sourceXPart,L,sourceThetaPart,sourceThetaPartFSA,sourceSpeciesPart)

    ! Here we must calculate the following variables:
    ! sourceXPart(Nx)
    ! sourceThetaPart(Ntheta)
    ! sourceSpeciesPart(Nspecies)
    ! sourceThetaPartFSA(Npsi)
    ! L

    integer, intent(in) :: isources
    PetscScalar, dimension (:), intent(out) :: sourceXPart, sourceThetaPart, &
         sourceSpeciesPart, sourceThetaPartFSA
    integer, intent(out) :: L    

    call initializeExtraSourcesVPart(isources,sourceXPart,L)
    
    call initializeExtraSourcesThetaPart(isources,sourceThetaPart,sourceThetaPartFSA)

    call initializeExtraSourcesSpeciesPart(isources,sourceSpeciesPart)
    
  end subroutine initializeExtraSources

  subroutine initializeExtraSourcesVPart(isources,sourceXPart,L)
    integer, intent(in) :: isources
    PetscScalar, intent(out) :: sourceXPart(:)
    integer, intent(out) :: L
    select case(extraSourcesVStructure(isources))
    case(1)
       L = 0
       sourceXPart = (x2-5/two)*exp(-x2)

    case(2)
       L = 0
       sourceXPart = (x2-3/two)*exp(-x2)

    case(3)
       L = 1
       sourceXPart = x*(x2-7/two)*exp(-x2)
       
    case default
       print *,"Error! Invalid extraSourcesVStructure. 1,2,3 implemented"
       stop
    end select

    
  end subroutine initializeExtraSourcesVPart
  
  subroutine initializeExtraSourcesThetaPart(isources,sourceThetaPart,sourceThetaPartFSA)
    integer, intent(in) :: isources
    PetscScalar, intent(out) :: sourceThetaPartFSA(:), sourceThetaPart(:)
    integer i

    select case(extraSourcesThetaStructure(isources))
    case (0)
       do i = 1,Ntheta
          sourceThetaPart(i) = one
       end do
    case (1)
       do i=1,Ntheta
          sourceThetaPart(i) = one + sourcePoloidalVariationStrength * cos(theta(i) + sourcePoloidalVariationPhase)           
       end do
    case default
       print *,"Error! Invalid extraSourcesThetaStructure. 0,1 implemented"
       stop
    end select

    ! FSA of theta structure
    do i=lowestEnforcedIpsi,highestEnforcedIpsi
       sourceThetaPartFSA(i - lowestEnforcedIpsi + 1) = dot_product(thetaWeights, sourceThetaPart/JHat(:,i)) / VPrimeHat(i)
    end do


  end subroutine initializeExtraSourcesThetaPart

  subroutine initializeExtraSourcesSpeciesPart(isources,sourceSpeciesPart)
    integer, intent(in) :: isources
    PetscScalar, intent(out) :: sourceSpeciesPart(:)
    
    select case(extraSourcesSpeciesStructure(isources))
    case(0)
       sourceSpeciesPart = masses(1:Nspecies)
    case(1)
       sourceSpeciesPart = zero
       sourceSpeciesPart(1) = one
    case(2)
       sourceSpeciesPart = masses(1:Nspecies)*nHats(:,1)
    case(3)
       sourceSpeciesPart = masses(1:Nspecies)**(1.5) &
       			 * nHats(:,1)/charges(1:Nspecies)				
    case(4)
       sourceSpeciesPart = masses(1:Nspecies) &
       			 * nHats(:,1)/charges(1:Nspecies)
    
    case default
       print *,"Error! Invalid extraSourcesSpeciesStructure. Currently supported values are: 0,1,2,3."
       stop
    end select


  end subroutine initializeExtraSourcesSpeciesPart

  !
  ! INITIALIZE CONSTANT SOURCES
  !

  subroutine initializeConstantSources(isources,sourceXPart,L,sourceThetaPart)

    ! Here we must calculate the following variables:
    ! sourceXPart(Nx)
    ! sourceThetaPart(Ntheta)
    ! L

    integer, intent(in) :: isources
    PetscScalar, dimension (:), intent(out) :: sourceXPart, sourceThetaPart
    integer, intent(out) :: L

    call initializeConstantSourcesVPart(isources,sourceXPart,L)
    
    call initializeConstantSourcesThetaPart(isources,sourceThetaPart)
    
  end subroutine initializeConstantSources

  subroutine initializeConstantSourcesVPart(isources,sourceXPart,L)
    integer, intent(in) :: isources
    PetscScalar, intent(out) :: sourceXPart(:)
    integer, intent(out) :: L
    select case(constantSourcesVStructure(isources))
    case(1)
       ! iconstantparticleSource = isources
       L = 0
       sourceXPart = (x2-5/two)*exp(-x2)
       
    case(2)
       ! iconstantheatSource = isources
       L = 0
       sourceXPart = (x2-3/two)*exp(-x2)

    case(3)
       ! iconstantmomentumSource = isources
       L = 1
       sourceXPart = x*(x2-7/two)*exp(-x2)
       
    case default
       print *,"Error! Invalid sourcesVStructure. 1,2,3 implemented"
       stop
    end select

    
  end subroutine initializeConstantSourcesVPart
  
  subroutine initializeConstantSourcesThetaPart(isources,sourceThetaPart)
    integer, intent(in) :: isources
    PetscScalar, intent(out) :: sourceThetaPart(:)
    integer i

    select case(constantSourcesThetaStructure(isources))
    case (0)
       do i = 1,Ntheta
          sourceThetaPart(i) = one
       end do
    case (1)
       do i=1,Ntheta
          sourceThetaPart(i) = one + sourcePoloidalVariationStrength * cos(theta(i) + sourcePoloidalVariationPhase)           
       end do
    case default
       print *,"Error! Invalid sourcesThetaStructure. 0,1 implemented"
       stop
    end select

    ! FSA of theta structure
    !do i=1,Npsi
    !   sourceThetaPartFSA(i) = dot_product(thetaWeights, sourceThetaPart/JHat(:,i)) / VPrimeHat(i)
    ! end do


  end subroutine initializeConstantSourcesThetaPart

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  speciesIndep sources
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine initializeSpeciesIndepSources(isources,&
       sourceXPart,L,sourceThetaPart,sourceThetaPartFSA,sourceSpeciesPart)

    ! Here we must calculate the following variables:
    ! sourceXPart(Nx)
    ! sourceThetaPart(Ntheta)
    ! sourceSpeciesPart(Nspecies)
    ! sourceThetaPartFSA(Npsi)
    ! L

    integer, intent(in) :: isources
    PetscScalar, dimension (:), intent(out) :: sourceXPart, sourceThetaPart, &
         sourceSpeciesPart, sourceThetaPartFSA
    integer, intent(out) :: L    

    call initializeSpeciesIndepSourcesVPart(isources,sourceXPart,L)
    
    call initializeSpeciesIndepSourcesThetaPart(isources,sourceThetaPart,sourceThetaPartFSA)

    call initializeSpeciesIndepSourcesSpeciesPart(isources,sourceSpeciesPart)
    
  end subroutine initializeSpeciesIndepSources

  subroutine initializeSpeciesIndepSourcesVPart(isources,sourceXPart,L)
    integer, intent(in) :: isources
    PetscScalar, intent(out) :: sourceXPart(:)
    integer, intent(out) :: L
    select case(speciesIndepSourcesVStructure(isources))
    case(1)
       L = 0
       sourceXPart = (x2-5/two)*exp(-x2)

    case(2)
       L = 0
       sourceXPart = (x2-3/two)*exp(-x2)

    case(3)
       L = 1
       sourceXPart = x*(x2-7/two)*exp(-x2)
       
    case default
       print *,"Error! Invalid speciesIndepSourcesVStructure. 1,2,3 implemented"
       stop
    end select

    
  end subroutine initializeSpeciesIndepSourcesVPart
  
  subroutine initializeSpeciesIndepSourcesThetaPart(isources,sourceThetaPart,sourceThetaPartFSA)
    integer, intent(in) :: isources
    PetscScalar, intent(out) :: sourceThetaPartFSA(:), sourceThetaPart(:)
    integer i

    select case(speciesIndepSourcesThetaStructure(isources))
    case (0)
       do i = 1,Ntheta
          sourceThetaPart(i) = one
       end do
    case (1)
       do i=1,Ntheta
          sourceThetaPart(i) = one + sourcePoloidalVariationStrength * cos(theta(i) + sourcePoloidalVariationPhase)           
       end do
    case default
       print *,"Error! Invalid speciesIndepSourcesThetaStructure. 0,1 implemented"
       stop
    end select

    ! FSA of theta structure
    do i=lowestEnforcedIpsi,highestEnforcedIpsi
       sourceThetaPartFSA(i - lowestEnforcedIpsi + 1) = dot_product(thetaWeights, sourceThetaPart/JHat(:,i)) / VPrimeHat(i)
    end do


  end subroutine initializeSpeciesIndepSourcesThetaPart

  subroutine initializeSpeciesIndepSourcesSpeciesPart(isources,sourceSpeciesPart)
    integer, intent(in) :: isources
    PetscScalar, intent(out) :: sourceSpeciesPart(:)
    
    select case(speciesIndepSourcesSpeciesStructure(isources))
    case(0)
       sourceSpeciesPart = masses(1:Nspecies)
    case(1)
       sourceSpeciesPart = zero
       sourceSpeciesPart(1) = one
    case(2)
       sourceSpeciesPart = masses(1:Nspecies)*nHats(:,1)
    case(3)
       sourceSpeciesPart = nHats(:,1)*masses(1:Nspecies)**(1.5)/charges(1:Nspecies)
    case default
       print *,"Error! Invalid speciesIndepSourcesSpeciesStructure. Currently supported values are: 0,1,2,3."
       stop
    end select


  end subroutine initializeSpeciesIndepSourcesSpeciesPart


  
end module
