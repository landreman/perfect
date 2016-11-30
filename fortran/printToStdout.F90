module printToStdout

  use globalVariables
  use petscsysdef

  implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

contains

  ! -----------------------------------------------------------------------------------

  subroutine printGreeting()

    implicit none

    if (masterProc) then
       print *,"****************************************************************************"
       print *,"PERFECT: Pedestal and Edge Radially-global Fokker-Planck Evaluation of Collisional Transport."
#if defined(PETSC_USE_REAL_SINGLE)
       print *,"Using single precision."
#else
       print *,"Using double precision."
#endif
       if (numProcs==1) then
          print *,"Serial job (1 process) detected."
       else
          print "(a, i4, a)", " Parallel job (",numProcs," processes) detected."
       end if
    end if
  end subroutine printGreeting

  ! -----------------------------------------------------------------------------------

  subroutine printInputs()

    implicit none

    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] ---- Physics parameters: ----"
       print *,"[",myCommunicatorIndex,"] number of species = ",Nspecies
       print *,"[",myCommunicatorIndex,"] masses = ",masses(1:Nspecies)
       print *,"[",myCommunicatorIndex,"] charges = ",charges(1:Nspecies)
       print *,"[",myCommunicatorIndex,"] desiredU = ", desiredU
       print *,"[",myCommunicatorIndex,"] epsilon = ", epsil
!       print *,"[",myCommunicatorIndex,"] nuStar  = ", nuStar
       print *,"[",myCommunicatorIndex,"] nu_r = ", nu_r
       print *,"[",myCommunicatorIndex,"] Geometry scheme: ", geometryToUse
       if (makeLocalApproximation) then
          print *,"[",myCommunicatorIndex,"] Making the local approximation"
       else
          print *,"[",myCommunicatorIndex,"] Not making the local approximation"
       end if
       print *,"[",myCommunicatorIndex,"] ---- Numerical parameters: ----"
       print *,"[",myCommunicatorIndex,"] NpsiPerDiameter    = ", NpsiPerDiameter
       print *,"[",myCommunicatorIndex,"] psiDiameter        = ", psiDiameter
       print *,"[",myCommunicatorIndex,"] widthExtender      = ", widthExtender
       print *,"[",myCommunicatorIndex,"] Npsi               = ", Npsi
       print *,"[",myCommunicatorIndex,"] Ntheta             = ", Ntheta
       print *,"[",myCommunicatorIndex,"] Nxi                = ", Nxi
       print *,"[",myCommunicatorIndex,"] NL                 = ", NL
       print *,"[",myCommunicatorIndex,"] Nx                 = ", Nx
       print *,"[",myCommunicatorIndex,"] NxPotentialsPerVth = ", NxPotentialsPerVth
       print *,"[",myCommunicatorIndex,"] xMax               = ",xMax
       print *,"[",myCommunicatorIndex,"] solverTolerance    = ",solverTolerance
       print *,"[",myCommunicatorIndex,"] xScaleFactor       = ",xScaleFactor
       select case (psiDerivativeScheme)
       case (0)
          print *,"[",myCommunicatorIndex,"] Psi derivative: Chebyshev"
       case (1)
          print *,"[",myCommunicatorIndex,"] Psi derivative: centered finite differences, 3-point stencil"
       case (2)
          print *,"[",myCommunicatorIndex,"] Psi derivative: centered finite differences, 5-point stencil"
       case (3)
          print *,"[",myCommunicatorIndex,"] Psi derivative: upwinded finite differences, 2-point stencil"
       case (4)
          print *,"[",myCommunicatorIndex,"] Psi derivative: upwinded finite differences, 3-point stencil"
       case default
          print *,"[",myCommunicatorIndex,"] Error! Invalid setting for psiDerivativeScheme"
          stop
       end select
       select case (thetaDerivativeScheme)
       case (0)
          print *,"[",myCommunicatorIndex,"] Theta derivative: spectral collocation"
       case (1)
          print *,"[",myCommunicatorIndex,"] Theta derivative: centered finite differences, 3-point stencil"
       case (2)
          print *,"[",myCommunicatorIndex,"] Theta derivative: centered finite differences, 5-point stencil"
       case default
          print *,"[",myCommunicatorIndex,"] Error! Invalid setting for thetaDerivativeScheme"
          stop
       end select
       if (useIterativeSolver) then
          print *,"[",myCommunicatorIndex,"] Using iterative solver"
       else
          print *,"[",myCommunicatorIndex,"] Using direct solver"
       end if
    end if
  end subroutine printInputs

  ! -----------------------------------------------------------------------------------

  subroutine printOutputs()

    implicit none

    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] Total elapsed time: ", elapsedTime, " seconds."
!!$       print *,"[",myCommunicatorIndex,"] q = ", q
!!$       print *,"[",myCommunicatorIndex,"] particleFlux = ", particleFlux
!!$       print *,"[",myCommunicatorIndex,"] k = ", k
    end if

  end subroutine printOutputs

end module printToStdout

