! Main program

#include <finclude/petscsysdef.h>
#include "PETScVersions.F90"
#include "perfectVersion.h"

program perfect
  use globalVariables
  use geometry
  use printToStdout
  use writeHDF5Output
  use readInput
  use scan
  use solveDKE
  use petscsysdef

  implicit none


  PetscErrorCode ierr
  PetscLogDouble :: time1, time2
  integer :: runNum


  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  CHKERRQ(ierr)

  call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)
  ! This program uses 1-based indices to number the MPI processes, consistent with the Fortran
  ! convention for array indices, not the 0-based indexing used by C and MPI.
  myRank = myRank + 1
  masterProc = (myRank==1)

  call printGreeting()
  
  ! Default values for all input parameters are set in parameters.F90.
  ! Defaults are replaced by any values given in the input namelist file:
  call readNamelistInput("input.namelist")

  select case (programMode)
  case (-1)
     print *,"Switching the sign of programMode to +1 would generate 1 run."
     stop
  case (-2,-3)
     call processConvergenceScanParameters()
     print *,"Switching the sign of programMode would yield a scan with",numRunsInScan," runs."
     stop
  case (-4)
     call setUpUScan()
     print *,"Switching the sign of programMode would yield a U scan with",numRunsInScan," runs."
     stop
  end select

  call chooseParallelDirectSolver()

  ! In case of Miller geometry or a geometry loaded from a file, do any required initialization:
  call initializeGeometry()

  select case (programMode)
  case (1)
     ! Single run

     ! There is no parallelization over runs, so set the sub-commmunicator 
     ! to be equal to the global communicator:
     MPIComm = PETSC_COMM_WORLD
     myRankInSubComm = myRank
     numProcsInSubComm = numProcs
     masterProcInSubComm = masterProc
     myCommunicatorIndex = 1

     call openOutputFile()
     call allocateArraysForSingleRun()
     call setupOutput(GIT_COMMIT)
     call solveDKEMain()
     call writeRunToOutputFile(1)

  case (2, 3, 4)
     ! Convergence scan or U scan

     select case (programMode)
     case (2, 3)
        call processConvergenceScanParameters()
     case (4)
        call setUpUScan()
     end select

     if (masterProc) then
        print *,"Beginning a scan involving ",numRunsInScan," runs."
        print *,"The numbers in brackets below indicate the MPI communicator."
     end if
     call setMPICommunicatorsForScan()
     call openOutputFile()
     call setupOutput(GIT_COMMIT)
     call PetscTime(time1, ierr)

     do runNum = minScanUnit,maxScanUnit
        if (masterProcInSubComm) then
           print *,"[",myCommunicatorIndex,"] --------------------------------------------------------------"
           print *,"[",myCommunicatorIndex,"] Run", runNum, "of", numRunsInScan
        end if
        makeLocalApproximation = makeLocalApproximationsforScan(runNum)
        NpsiPerDiameter = NpsiPerDiametersForScan(runNum)
        Npsi = NpsisForScan(runNum)
        psiDiameter = psiDiametersForScan(runNum)
        widthExtender = widthExtendersForScan(runNum)
        Ntheta = NthetasForScan(runNum)
        Nxi = NxisForScan(runNum)
        NL = NLsForScan(runNum)
        Nx = NxsForScan(runNum)
        NxPotentialsPerVth = NxPotentialsPerVthsForScan(runNum)
        xMax = xMaxsForScan(runNum)
        solverTolerance = solverTolerancesForScan(runNum)
        desiredU = desiredUsForScan(runNum)
        call solveDKEmain()
        call writeRunToOutputFile(runNum)
        call deallocateArrays()
     end do
     call PetscTime(time2, ierr)
     if (masterProcInSubComm) then
        print *,"[",myCommunicatorIndex,"] --------------------------------------------------------------"
        print *,"[",myCommunicatorIndex,"] Total time for scan on this communicator: ", &
             time2-time1, "seconds."
     end if

  case default
     if (masterProc) then
        print *,"Error: invalid programMode"
     end if
     stop
  end select

  call closeOutputFile()
  if (masterProc) then
     print *,"Output has been written. Will now crash at PetscFinalize."
  end if
  call PetscFinalize(ierr)

end program perfect
