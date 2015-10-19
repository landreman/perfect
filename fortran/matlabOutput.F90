module matlabOutput

  use DKEMatrices
  use DKERhs
  use globalVariables

#include <finclude/petsckspdef.h>

  implicit none

contains

  ! *********************************************************
  ! Create a PETSc viewer to record output for Matlab
  ! Note: MatView calls may take a very long time to complete,
  ! As an example a 13 second simulation took 650 seconds to
  ! Write, and the file size ended up being 110MB.
  ! *********************************************************
  subroutine writeMatlabOutput(soln,time1)
    Vec, intent(in) :: soln

    PetscErrorCode :: ierr
    PetscLogDouble, intent(in) :: time1
    PetscLogDouble :: time2
    PetscViewer MatlabOutputViewer

       print *,"[",myCommunicatorIndex,"] Creates Matlab output. Note that this may take a very long time."

    call PetscViewerASCIIOpen(MPIComm, &
         & MatlabOutputFilename,&
         & MatlabOutputViewer, ierr)
    CHKERRQ(ierr)

    call PetscViewerSetFormat(MatlabOutputViewer, PETSC_VIEWER_ASCII_MATLAB, ierr)
    CHKERRQ(ierr)


    call PetscObjectSetName(rhs, "rhs", ierr)
    CHKERRQ(ierr)
    call VecView(rhs, MatlabOutputViewer, ierr)
    CHKERRQ(ierr)
    call PetscObjectSetName(soln, "soln", ierr)
    CHKERRQ(ierr)
    call VecView(soln, MatlabOutputViewer, ierr)
    CHKERRQ(ierr)

    if (useIterativeSolver) then
       call PetscObjectSetName(preconditionerMatrix, "preconditionerMatrix", ierr)
       CHKERRQ(ierr)
       call MatView(preconditionerMatrix, MatlabOutputViewer, ierr)
       CHKERRQ(ierr)
          print *,"[",myCommunicatorIndex,"] Done writing preconditioner matrix. Start writing the main matrix..."
    end if
    call PetscObjectSetName(matrix, "matrix", ierr)
    CHKERRQ(ierr)
    call MatView(matrix, MatlabOutputViewer, ierr)
    CHKERRQ(ierr)



    call PetscTime(time2, ierr)
    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] Time to write output: ", time2-time1, " seconds."
    end if
    call PetscTime(time1, ierr)

    call PetscViewerDestroy(MatlabOutputViewer, ierr)
    CHKERRQ(ierr)

  end subroutine writeMatlabOutput

end module matlabOutput
