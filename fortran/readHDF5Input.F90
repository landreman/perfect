module readHDF5Input

  use HDF5

  implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

  interface readVariable
    module procedure readVariable_integer
    module procedure readVariable_scalar
    module procedure readVariable_1d
    module procedure readVariable_2d
  end interface readVariable

  integer(HID_T), private :: fileID, groupID

contains

  subroutine openInputFile(fileName, groupName)

    character(len=*), intent(in) :: fileName, groupName
    integer :: HDF5Error

    ! Open input file
    call h5fopen_f(trim(fileName), H5F_ACC_RDONLY_F, fileID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening file:", fileName
      stop
    end if
       
    ! Open group
    call h5gopen_f(fileID, trim(groupName), groupID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening group:",groupName
      stop
    end if

  end subroutine openInputFile

  subroutine closeInputFile()

    integer :: HDF5Error

    ! Close the group
    call h5gclose_f(groupID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing input group"
      stop
    end if

    ! Close the file
    call h5fclose_f(fileID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing input file"
      stop
    end if

  end subroutine closeInputFile

  subroutine readVariable_integer(variable,varname)

    integer, intent(inout) :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(0) :: dimensions
    integer :: HDF5Error

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable:",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_INTEGER, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable:",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable:",varname
      stop
    end if

  end subroutine readVariable_integer

  subroutine readVariable_scalar(variable,varname)

    PetscScalar, intent(inout) :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(0) :: dimensions
    integer :: HDF5Error

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable:",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable:",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable:",varname
      stop
    end if

  end subroutine readVariable_scalar

  subroutine readVariable_1d(variable,varname)

    PetscScalar, intent(inout), dimension(:), allocatable :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(1) :: dimensions
    integer :: HDF5Error

    if (.not. allocated(variable)) then
      print *,"Tried to read into unallocated array:",varname
      stop
    end if

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable:",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable:",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable:",varname
      stop
    end if

  end subroutine readVariable_1d

  subroutine readVariable_2d(variable,varname)

    PetscScalar, intent(inout), dimension(:,:), allocatable :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(2) :: dimensions
    integer :: HDF5Error

    if (.not. allocated(variable)) then
      print *,"Tried to read into unallocated array:",varname
      stop
    end if

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable:",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable:",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable:",varname
      stop
    end if

  end subroutine readVariable_2d

end module readHDF5Input

