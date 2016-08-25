module indices

  use globalVariables

  implicit none

  integer, dimension(:), allocatable :: first_index_for_x

contains

  integer function getIndex(iSpecies,ix,L,itheta,ipsi)

    integer,intent(in) :: iSpecies,ix,L,itheta,ipsi 

    ! Validate inputs:

    if (iSpecies < 1) then
      print *,"Error: iSpecies < 1"
      stop
    end if

    if (iSpecies > numSpecies) then
      print *,"Error: iSpecies > numSpecies"
      stop
    end if

    if (ix < 1) then
      print *,"Error: ix < 1"
      stop
    end if

    if (ix > Nx) then
      print *,"Error: ix > Nx"
      stop
    end if

    if (L < 0) then
      print *,"Error: L < 0"
      stop
    end if

    if (L > Nxi_for_x(ix)-1) then
      print *,"Error: L > Nxi_for_x(ix)-1"
      print *,"L:",L
      print *,"ix:",ix
      stop
    end if

    if (itheta < 1) then
      print *,"Error: itheta < 1"
      stop
    end if

    if (itheta > Ntheta) then
      print *,"Error: itheta > Ntheta"
      stop
    end if

    if (ipsi < 1) then
      print *,"Error: ipsi < 1"
      stop
    end if

    if (ipsi > Npsi) then
      print *,"Error: ipsi > Npsi"
      stop
    end if

    ! Done with validation.

    getIndex = (ipsi-1)*localMatrixSize + (ispecies-1)*localDKEMatrixSize &
        + first_index_for_x(ix)*Ntheta + L*Ntheta + itheta - 1

  end function getIndex

  subroutine calculate_first_index_for_x

    integer :: ix

    allocate(first_index_for_x(Nx))
    first_index_for_x(1)=0
    do ix=2,Nx
      first_index_for_x(ix) = first_index_for_x(ix-1)+Nxi_for_x(ix-1)
    end do

  end subroutine calculate_first_index_for_x

end module indices
