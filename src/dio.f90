!==============================================================================!
! PROGRAM DIO                                                                  !
!                                                                              !
! This program ...                                                             !
!                                                                              !
! Licence:                                                                     !
! Article describing the code:                                                 !
! GitHub repository:                                                           !
!==============================================================================!
PROGRAM DIO
use Globals
use Inout, only: set_output_filename
use Preparation, only: do_preparation
use Wavefunctions, only: choose_basis_wavefuction
use Field, only: generate_initial_fields
use Matrix, only: calculate_sigma_nabla,calculate_meson_propagators
implicit none

integer :: index
call do_preparation(.True.)
call choose_basis_wavefuction(.True.)
call generate_initial_fields(.True.)
call calculate_sigma_nabla(.True.)
call calculate_meson_propagators(.True.)
! do index = 1, constraint%length
!     constraint%index = index
!     ! input_constraint中的betac、bet3c 其实就是woodssaxon中的 beta2、beta3 ? Not exactly.
!     if(constraint%icstr > 0) then
!         woodssaxon%beta2 = constraint%betac(index) ! set beta2 of WoodsSaxon Potential as beta2 of constraint
!     end if
!     call set_output_filename(constraint%betac(index),constraint%bet3c(index))
!     ! print *,index,outputfile
! end do

! print *, woodssaxon
! print *, basis_harosc
! print *, force
END PROGRAM DIO

!==============================================================================!
! End of file                                                                  !
!==============================================================================!