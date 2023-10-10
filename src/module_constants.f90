!==============================================================================!
! MODULE Constants                                                             !
!                                                                              !
! This module contains the variables related to the physical, mathematical and !
! numerical constants.                                                         !
!==============================================================================!
MODULE Constants

use Iso_fortran_env

implicit none
public
!!! Definition of input/output units (for portability)
integer, parameter :: u_input  = input_unit,  &
                      u_output = output_unit, &
                      u_start  = u_input + u_output, &!the beginning of the unit of files
                      u_config = u_start + 999 ! unit of config.out

!!! Definition of kind parameters (for portability)
integer, parameter :: i8  = int8,   & ! integer  8 bits
                      i16 = int16,  & ! integer 16 bits 
                      i32 = int32,  & ! integer 32 bits 
                      i64 = int64,  & ! integer 64 bits 
                      r32 = real32, & ! real 32 bits (single precision)
                      r64 = real64    ! real 64 bits (double precision)

!!! Definition of simple names for numerical values
real(r64), parameter :: zero  = 0.0d0, &
                        one   = 1.0d0, &
                        two   = 2.0d0, &
                        half  = 0.5d0, &
                        third = 0.333333333333333333d0
                        


!!! Definition of physical constants
real(r64), parameter :: pi  = 4.0d0 * atan(1.0d0), & ! 3.14159265...
                        hbc = 197.328284d0,        & ! hbar*c in Mev.fm
                        radius_r0 = 1.2d0          ! radius factor
                        !alphi = 137.03602d0       & 
                        ! mass_mp = 938.27208816d0, & ! proton mass
                        ! mass_mn = 939.56542052d0, & ! neutron mass
                        ! mass_ma = (mass_mp + mass_mn)/2, & ! nucleon mass
                        ! hbarmass = hbarc**2 / (2*mass_ma)  ! factor kin. energy

!!! Constants used in some subroutine or functions(to declare the length of arrays)
! maximal number for GFV                                                      
integer, parameter :: igfv   = 100, &
                      igfvbc = 100                                            
! number of gauss-meshpoints
integer, parameter :: ngh   = 32,&
                      ngl   = 16,&
                      nghl  = ngh * ngl 

! number for basis
integer, parameter :: nt_max = 2023, &! ntx! maximal number of all levels for protons or neutrons(fermions)
                      ntb_max = 132, & !nox !maximal number of all levels for bosons
                      nb_max = 42, &! nbx! maximal number of nb(nb is the number of K-parity-blocks)
                      nz_max = 21, &!nzx ! maximal nz-quantum number for fermions
                      nr_max = 10, &!nrx ! maximal nr-quantum number for fermions
                      ml_max = 21, &!mlx ! maximal ml-quantum number for fermions 
                      nz_max_boson= 20, &!nzbx! maximal nz-quantum number for bosons  
                      nr_max_boson= 10, &!nrbx! maximal nr-quantum number for bosons
                      nfx = 242, & ! maximal dimension F of one k-block
                      ngx = 264, & ! maximal dimension G of one k-block
                      nhx = nfx  + ngx, &
                      nfgx = nfx * ngx

! fixed text
character(len=1) :: tp(2) = ['+','-'], &
                    tis(2) = ['n','p'], &
                    tl(0:20) = ['s','p','d','f','g','h','i','j','k','l','m','n','o','P','q','r','S','t','u','v','w']
character(len=8) :: tit(2) = ['Neutron:','Proton: ']




!!!Option
integer :: icou = 1, & ! Coulomb-field: not at all (0), direct term (1), plus exchange (2), ususally set to 1.
           icm  = 0, & ! Center of mass:  not the usual sense of center of mass correction, make sure icm=0 is used
           itx  = 2    ! 1: only neutrons;  2: neutrons and protons                 
END MODULE Constants
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
