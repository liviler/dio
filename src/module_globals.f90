!==============================================================================!
! MODULE Globals                                                               !
!                                                                              !
! This module defines global variables that are shared by subroutines.         !                                                  !
!==============================================================================!
MODULE Globals

use Constants, only: r32,r64,i8,i16,i32,igfv,igfvbc,ngh,ngl,nghl,&
                     nt_max,ntb_max,nb_max,nz_max,nr_max,ml_max,nz_max_boson,nr_max_boson,&
                     nfgx
implicit none

public

!-------------------------Initialization---------------------------------------
type Constraint_Parameter
    integer(i32) :: length ! length of data store in betac,bet3c,clam2
    integer(i32) :: index ! index of betac,bet3c,clam2
    real(r64), dimension(:), allocatable :: betac ! Constraint beta2-values 
    real(r64), dimension(:), allocatable :: bet3c ! Constraint beta3-values
    real(r64), dimension(:), allocatable :: clam2
    ! quadratic constraint
    integer(i8) :: icstr ! quadratic constraint (no 0; beta2 1; b2+b3 2)
    real(r64) :: cspr ! spring constant
    real(r64) :: cmax ! cutoff for dE/db

    real(r64) :: c1x
    real(r64) :: c2x
    real(r64) :: c3x
    real(r64) :: c1xold
    real(r64) :: c2xold
    real(r64) :: c3xold

    real(r64) :: calq1
    real(r64) :: calq2
    real(r64) :: calq3

    real(r64), dimension(nghl,3) :: vc

end type
type(Constraint_Parameter) :: constraint

type Harmonic_Oscillator_Basis !!!!
    integer(i8) :: n0f,n0b  !n0f n0b  ! number of oscillator shells
    real(r64)   :: b0       ! b_0 oscillator parameter(fm**-1) of basis
    real(r64)   :: beta0    !deformation parameter of basis
    real(r64)   :: q        !momentum?
    real(r64)   :: hb0      ! hbc/(2* amu)
    real(r64)   :: hom      ! hbar*omega_0
    real(r64)   :: bp       ! b_p/b_0
    real(r64)   :: bz       ! b_z/b_0
    
    integer,dimension(nt_max) :: nz,nr,ml,ms
    character(len=8),dimension(nt_max) :: tb
    integer :: nzm,nrm,mlm ! max of nz,nr,ml
    integer,dimension(nb_max,2) :: ia  ! ia(b,1):begin of the large components of block b is ia(b,1)+1 ;ia(b,2):begin of the small components of block b is ia(b,2)+1
    integer,dimension(nb_max,2) :: id  ! id(b,1):dimension large components of block b;id(b,2):dimension small components of block b
    integer,dimension(nb_max) :: iag0  ! begin of the spurious components of block b is iag0(b)+1 
    integer,dimension(nb_max) :: idg0  ! number of spurious components of block b
    integer,dimension(nb_max) :: ikb   ! K-quantum number of each block (K+1/2) 
    integer,dimension(nb_max) :: ipb   ! parity of each block ( 1 for +, 2 for -)
    character(len=25),dimension(nb_max) :: txb
    integer :: nb ! number of K-parity-blocks 
    integer :: nt ! number of total levels

    integer,dimension(ntb_max) :: nzb,nrb
    integer :: nzbm,nrbm ! max of nzb,nrb
    integer :: ntb !no ! number of total levels for bosons

end type
type(Harmonic_Oscillator_Basis) :: HO ! harmonic oscillator basis


type WoodsSaxon_Parameters
    real(r64) :: qs !qs = exp(beta2 * 3*sqrt(5/(16*pi)))
    real(r64) :: beta2!betas ! deformation beta2 of WoodsSaxon potential, beta2 = ln(qs)/(3*sqrt(5/(16*pi)))
    real(r64) :: beta3!bet3s ! deformation beta3 of WoodsSaxon potential
    
    real(r64) :: v0
    real(r64) :: akv
    real(r64), dimension(2) :: r0v
    real(r64), dimension(2) :: av
    real(r64), dimension(2) :: vso
    real(r64), dimension(2) :: rso
    real(r64), dimension(2) :: aso
end type
type(WoodsSaxon_Parameters) :: woodssaxon



type Fields_
    real(r64),dimension(ngh,ngl,2) :: vps !V+S
    real(r64),dimension(ngh,ngl,2) :: vms !V-S
    real(r64),dimension(ngh,ngl,2) :: vpstot
    real(r64),dimension(ngh,ngl,2) :: vmstot
 
    real(r64),dimension(ngh,ngl) :: coulomb
    real(r64),dimension(ngh,ngl) :: sigma
    real(r64),dimension(ngh,ngl) :: omega
    real(r64),dimension(ngh,ngl) :: rho
end type
type(Fields_) :: fields


type Iteration_Parameters 
    integer(r32) :: iteration_max !maxi ! max number of iterations
    real(r32)    :: xmix  ! maxing parameter
end type
type(Iteration_Parameters) :: iteration

type Pairing_Parameters
    integer(i8), dimension(2) :: ide
    real(r32), dimension(2) :: dec !dec !frozen gaps parameter(neutrons and protons)
    real(r32), dimension(2) :: ga ! Pairing-Constants GG = GA/A
    real(r32), dimension(2) :: del !Initial values for gap-parameters
    real(r64), dimension(2) :: gg
    real(r64), dimension(2) :: spk
    real(r64) :: pwi
    real(r64), dimension(2) :: vpair !pairing strength for delta force
    real(r64), dimension(nghl,2) :: delq
end type
type(Pairing_Parameters) :: pairing

type Nucleus_Parameters
    character(len=3) :: name = ''!nucnam ! nucleus name
    integer(i16) :: mass_number_int!nama ! mass number
    real(r64) :: mass_number !amas ! float type of mass_number
    real(r64) :: proton_number = 0!npro ! proton number
    real(r64) :: neutron_number!nneu ! neutron number
end type
type(Nucleus_Parameters) :: nucleus_attributes

type Output_FileName
    character(len=28) :: config = './output/dio_config.out'
    character(len=28) :: outputf
    character(len=28) :: outputw
    character(len=28) :: outdel
    character(len=28) :: outputwf
end type
type(Output_FileName) :: outputfile

! -----define force parameters
type mass_parameters
    real(r64) :: amu
    real(r64) :: amsig
    real(r64) :: amome
    real(r64) :: amdel
    real(r64) :: amrho
    real(r64) :: ampi
end type
type couplg_parameters
    ! quadratic terms
    real(r64) :: ggsig
    real(r64) :: ggome
    real(r64) :: ggdel
    real(r64) :: ggrho
    integer(i8), dimension(4) :: lmes
end type
type coupld_parameters
    ! derivative terms
    real(r64) :: ddsig
    real(r64) :: ddome
    real(r64) :: dddel
    real(r64) :: ddrho
end type
type coupnl_parameters
    ! non-linear terms
    real(r64) :: ggbet
    real(r64) :: gggams
    real(r64) :: gggamv
end type
type couplm_parameters
    real(r64) :: gsig
    real(r64) :: gome
    real(r64) :: gdel
    real(r64) :: grho
    real(r64) :: gpi
end type
type nonlin_parameters
    real(r64) :: g2
    real(r64) :: g3
    real(r64) :: c3
end type
type dtypel_parameters
    real(r64) :: a_s
    real(r64) :: b_s
    real(r64) :: c_s
    real(r64) :: d_s
    real(r64) :: a_v
    real(r64) :: b_v
    real(r64) :: c_v
    real(r64) :: d_v
    real(r64) :: a_tv
    real(r64) :: b_tv
    real(r64) :: c_tv
    real(r64) :: d_tv
    real(r64) :: dsat
end type
type dpolyn_parameters
    real(r64), dimension(4) :: alf
    real(r64), dimension(4) :: bet
    real(r64), dimension(4) :: gam
end type
type option_parameters
    integer(i8) :: inl !non-linear meson-coupling
    integer(i8) :: ipc !density-dependent meson-coupling
    integer(i8) :: idd !point-coupling
end type

type Force_Parameter_Sets
    character(len=10) :: parname ! Parameterset name of the Lagrangian
    type(mass_parameters) :: masses
    type(couplg_parameters) :: couplg
    type(coupld_parameters) :: coupld
    type(coupnl_parameters) :: coupnl
    type(couplm_parameters) :: couplm
    type(nonlin_parameters) :: nonlin
    type(dtypel_parameters) :: dtypel
    type(dpolyn_parameters) :: dpolyn
    type(option_parameters) :: option
    real(r64) :: rosat
end type
type(Force_Parameter_Sets) :: force
! -----end define force parameters

type math_gfv
    integer,dimension(0:igfv) :: iv
    real(r64),dimension(0:igfv) :: sq
    real(r64),dimension(0:igfv) :: sqi
    real(r64),dimension(0:igfv) :: sqh
    real(r64),dimension(0:igfv) :: shi
    real(r64),dimension(0:igfvbc,0:igfvbc) :: ibc
    real(r64),dimension(0:igfv) :: fak
    real(r64),dimension(0:igfv) :: fad
    real(r64),dimension(0:igfv) :: fi
    real(r64),dimension(0:igfv) :: fdi
    real(r64),dimension(0:igfv) :: wf
    real(r64),dimension(0:igfv) :: wfi
    real(r64),dimension(0:igfv) :: wfd
    real(r64),dimension(0:igfv) :: gm2
    real(r64),dimension(0:igfv) :: gmi
    real(r64),dimension(0:igfv) :: wg
    real(r64),dimension(0:igfv) :: wgi
end type
type(math_gfv) :: gfv

type Gauss_Integral
    ! GaussHermite
    integer :: nh ! points length
    real(r64),dimension(ngh) :: xh ! store guass points, actually is zeta
    real(r64),dimension(ngh) :: w_hermite ! store the weights of the related points
    real(r64),dimension(ngh) :: zb ! zb = xh * b_z, actually is z
    real(r64),dimension(ngh) :: wh ! wh = w_hermite * exp(xh**2)
    !GaussLaguerre
    integer :: nl
    real(r64),dimension(ngl) :: xl ! actually is eta
    real(r64),dimension(ngl) :: w_laguerre
    real(r64),dimension(ngl) :: rb ! rb = sqrt(xl) * b_\prep, actually is r_prep
    real(r64),dimension(ngl) :: wl ! wl = w_laguerre * exp(xl)

    real(r64),dimension(ngh,ngl) :: wdcor
end type
type(Gauss_Integral) :: gauss

type Wavefunction_HO
    real(r64),dimension(:,:),allocatable :: qh  ! z-function at meshpoint xh(ih) with sqrt(wh)
    real(r64),dimension(:,:),allocatable :: qh1 ! the first-order derivative of z-funtion with sqrt(wh)
    real(r64),dimension(:,:,:),allocatable :: ql ! r-function at meshpoint xl(il) with sqrt(wl/2) 
    real(r64),dimension(:,:,:),allocatable :: ql1 ! the first-order derivative of r-funtion with r*sqrt(wl/2);NOTE: r d/dr psi_(nr,L)(rr) = QL1(nr,L,i) * sqrt( 2 / WL(i) )
    real(r64),dimension(:,:),allocatable :: qhb ! the z-function for the expansion of the mesonfields
    real(r64),dimension(:,:),allocatable :: qlb  ! the r-function for the expansion of the mesonfields
end type

type Wavefunction
    type(Wavefunction_HO)  HO
end type
type(Wavefunction) WF 

type Matrix_
    real(r64),dimension(nfgx,nb_max) :: sp
    real(r64),dimension(ntb_max*ntb_max,4) :: meson_propagators 
end type
type(Matrix_) matrices

integer(i8) :: inin !Initialization of wavefunctions. 1(calc. from beginning);0(read saved pot.)

type erwar_ ! what's the meaning of erwar?
    real(r64) :: ea
    real(r64) :: rms
    real(r64) :: qp
end type
type(erwar_) :: erwar

type fermi_
    real(r64) :: ala
    real(r64) :: tz
end type
type(fermi_) :: fermi








END MODULE Globals