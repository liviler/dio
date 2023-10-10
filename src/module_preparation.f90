!===============================================================================!
! MODULE Preparation                                                            !
!   In this module, we call some subroutine(from other module or in this module)!
!to set or caculate                                                                     !
!==============================================================================!
MODULE Preparation
use Constants, only:r64,hbc,radius_r0,pi,zero,two,third,ngh,ngl,nghl,pi,u_config, icou, icm
use Globals , only: HO,force,gauss,pairing,constraint,nucleus_attributes,outputfile,iteration
use Inout, only: read_file_b23,read_file_dio
use Forces, only : set_force_parameters
use Nucleus, only: set_nucleus_attributes
use Basis, only : set_HO_basis_parameters
use MathMethods, only: GaussLaguerre, GaussHermite, math_gfv
implicit none

contains


subroutine do_preparation(ifPrint)
    logical,intent(in),optional :: ifPrint
    call read_file_b23
    call read_file_dio(ifPrint .and. .True.)
    call set_force_parameters(ifPrint .and. .True.)
    call set_nucleus_attributes(ifPrint .and. .True.)  ! set_nucleus_attributes should after set_force_parameters
    call set_HO_basis_parameters(ifPrint .and. .True.) ! set_HO_basis_parameters should after set_nucleus_attributes 
    call set_gauss_parameters(ifPrint .and. .True.)    ! set_gauss_parameters should after set_HO_basis_parameters
    call set_pairing_parameters(ifPrint .and. .True.)  ! set_pairing_parameters should after set_nucleus_attributes and set_HO_basis_parameters
    call set_constraint_parameters                     ! set_constraint_parameters should after set_nucleus_attributes and set_gauss_parameters
    call math_gfv

    call printOthersInPreparation(ifPrint .and. .True.)
    contains 
    subroutine printOthersInPreparation(ifPrintSub)
        logical,intent(in),optional :: ifPrintSub
        if(.not.ifPrintSub) return ! end this subroutine if not print it.
        open(u_config, file=outputfile%config, status='unknown', position='append')
        write(u_config,*) '*************************BEGIN print_others_in_preparation ********************'
        if(icou == 0) then
            write(u_config,*) 'Without Coulomb force'
        else if (icou == 1) then
            write(u_config,*) 'With Coulomb force'
        else if (icou == 2) then
            write(u_config,*) 'With Coulomb force with exchange'
        else 
            write(u_config,*) '[Preparation]: Wrong icou !'
        endif
        if(icm == 1) then
            write(u_config,*)'With center of mass correction'
        endif
        write(u_config,*)'Mixing- Parameter xmix: ',iteration%xmix
        write(u_config,"(a,/)") '*************************END print_others_in_preparation ********************'
        close(u_config)
    end subroutine printOthersInPreparation
end subroutine do_preparation

subroutine set_gauss_parameters(ifPrint)
    logical, intent(in), optional :: ifPrint
    integer :: ih,il
    gauss%nh = ngh
    gauss%nl = ngl
    ! set causs integral parameters
    call GaussHermite(gauss%nh,gauss%xh,gauss%w_hermite)
    call GaussLaguerre(gauss%nl,gauss%xl,gauss%w_laguerre)
    !set z-axis for both - and +
    do ih=1,gauss%nh
        gauss%zb(ih) = gauss%xh(ih)*HO%b0*HO%bz
        gauss%wh(ih) = gauss%w_hermite(ih)*exp(gauss%xh(ih)**2)
    enddo
    !set r-axis
    do il=1,gauss%nl
        gauss%rb(il) = sqrt(gauss%xl(il))*HO%b0*HO%bp
        gauss%wl(il) = gauss%w_laguerre(il)*exp(gauss%xl(il))
        do ih=1,gauss%nh
            gauss%wdcor(ih,il) = HO%b0**3 * HO%bz*HO%bp**2 * pi * gauss%wh(ih) * gauss%wl(il)
        enddo
    enddo

    if(ifPrint) call printGaussParameters
    contains
    subroutine printGaussParameters
        character(len=*), parameter :: format1= "(a30,35f15.10)"
        open(u_config, file=outputfile%config, status='unknown', position='append')
        write(u_config,*) '*************************BEGIN set_gauss_parameters********************'
        write(u_config,*) 'Gauss-Hermite:  b_z = ',HO%b0*HO%bz
        write(u_config,format1) 'points(zeta)               :',gauss%xh
        write(u_config,format1) 'weight(w_h)                :',gauss%w_hermite
        write(u_config,format1) 'zb(z=zeta*b_z)             :',gauss%zb
        write(u_config,format1) 'wh(w_h*e^{zeta^2})         :',gauss%wh
        write(u_config,*) 'Gauss-Lagueere: b_prep = ',HO%b0*HO%bp
        write(u_config,format1) 'points(eta)                :',gauss%xl
        write(u_config,format1) 'weight(w_l)                :',gauss%w_laguerre
        write(u_config,format1) 'rb(r_prep=sqrt(eta)*b_prep):',gauss%rb
        write(u_config,format1) 'wl(w_l*e^{eta})            :',gauss%wh    
        write(u_config,"(a,/)") '*************************END set_gauss_parameters********************'
    end subroutine printGaussParameters
end subroutine set_gauss_parameters

subroutine set_pairing_parameters(ifPrint)
    logical,intent(in),optional :: ifPrint
    integer :: it,ii
    ! pairing force
    do it = 1,2
        pairing%gg(it) = pairing%ga(it)/nucleus_attributes%mass_number + 1.d-10
        if (pairing%ide(it) == 4) then
            do ii =1,nghl
                pairing%delq(ii,it) = pairing%del(it)
            enddo
        endif
    enddo

    if(ifPrint) call printPairingProperties
    contains
    subroutine printPairingProperties
        character(len=*), parameter :: format1 = "(a20,5x,2f10.6)", &
                                       format2 = "(a,2(f8.4,2h/A),2f10.6,/)"
        open(u_config, file=outputfile%config, status='unknown', position='append')
        write(u_config,*) '*************************BEGIN set_pairing_parameters ********************'
        write(u_config,format1) 'Gap parameter(dec) :', pairing%dec
        write(u_config,format1) 'Gap parameter(del) :', pairing%del
        write(u_config,format1) 'Pairing Window     :', (pairing%pwi+7.d0)/HO%hom,pairing%pwi
        write(u_config,format2) 'Pairing const.     :', pairing%ga,pairing%gg
        write(u_config,"(a,/)") '*************************END set_pairing_parameters ********************'
        close(u_config)
    end subroutine printPairingProperties
end subroutine set_pairing_parameters

subroutine set_constraint_parameters
    integer :: ih,il,ii
    real(r64) :: R00,tmp_fac,tmp_c1,tmp_c2,tmp_c3
    !constraining field
    constraint%c1x = zero
    constraint%c2x = zero
    constraint%c3x = zero
    R00 = radius_r0 * nucleus_attributes%mass_number**third
    tmp_fac = 4*pi/3/nucleus_attributes%mass_number
    tmp_c1 = tmp_fac*dsqrt(3/(4*pi)) /R00
    tmp_c2 = tmp_fac*dsqrt(5/(16*pi))/(R00**2)
    tmp_c3 = tmp_fac*dsqrt(7/(16*pi))/(R00**3)

    constraint%cspr = constraint%cspr * nucleus_attributes%mass_number
    constraint%cmax = constraint%cmax * nucleus_attributes%mass_number

    do ii =1,nghl
        il = 1 + (ii-1)/ngh
        ih = ii - (il-1)*ngh
        constraint%vc(ii,1) = tmp_c1 * gauss%zb(ih)
        constraint%vc(ii,2) = tmp_c2 * (two*gauss%zb(ih)**2 - gauss%rb(il)**2)
        constraint%vc(ii,3) = tmp_c3 * (two*gauss%zb(ih)**3 - 3.d0*gauss%zb(ih)*gauss%rb(il)**2)
    enddo
end subroutine set_constraint_parameters

END MODULE Preparation