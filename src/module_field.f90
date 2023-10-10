!====================================================================================!
! MODULE Field                                                                      !
!                                                                                    !
! This module contains subroutines related to  
!                                                                                    !
!                                                                                    !
! List of subroutines and functions:                                                 !
! - subroutine                                                                       !
!====================================================================================!
MODULE Field
use Constants, only: r64,zero,one,third,half,pi,icou,hbc,u_config
use Globals, only: inin,woodssaxon,gauss,constraint,nucleus_attributes,fields,&
                   outputfile
use Inout, only: read_fields
implicit none    
 

contains

subroutine generate_initial_fields(ifPrint)
    !-------------------------------------------------------------------------------------!
    ! given_initial_field                                                                       !
    !Purpose : initializes the potential.                                                 !
    !          init=0: reads fields from a file that already exists.                      !
    !          init=1: calculates fields in saxon-woods form.                             !
    !                                                                                     !
    !-------------------------------------------------------------------------------------!
    logical, intent(in), optional :: ifPrint
    if (inin==0) then 
        ! read from input files
        call read_fields
    else if(inin==1) then
        ! saxon-woods fields
        call calculate_saxonwoods_fields(ifPrint .and. .True.)
    else
        stop '[Potential]: inin must be 0 or 1'
    endif 
end subroutine generate_initial_fields

subroutine calculate_saxonwoods_fields(ifPrint)
    logical, intent(in), optional :: ifPrint
    integer :: it,ih,il
    real(r64) :: z,zz,r,rr,cos_theta,Y20,Y30,facb, u,w,tmp_c
    real(r64),dimension(2) :: rrv, rls, vp, vls 

    ! Saxon-Woods parameter von Koepf und Ring, Z.Phys. (1991)
    ! W. Koepf and P. Ring Z. Phys. A 339, 81 (1990)
    ! 1: neutron and 2: proton
    woodssaxon%v0  = -71.28d0
    woodssaxon%akv = 0.4616d0
    woodssaxon%r0v = [1.2334d0,1.2496d0]
    woodssaxon%av  = [0.615d0,0.6124d0]
    woodssaxon%vso = [11.1175d0,8.9698d0]
    woodssaxon%rso = [1.1443d0,1.1401d0]
    woodssaxon%aso = [0.6476d0,0.6469d0]

    do it =1,2
        rrv(it) = woodssaxon%r0v(it)*nucleus_attributes%mass_number**third 
        rls(it) = woodssaxon%rso(it)*nucleus_attributes%mass_number**third 
        vp(it)  = woodssaxon%v0 * (one - woodssaxon%akv &
                *(-1)**it * ( nucleus_attributes%proton_number - nucleus_attributes%neutron_number) &
                /nucleus_attributes%mass_number)
        vls(it) = vp(it) * woodssaxon%vso(it)
        do ih =1,gauss%nh
            z  = gauss%zb(ih)
            zz = z**2
            do il = 1,gauss%nl
                rr = (gauss%rb(il)**2 + zz)
                r  = sqrt(rr)
                cos_theta = z/r ! cos(\theta)
                Y20 = sqrt(5/(16*pi))*(3*cos_theta**2-one)
                Y30 = sqrt(7/(16*pi))*(5*cos_theta**3-3*cos_theta)
                if(constraint%icstr == 0) then
                    facb = one + 0.75*(woodssaxon%beta2*Y20 + woodssaxon%beta3*Y30)
                elseif(constraint%icstr == 1) then
                    facb = one + 0.75*(constraint%betac(constraint%index)*Y20 + woodssaxon%beta3*Y30)
                elseif(constraint%icstr == 2) then
                    facb = one + 0.75*(constraint%betac(constraint%index)*Y20 + constraint%bet3c(constraint%index)*Y30)
                else 
                    stop "[Potential]: constraint%icstr must be 0,1,2"
                endif
                 u= vp(it)/(one+exp((r-rrv(it)*facb)/woodssaxon%av(it))) 
                 w= -vls(it)/(one+exp((r-rls(it)*facb)/woodssaxon%aso(it)))
                
                 ! set fields
                 fields%vps(ih,il,it) = u
                 fields%vms(ih,il,it) = w 
                 fields%vpstot(ih,il,it) = u
                 fields%vmstot(ih,il,it) = w
                !  pairing%delq(ih,il,it) = pairing%del(it) ! samething have been done in moudule_preparation's set_pairing_parameters

                 ! Coulomb, approximate the coulomb potential
                 fields%coulomb(ih,il) = zero
                 if(icou /= 0 .and. it==2) then ! only for proton
                    if(r < rrv(2)) then
                        tmp_c = half * (3/rrv(2)-rr/(rrv(2)**3))
                    else 
                        tmp_c = one/r
                    endif
                    fields%coulomb(ih,il) = tmp_c * nucleus_attributes%proton_number
                    fields%vps(ih,il,2) = fields%vps(ih,il,2) + hbc*fields%coulomb(ih,il)
                    fields%vms(ih,il,2) = fields%vms(ih,il,2) + hbc*fields%coulomb(ih,il)
                    fields%vpstot(ih,il,2) = fields%vpstot(ih,il,2) + hbc*fields%coulomb(ih,il)
                    fields%vmstot(ih,il,2) = fields%vmstot(ih,il,2) + hbc*fields%coulomb(ih,il)
                 endif
            enddo
        enddo
    enddo

    if(ifPrint) call printSaxonWoodsField
    contains
    subroutine printSaxonWoodsField
        character(len=*),parameter :: format1 = "(a,2f10.4)"
        open(u_config,file=outputfile%config,status='unknown',position='append')
        write(u_config,*)  '*************************BEGIN calculate_saxonwoods_fields ********************'
        write(u_config,*) 'Calculates fields in saxon-woods form !'
        ! coulomb force
        if(icou == 0) then
            write(u_config,*) 'Without Coulomb force'
        else if (icou == 1) then
            write(u_config,*) 'With Coulomb force'
        else if (icou == 2) then
            write(u_config,*) 'With Coulomb force with exchange'
        else 
            write(u_config,*) '[Preparation]: Wrong icou !'
        endif
        ! print parameters
        write(u_config,"(/,a)") 'Saxon-Woods Parameters : '
        write(u_config,format1) ' v0     = ',woodssaxon%v0                                  
        write(u_config,format1) ' kappa  = ',woodssaxon%akv                                 
        write(u_config,format1) ' lambda = ',woodssaxon%vso                                 
        write(u_config,format1) ' r0     = ',woodssaxon%r0v                                 
        write(u_config,format1) ' a      = ',woodssaxon%av                                  
        write(u_config,format1) ' r0-so  = ',woodssaxon%rso                                 
        write(u_config,format1) ' a-so   = ',woodssaxon%aso                                 
        write(u_config,format1) ' beta2  = ',woodssaxon%beta2                              
        write(u_config,format1) ' beta3  = ',woodssaxon%beta3
        ! print fields
        write(u_config,"(/,a)") 'Caculated Fields:'
        call prigh(u_config,1,fields%vps,one,'VPS-n')
        call prigh(u_config,1,fields%vps(1,1,2),one,'VPS-p')
        call prigh(u_config,1,fields%vms,one,'VMS-n')
        call prigh(u_config,1,fields%vms,one,'VMS-p')
        write(u_config,"(a,/)")  '*************************END calculate_saxonwoods_fields *********************'
        close(u_config)
    end subroutine printSaxonWoodsField
end subroutine calculate_saxonwoods_fields
subroutine prigh(u_config,is,ff,f,text) 
    !---------------------------------------------------------------!
    ! IS = 1:  prints f*ff(z,r) at gauss-meshpoints                 !
    ! IS = 2:  prints f*ff(z,r)/wdcor at gauss-meshpoints           !
    ! (is=2) is used to output density                              !
    !  wdcor is the overall normalization factor equal to,          !
    !           wdcor=pi*WH(ih)*WL(il)* b0**3 * bz*bp*bp            !
    !---------------------------------------------------------------!
    integer :: u_config,is
    real(r64) :: f
    character(len=*) :: text
    real(r64),dimension(gauss%nh,gauss%nl) :: ff

    integer :: ixh =7, ixl = 3 ! change it to print more data; ixh should less than gauss%nh, ixl should less than gauss%nl
    integer :: ih,il
    real(r64) :: r 
    write(u_config,'(/,1x,a6,12f10.3)') text, (gauss%zb(ih),ih=1,ixh)
    do il = 1,ixl
        r = gauss%rb(il)
        if(is==1) then
            write(u_config,'(1x,f6.3,12f10.4)') r, (f*ff(ih,il),ih=1,ixh)
        elseif(is==2) then
            write(u_config,'(1x,f6.3,12f10.4)') r, (f*ff(ih,il)/gauss%wdcor(ih,il),ih=1,ixh)
        else
            write(u_config,*)'Wrong is in prigh()'
        endif
    enddo
end subroutine prigh

END MODULE Field