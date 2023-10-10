!==============================================================================!
! MODULE Potential                                                             !
!                                                                              !
! This module  initializes potentials                                          !                                                                      !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE Potential
use Constants, only: r64, third,zero,one,pi,half,icou,hbc
use Globals, only: inin, woodssaxon,nucleus_attributes,gauss,constraint,&
                   potential
use Inout, only: read_fields
implicit none

contains


!-------------------------------------------------------------------------------------!
! set_potential                                                                       !
!Purpose : initializes the potential. init=0: reads fields from a file that already   ! 
!          exists. init=1: calculates fields in saxon-woods form                      !                               !
!                                                                                     !
!-------------------------------------------------------------------------------------!
subroutine set_potential
    if (inin==0) then 
        call read_fields
    else if(inin==1) then
        call calculate_fields
    else
        stop '[Potential]: inin must be 0 or 1'
    endif 
end subroutine set_potential

subroutine calculate_fields
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
        vp(it)  = woodssaxon%v0 * (one - woodssaxon%akv *(-1)**it * ( nucleus_attributes%proton_number - nucleus_attributes%neutron_number)/nucleus_attributes%mass_number)
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
                 potential%vps(ih,il,it) = u
                 potential%vms(ih,il,it) = w 
                 potential%vpstot(ih,il,it) = u
                 potential%vmstot(ih,il,it) = w
                 pairing%delq(ih,il,it) = del(it)

                 ! Coulomb, approximate the coulomb potential
                 fields%coulomb(ih,il) = zero
                 if(icou /= 0 .and. it==2) then ! only for proton
                    if(r < rrv(2)) then
                        tmp_c = half * (3/rrv(2)-rr/(rrv(2)**3))
                    else 
                        tmp_c = one/r
                    endif
                    fields%coulomb(ih,il) = tmp_c * nucleus_attributes%proton_number
                    potential%vps(ih,il,2) = potential%vps(ih,il,2) + hbc*fields%coulomb(ih,il)
                    potential%vms(ih,il,2) = potential%vms(ih,il,2) + hbc*fields%coulomb(ih,il)
                    potential%vpstot(ih,il,2) = potential%vpstot(ih,il,2) + hbc*fields%coulomb(ih,il)
                    potential%vmstot(ih,il,2) = potential%vmstot(ih,il,2) + hbc*fields%coulomb(ih,il)
                 endif
            enddo
        enddo
    enddo

end subroutine calculate_fields

END MODULE Potential