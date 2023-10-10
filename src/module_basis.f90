!==============================================================================!
! MODULE Basis                                                                 !
!                                                                              !
! This module contains the variables and routines related to the Harmonic Osc- !
! illator model space.                                                         !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_basis                                                       !
!==============================================================================!
MODULE Basis
use Constants, only: zero, two, third, hbc,&
                     nb_max,tp,nhx,nt_max,nz_max,nr_max,ml_max,ntb_max,&
                     nr_max_boson, nz_max_boson, u_config
use Globals, only: HO, force, nucleus_attributes,outputfile
implicit none

contains
subroutine set_HO_basis_parameters(ifPrint)
    !------------------------------------------------------------------------------!
    ! subroutine set_HO_basis                                                         !
    !                                                                              !
    ! Set number of levels, quantum numbers, oscillator parameters, ...            !
    !                                                                              !
    !                                                                              !
    !------------------------------------------------------------------------------!
    logical,intent(in),optional :: ifPrint
    ! basis parameters
    HO%hb0 = hbc/(two * force%masses%amu)
    HO%hom = 41.0* nucleus_attributes%mass_number**(-third)
    if(HO%b0 <= zero) then
        HO%b0 = sqrt(two * HO%hb0 /HO%hom)
    endif
    HO%bp = HO%q ** ((-1.d0/6.d0))
    HO%bz = HO%q ** ((1.d0/3.d0))
    call set_HO_basis_Fermions
    call set_HO_basis_Bosons
    
    if(ifPrint) call printBasisParameters
    contains
    subroutine printBasisParameters
        character(len=*),parameter :: format1="(2(a16,i5))", &
                                      format2="(a16,f12.5)"
        open(u_config, file=outputfile%config, status='unknown', position='append')
        write(u_config,*)'*************************BEGIN set_HO_basis_parameters********************'
        write(u_config,format1) 'n0f    =',HO%n0f, 'n0b = ', HO%n0b
        write(u_config,format2) 'hb0    =',HO%hb0
        write(u_config,format2) 'hom    =',HO%hom
        write(u_config,format2) 'b0     =',HO%b0
        write(u_config,format2) 'beta0  =',HO%beta0
        write(u_config,format2) 'b_z/b0(bz) =',HO%bz
        write(u_config,format2) 'b_z    =',HO%b0*HO%bz
        write(u_config,format2) 'b_prep/b0(bp) =',HO%bz
        write(u_config,format2) 'b_prep =',HO%b0*HO%bp
        write(u_config,*)'The quantum number of oscillator base has been caculated!'
        write(u_config,"(a,/)")'*************************END set_HO_basis_parameters********************'
        close(u_config)
    end subroutine printBasisParameters

end subroutine set_HO_basis_parameters

subroutine set_HO_basis_Fermions
    !------------------------------------------------------------------------------!
    ! Oscillator-Base for Fermions                                                 !
    !------------------------------------------------------------------------------!
    integer :: il,k,ib,ip,ilk,im,ir,iz,nn,ipf,ipg,i1,i2,nz1,nr1,ml1,nz2,nr2,ml2,&
               nn1,nza2,nze2,ipa,iba,ie,i3
    integer,dimension(nhx,2) :: nzz,nrr,mll
    integer,dimension(2) :: idd,idds
    logical :: loz,lor,lom 
    HO%nzm = 0
    HO%nrm = 0
    HO%mlm = 0
    il = 0 ! all level numbers
    do k = 1,nb_max !loop over parity  
        !-----------------------------------------------------------------------
        ! Basis for the large components f  
        ! 2*n_r + n_z +  m_l  should not be larger than N_F( HO%n0f )                        
        ! k=Omega+1/2, ip=1 for parity '+', ip=2 for parity '-'
        ! Omega = m_l + m_s
        ! ib is an index to order different compinations of Omega and parity                                              
        !-----------------------------------------------------------------------
        ib = k 
        ip = 1
        HO%ikb(ib) = k 
        HO%ipb(ib) = ip
        write(HO%txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  K = ',k+k-1,'/2',tp(ip)
        ilk = 0
        ! the following three loops make sure that the quantum numbers satisfy nn=2*ir+im+iz 
        ! and pi=(-1)^{n} and k=Omega+1/2=m or m+1
        do im = k-1,k ! loop over quantum number ml
            do ir=0, (HO%n0f -im)/2 ! loop over quantum number nr
                do iz=0, HO%n0f  ! loop over quantum number nz
                    nn = 2*ir + im +iz
                    if(nn > HO%n0f) exit 
                    ilk = ilk + 1
                    if(ilk >nhx) stop '[Basis]: nhx too small'
                    nzz(ilk,ip) = iz
                    nrr(ilk,ip) = ir
                    mll(ilk,ip) = im
                enddo
            enddo
        enddo
        idd(ip) = ilk
        !----------------------------------------------------------------------
        ! Basis for the small components g                            
        ! N_F(for small components) should be larger than N_F(for larger components) by 1, 
        ! the extra states are called spurious states                                            
        !------------------------------------------------
        ipf = ip
        ipg = ip
        ilk = idd(ipg)
        do i1 = 1,idd(ipf)
            nz1 = nzz(i1,ipf)
            nr1 = nrr(i1,ipf)
            ml1 = mll(i1,ipf)

            loz = .true.
            lor = .true.
            lom = .true.

            do i2 = 1, ilk
                nz2 = nzz(i2,ipg)
                nr2 = nrr(i2,ipg)
                ml2 = mll(i2,ipg)
                !nz1+1
                if(loz .and. nz2.eq.nz1+1 .and. nr2.eq.nr1 .and. ml2.eq.ml1) loz = .false.
                ! nr1+1, ml1-1
                if(lor .and. nr2.eq.nr1+1 .and. nz2.eq.nz1 .and. ml2.eq.ml1-1) lor = .false.
                !ml1+1
                if(lom .and. ml2.eq.ml1+1 .and. nr2.eq.nr1 .and. nz2.eq.nz1) lom = .false. 
            enddo
            !add them to the basis of small components
            !the following three 'if' statements check whether the above mentioned relations among the quantum numbers are still satisfied
            if(loz) then
                nn1 = nz1 + 2*nr1 + ml1 + 1
                nza2= mod(nn1-ml1,2)/2                                
                nze2= (nn1-ml1)/2 
                if (nza2.le.(nz1+1)/2.and.(nz1+1)/2.le.nze2) then ! always true
                    ilk = ilk + 1                                      
                    nzz(ilk,ipg) = nz1 + 1                             
                    nrr(ilk,ipg) = nr1                                 
                    mll(ilk,ipg) = ml1
                endif
            endif
            if (lor) then                                            
                 nn1= nz1 + 2*nr1 + ml1 +1                             
                 if (k-1.le.ml1-1.and.ml1-1.le.min(k,nn1)) then  ! ml1==k   
                     ilk = ilk + 1                                      
                     nzz(ilk,ipg) = nz1                                 
                     nrr(ilk,ipg) = nr1 + 1                             
                     mll(ilk,ipg) = ml1 - 1                             
                 endif                                                 
            endif
            if (lom) then                                            
                nn1= nz1 + 2*nr1 + ml1 +1                             
                if (k-1.le.ml1+1.and.ml1+1.le.min(k,nn1)) then  !ml1==k−1      
                    ilk = ilk + 1                                      
                    nzz(ilk,ipg) = nz1                                 
                    nrr(ilk,ipg) = nr1                                 
                    mll(ilk,ipg) = ml1 + 1    
                endif                                                 
            endif
            
            if(ilk > nhx) stop '[Basis]: nhx too small'
        enddo
        idds(ipg)= ilk - idd(ipg)
        !--------------------------------------------------------------------------------------------
        ! Reordering and construction of the fields ia, id, iag0
        ! Attention:  note that through ipa=3-ip; ib=2*(k-1)+ip, iba=2*(k-1)+ipa ???
        !             and id(ib,2)=idd(ipa)+idds(ipa)
        !             and idg0(ib)=idds(ipa)
        !             and ia(iba,2)=il
        !             and iag0(iba)=il+idd(ip)
        !             the index has been redirected, so the ordering described at the head is ahieved.
        !          id(ib,1) ia(ib,1) for large component f 
        !          id(ib,2) ia(ib,2) for small component g
        !--------------------------------------------------------------------------------------------
        ipa = ip
        iba = ib
        HO%id(ib,1) = idd(ip)
        HO%id(ib,2) = idd(ipa) + idds(ipa)
        HO%idg0(ib) = idds(ipa)

        HO%ia(ib,1) = il
        HO%ia(iba,2) = il
        HO%iag0(iba) = il + idd(ip)

        ie = idd(ip) + idds(ip)
        if(il+ie > nt_max) stop '[Basis]: nt_max too small'
        do i3=1,ie
            il = il + 1
            HO%nz(il) = nzz(i3,ip)
            HO%nr(il) = nrr(i3,ip)
            HO%ml(il) = mll(i3,ip)                                       
            HO%ms(il) = (-1)**(k-HO%ml(il)+1) !ms=1: spin up Omega=m+1/2, ms=-1: spin down, Omega=m-1/2
            nn = HO%nz(il) +2*HO%nr(il) + HO%ml(il)
            if(nn < 10) then
                write(HO%tb(il),"(i2,a1,'[',3i1,']')") 2*k-1,tp(ip),nn,HO%nz(il),HO%ml(il)
            else
                write(HO%tb(il),"(4i2)") 2*k-1,nn,HO%nz(il),HO%ml(il)
            endif
            HO%nzm = max0(HO%nzm,HO%nz(il))
            HO%nrm = max0(HO%nrm,HO%nr(il))
            HO%mlm = max0(HO%mlm,HO%ml(il))
        enddo
    enddo 
    HO%nt = il
    HO%nb = ib
    if(HO%mlm > ml_max) stop '[Basis]: ml_max too small'
    if(HO%nrm > nr_max) stop '[Basis]: nr_max too small'
    if(HO%nzm > nz_max) stop '[Basis]: nz_max too small'

end subroutine set_HO_basis_fermions

subroutine set_HO_basis_Bosons
    !------------------------------------------------------------------------------!
    ! Oscillator-Base for Bosons                                                   !
    !------------------------------------------------------------------------------!
    
    ! if (mod(HO%n0b,2).ne.0) stop ' [Basis]: n0b must be even'
    integer :: il,iz,ir
    HO%nzbm = 0
    HO%nrbm = 0
    il = 0
    do iz =0,HO%n0b !loop over nz-quantum number
        do ir=0, (HO%n0b-iz)/2 !loop over nz-quantum number
            if (iz + 2*ir > HO%n0b) exit
            il = il + 1
            if(il > ntb_max) stop '[Basis]: ntb_max too small'
            HO%nrb = ir
            HO%nzb = iz
            HO%nrbm = max0(HO%nrbm,ir)
            HO%nzbm = max0(HO%nzbm,iz)
        enddo
    enddo
    HO%ntb = il
    if(HO%nrbm > nr_max_boson) stop '[Basis]: nr_max_boson too small'
    if(HO%nzbm > nz_max_boson) stop '[Basis]: nz_max_boson too small'
end subroutine set_HO_basis_Bosons

END MODULE Basis