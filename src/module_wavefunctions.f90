!====================================================================================!
! MODULE Wavefunction                                                                !
!                                                                                    !
! This module contains subroutines related to the wave function and their densities. !
!                                                                                    !
!                                                                                    !
! List of subroutines and functions:                                                 !
! - subroutine                                                                       !
!====================================================================================!
MODULE Wavefunctions
use Constants, only: r64,pi,half,one,nz_max,nz_max_boson,ml_max,nr_max,nr_max_boson,u_config
use Globals, only: gauss,HO,WF,gfv,outputfile
implicit none

contains
subroutine choose_basis_wavefuction(ifPrint)
    logical,intent(in),optional :: ifPrint
    ! choose the cylindrical oscillator wavefunction
    call calculate_HO_wavefuction(ifPrint .and. .True.)
    
    ! choose another wavefunction as the basis

end subroutine choose_basis_wavefuction

subroutine calculate_HO_wavefuction(ifPrintHO)
    !--------------------------------------------------------------------------------!
    ! Calculates the wavefunctions for the cylindrical oscillator.      
    !     they are given as:                                                 
    !                                                                       
    !     Phi(zz,rr,phi) = 1/((b0**3 bz*bp*bp)^(1/2)) * psi_nz(zz) *        
    !                      psi_(L,nr)(rr) * exp(i*L*phi)/sqrt(2*pi)         
    !                                                                       
    !     zz is z-coordinate in units of fm                                 
    !     rr is perpendicular coordinate in units of fm                    
    !                                                                       
    !     psi_nz(zz)     = N_nz * H_nz(z) * exp(-z*z/2)                                                                                            
    !     psi_(L,nr)(rr) = N_(nr,L) * sqrt(2) * eta^(L/2) * L_nr^L(eta) * exp(-eta/2)   
    !                                                                       
    !     z = zz/(bz*b0),    r = rr/(bp*b0),       eta = r*r                                                                                      
    !     N_nz     = 1/sqrt(sqrt(pi) * 2^nz * nz!)                                                                                                 
    !     N_(nr,L) = sqrt( nr! / (nr+L)! )                                  
    !                                                                       
    !                                                                       
    !     the contribution to the density from the level i is:                                                                                     
    !     rho_k(zz,rr)= 1/(2*pi b0*3 bz*bp*bp) * ( psi_nz(zz) * psi_(L,nr)(rr) )^2                 
    !                                                                       
    !---- The z-function at meshpoint xh(ih) is stored in QH(nz,ih)         
    !     such that QH is normalized in such way that the                   
    !     norm integral reads                                               
    !                                                                       
    !     \int dz_(-inf)^(inf) (psi_nz)**2 = 1 = \sum_i QH(nz,i)**2         
    !                                                                       
    !     this means, that QH contains the following factors:               
    !                                                                       
    !     a)  the z-part of the wavefunction psi_nz(zz)                     
    !     b)  the Gaussian weight sqrt( WH(i) ):                            
    !         \inf_(-inf)^inf f(z) dz = \sum_i f(x_i) * WH(i)               
    !                                                                       
    !     having QH(nz,i) we get the z-wavefunction:                                                                                               
    !     psi_nz(zz) =  QH(nz,i) / sqrt( WH(i) )                            
    !                                                                       
    !---- The r-function at meshpoint XL(il) is stored in QL(nr,L,il)       
    !     such that QL is normalized in such way that the                   
    !     2-dimensional norm integral reads                                 
    !                                                                       
    !     \int_0^inf r dr (phi_(L,nr)**2 = 1 = \sum_i QL(nr,L,i)**2         
    !                                                                       
    !     this means, that QL contains the following factors:              
    !                                                                       
    !     a)  the part of the wavefunction psi_(nr,L)(rr)                   
    !     b)  a factor sqrt(1/2) from the transformation from r to eta      
    !     c)  the Gaussian weight sqrt( WL(i) ):                            
    !         \inf_0^inf f(eta) d(eta) = \sum_i f(XL(i)) * WL(i)            
    !                                                                       
    !     having QL(nr,L,i) we get the r-wavefunction:                                                                                             
    !     psi_(nr,L)(rr) =  QL(nr,L,i) * sqrt( 2 / WL(i) )                  
    !                                                                       
    !                                                                       
    !---- the density contribution from the level k                         
    !                                                                       
    !     rho_k(zz,rr)= (QH(nz,ih) * QL(nr,L,il))^2 /                       
    !                         (  pi * WH(ih)*WL(il)* b0**3 * bz*bp*pb)      
    !                                                                       
    !----------------------------------------------------------------------
    !                                                                       
    !     QH1 contains the z-derivatives in the following form:             
    !                                                                       
    !     d/dz psi_nz(zz) = QH1(nz,i) / sqrt( WH(i) )                       
    !                                                                      
    !     QL1 contains the r-derivatives in the following form:             
    !                                                                       
    !     r d/dr psi_(nr,L)(rr) = QL1(nr,L,i) * sqrt( 2 / WL(i) )           
    !     note this difference with the text                                
    !----------------------------------------------------------------------
    !                                                                       
    !     QHB(nz,i) is the z-function for the expansion of the mesonfields  
    !     QLB(nr,i) is the r-function for the expansion of the mesonfields  
    !                                                                       
    !     QHB(nz,i) = psi_nz(zz) / (b0b*bz)^(1/2)                          
    !     QLB(nr,i) = psi_(nr,L=0)(rr) / ( sqrt(2*pi) * b0b*bp)             
    !                                                                       
    !----------------------------------------------------------------------
    ! 
    !Output:
    ! QH (0:nz_max,  1:ngh)  : z-function at meshpoint xh(ih) with sqrt(wh)
    ! QH1(0:nz_max,  1:ngh)  : the first-order derivative of z-funtion with sqrt(wh)
    ! QL (0:nr_max, 0:ml_max, 1:ngl)  : r-function at meshpoint xl(il) with sqrt(wl/2) 
    ! QL1(0:nr_max, 0:ml_max, 1:ngl)  : the first-order derivative of r-funtion with r*sqrt(wl/2);NOTE: r d/dr psi_(nr,L)(rr) = QL1(nr,L,i) * sqrt( 2 / WL(i) )
    ! QHB(0:nz_max_boson, 1:ngh)  : the z-function for the expansion of the mesonfields
    ! QLB(0:nr_max_boson, 1:ngl)  : the r-function for the expansion of the mesonfields
    !----------------------------------------------------------------------
    logical, intent(in), optional :: ifPrintHO
    real(r64) :: w4pii,z,psiH_0,qh_0,cb,zb,psiHB_0, &
                 x,psiL_00,ql_l0,tmp0,tmp1,tmp2,tmp3,tmp4, xb,psiLB_00
    integer :: ih,n,il,l
    allocate(WF%HO%qh (0:nz_max, 1:gauss%nh))
    allocate(WF%HO%qh1(0:nz_max, 1:gauss%nh))
    allocate(WF%HO%ql (0:nr_max, 0:ml_max, 1:gauss%nl))
    allocate(WF%HO%ql1(0:nr_max, 0:ml_max, 1:gauss%nl))
    allocate(WF%HO%qhb (0:nz_max_boson, 1:gauss%nh))
    allocate(WF%HO%qlb (0:nr_max_boson, 1:gauss%nl))

    !----------------------------------------------------------------------
    ! z-dependence
    !----------------------------------------------------------------------
    w4pii = pi**(-0.25d0)
    cb = one
    do ih = 1, gauss%nh
        !------------basis for the fermions----------------
        z = gauss%xh(ih)
        psiH_0 = w4pii*exp(-half*z**2)
        qh_0 = psiH_0*sqrt(gauss%wh(ih))
        WF%HO%qh(0,ih) = qh_0
        WF%HO%qh(1,ih) = gfv%sq(2) * qh_0 * z
        WF%HO%qh1(0,ih) = - qh_0 * z
        WF%HO%qh1(1,ih) = gfv%sq(2) * qh_0 * (one -z**2)
        do n = 2, nz_max
            WF%HO%qh(n,ih) = gfv%sqi(n)*(gfv%sq(n)*z*WF%HO%qh(n-1,ih)-gfv%sq(n-1)*WF%HO%qh(n-2,ih))
            WF%HO%qh1(n,ih) = gfv%sq(n+n) * WF%HO%qh(n-1,ih) - z * WF%HO%qh(n,ih)
        enddo
        !------------basis for the bosons---------------- 
        zb = z*cb
        psiHB_0 = w4pii*exp(-half*zb**2)/sqrt(HO%b0*HO%bz)
        WF%HO%qhb(0,ih) = psiHB_0
        WF%HO%qhb(1,ih) = gfv%sq(2) * psiHB_0 * zb
        do n = 2, nz_max_boson
            WF%HO%qhb(n,ih) = gfv%sqi(n)*(gfv%sq(2)*zb*WF%HO%qhb(n-1,ih)-gfv%sq(n-1)*WF%HO%qhb(n-2,ih))
        enddo
    enddo
    !----------------------------------------------------------------------
    ! r-dependence
    ! note here x is actually eta=r*r
    !----------------------------------------------------------------------
    do il = 1, gauss%nl
        !------------basis for the fermions----------------
        x = gauss%xl(il)
        psiL_00 = gfv%sq(2) * exp(-half*x)
        do l = 0,ml_max
              ql_l0 = psiL_00*sqrt(half*gauss%wl(il)*x**l)*gfv%wfi(l)
              WF%HO%ql(0,l,il) = ql_l0
              WF%HO%ql(1,l,il) = ql_l0 *(l+1-x) * gfv%wfi(l+1)/gfv%wfi(l)
              WF%HO%ql1(0,l,il) = ql_l0 * (l-x) ! why not ql_10 * (l-x)/(sqrt(x)) ?
              WF%HO%ql1(1,l,il) = ql_l0 * (l*l+l-x*(l+l+3)+x*x) * gfv%wfi(l+1)/gfv%wfi(l)
            do n = 2,nr_max
                tmp0 = gfv%sqi(n) * gfv%sqi(n+1)
                tmp1 = n+n+l-1-x
                tmp2 = gfv%sq(n-1) * gfv%sq(n-1+l)
                tmp3 = n+n+l-x
                tmp4 = 2 * gfv%sq(n) * gfv%sq(n+l)
                WF%HO%ql(n,l,il) = (tmp1*WF%HO%ql(n-1,l,il) - tmp2*WF%HO%ql(n-2,l,il))*tmp0
                WF%HO%ql1(n,l,il) =  tmp3*WF%HO%ql(n,l,il)  -tmp4*WF%HO%ql(n-1,l,il)
            enddo
        enddo
        !------------basis for the bosons----------------
        xb = gauss%xl(il) * cb **2 ! bosons
        psiLB_00 = gfv%sq(2)*exp(-half*xb) / (sqrt(2*pi)*HO%b0*HO%bp) 
        WF%HO%qlb(0,il) = psiLB_00
        WF%HO%qlb(1,il) = psiLB_00 * (1-xb)
        do n = 2,nr_max_boson
            WF%HO%qlb(n,il) = ((2*n-1-xb) * WF%HO%qlb(n-1,il) - (n-1) * WF%HO%qlb(n-2,il)) / n
        enddo
    enddo

    if(ifPrintHO) call printHOWavefunction
    contains
    subroutine printHOWavefunction
        character(len=*),parameter :: format1 = "(a,i3,30f15.8)", &
                                      format2 = "(a,i4,i3,30f15.8)"
        integer :: n,l
        integer :: ix = 10
        open(u_config,file=outputfile%config,status='unknown',position='append')
        write(u_config,*) '*************************BEGIN calculate_HO_wavefuction ********************'
        do n = 0,min(HO%nzm,HO%nzbm)
            write(u_config,'(a,i2,a)') ' QH(nz=',n,',ih=1...)'
            write(u_config,format1) ' H     xh',n,(WF%HO%qh(n,ih),ih=1,ix)              
            write(u_config,format1) ' dH/dx xh',n,(WF%HO%qh1(n,ih),ih=1,ix)             
            write(u_config,format1) ' H  cb*xh',n,(WF%HO%qhb(n,ih)*sqrt(gauss%wh(ih)),ih=1,ix)
        enddo
        do l=0, HO%mlm
            write(u_config,*) '      nr ml    QL(nr,l,il=1,...)'
            do n=0,HO%nrm
                write(u_config,format2)'ql  ', n,l,(WF%HO%ql(n,l,il),il=1,ix) 
                write(u_config,format2)'ql1 ', n,l,(WF%HO%ql1(n,l,il),il=1,ix)
                if(l==0) then
                write(u_config,format2)'qlb ', n,l,(WF%HO%qlb(n,il),il=1,ix) 
                endif
            enddo
        enddo
        write(u_config,"(a,/)") '*************************END calculate_HO_wavefuction ********************'
        close(u_config)
    end subroutine printHOWavefunction
end subroutine calculate_HO_wavefuction

subroutine test_integration_HO_wavafuction
    real(r64):: s,sb,s2
    integer :: n1,n2,ih,l,il
    character(len=*),parameter :: format1 = "(' G.Hermit: n1 =',i3,' n2 =',i3, 3f15.8)",&
                                  format2 = "(' G.Lague.: l =' ,i2,' n1 =',i3,' n2 =',i3, 3f15.8)"
    ! test for hermite integration
    do n1 = 0,min(nz_max,nz_max_boson)                                        
         do n2 = 0,n1                                                                                 
               s  = 0.0d0                                               
               sb = 0.0d0                                                
               s2 = 0.0d0                                                 
               do ih = 1,gauss%nh                                            
                  s = s  + WF%HO%qh(n1,ih)*WF%HO%qh(n2,ih)                         
                  if (n1.le.nz_max_boson) sb = sb + WF%HO%qhb(n1,ih)*WF%HO%qhb(n2,ih)*gauss%wh(ih)             
               enddo                                                    
               write(*, format1) n1,n2,s,sb*HO%b0*HO%bz                                                                                        
         enddo                                                          
    enddo
    ! test for Laguerre integration
         do l = 0,ml_max                                                   
            do n1 = 0,nr_max                                               
                do n2 = 0,n1                                                
                s  = 0.0d0                                               
                sb = 0.0d0                                                
                s2 = 0.0d0                                                
                do il = 1,gauss%nl                                           
                    s = s + WF%HO%ql(n1,l,il)*WF%HO%ql(n2,l,il)                      
                    if (l.eq.0) sb = sb + WF%HO%qlb(n1,il)*WF%HO%qlb(n2,il)*gauss%wl(il)                
                enddo                                                    
                write(*,format2) l,n1,n2,s,sb*pi*(HO%b0*HO%bp)**2                                  
                enddo                                                       
            enddo                                                       
         enddo 
end subroutine test_integration_HO_wavafuction

END MODULE Wavefunctions