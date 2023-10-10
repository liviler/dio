!====================================================================================!
! MODULE Matrix                                                                      !
!                                                                                    !
! This module contains subroutines related to  
!                                                                                    !
!                                                                                    !
! List of subroutines and functions:                                                 !
! - subroutine                                                                       !
!====================================================================================!
MODULE Matrix
use Constants, only: r64,zero,one,half,hbc,u_config,nfx,ngx,ntb_max
use Globals, only: HO,gfv,gauss,WF,matrices,force,outputfile

implicit none

contains

subroutine calculate_sigma_nabla(ifPrint)
    !-----------------------------------------------------------------------------------------------------
    ! calculates single particle matrix elements:
    ! \int d\boldsymbol{r} \Psi_{\alpha} \boldsymbol{\sigma}\cdot\boldsymbol{\nabla} \Psi_{\tilde{\alpha}}
    ! for Fermions in the cylindrical oscillator basis
    !-----------------------------------------------------------------------------------------------------
    logical,intent(in),optional :: ifPrint
    integer :: ib,nf,ng,i0f,i0g,i2,i1, nz2,nr2,ml2,ms2,nz1,nr1,ml1,ms1, m,il
    real(r64) :: fz,fp,s
    fz = hbc/(HO%b0*HO%bz)
    fp = hbc/(HO%b0*HO%bp)
    do ib =1,HO%nb
        nf = HO%id(ib,1) ! dimension of large components in block ib
        ng = HO%id(ib,2) ! dimension of small components in block ib
        i0f = HO%ia(ib,1) ! begin of large components in block ib
        i0g = HO%ia(ib,2) ! begin of small components in block ib
        do i2 = 1,nf
            nz2 = HO%nz(i0f+i2)
            nr2 = HO%nr(i0f+i2)
            ml2 = HO%ml(i0f+i2)
            ms2 = HO%ms(i0f+i2)
            do i1 = 1,ng
                nz1 = HO%nz(i0g+i1)
                nr1 = HO%nr(i0g+i1)
                ml1 = HO%ml(i0g+i1)
                ms1 = HO%ms(i0g+i1)
                
                s = zero
                if(ms1.eq.ms2) then ! also ml1 equal ml2 with the same k
                  if (nr1.eq.nr2) then                                     
                    if (nz1.eq.nz2+1) s = -ms1*fz*gfv%sq(nz1)*gfv%sqi(2)          
                    if (nz1.eq.nz2-1) s = +ms1*fz*gfv%sq(nz2)*gfv%sqi(2)
                  endif
                else
                    if (nz1.eq.nz2) then
                        if (ml1.eq.ml2+1) then                                
                            m = -ml2                                           
                        else                                                  
                            m = +ml2                                           
                        endif                                                 
                        do il = 1,gauss%nl                                         
                            s = s + WF%HO%ql(nr1,ml1,il)*(WF%HO%ql1(nr2,ml2,il)+m*WF%HO%ql(nr2,ml2,il))/sqrt(gauss%xl(il))   
                        enddo                                                
                        s = s*fp
                    endif
                endif
                matrices%sp(i1+(i2-1)*ng,ib) = -s 
            enddo
        enddo
    enddo

    if(ifPrint) call printSIgmaNablaMatrices
    contains

    subroutine printSIgmaNablaMatrices
        integer :: ib,i,nf,ng,i0f,i0g
        character(len=8) :: tp(nfx),tm(ngx)
        open(u_config,file=outputfile%config,status='unknown',position='append')
        write(u_config,*) '*************************BEGIN calculate_sigma_nabla ********************'
        do ib = 1,HO%nb
            nf = HO%id(ib,1) ! dimension of large components in block ib
            ng = HO%id(ib,2) ! dimension of small components in block ib
            i0f = HO%ia(ib,1) ! begin of large components in block ib
            i0g = HO%ia(ib,2) ! begin of small components in block ib
            write(u_config,"(/,a)") HO%txb(ib)
            do i= 1,nf
                tp(i) = HO%tb(i+i0f)
            enddo
            do i = 1,ng
                tm(i) = HO%tb(i+i0g)
            enddo
            call aprint(u_config,1,3,6,ng,ng,nf,matrices%sp(1,ib),tm,tp,'Sigma * P')
        enddo
        write(u_config,"(a,/)") '*************************BEGIN calculate_sigma_nabla ********************'
        close(u_config)
    end subroutine printSIgmaNablaMatrices
end subroutine calculate_sigma_nabla

subroutine aprint(u_write,is,it,ns,ma,n1,n2,a,t1,t2,text)
    !----------------------------------------------------------------!
    !                                                                !
    !     IS = 1    Full matrix                                      !
    !          2    Lower diagonal matrix                            !
    !          3    specially stored symmetric matrix                !
    !                                                                !
    !     IT = 1    numbers for rows and columns                     !
    !          2    text for rows and numbers for columns            !
    !          3    text for rows and columns                        !
    !                                                                !
    !     NS = 1     FORMAT   8F8.4     80 Coulums                   !
    !     NS = 2     FORMAT   8f8.2     80 Coulums                   !
    !     NS = 3     FORMAT  17F4.1     80 Coulums                   !
    !     NS = 4     FORMAT  30F4..1    120 Coulums                  !
    !     NS = 5     FORMAT  5F12.8     80 Coulums                   !
    !     NS = 6     FORMAT  5F12.4     80 Coulums                   !
    !     NS = 7     FORMAT  4E13.6     80 Coulums                   !
    !     NS = 8     FORMAT  8E15.8    130 Coulums                   !
    !----------------------------------------------------------------!
    implicit double precision (a-h,o-z)
    implicit integer (i-n)
    integer :: u_write,is,it,ns,ma,n1,n2
    real(r64),dimension(ma*n2),intent(in) :: a
    character(len=8) :: t1(n1), t2(n2)
    character(len=*) :: text

    character(len=30) :: fmt1, fmt2
    character*20 :: fti,ftt,fmt(8),fmti(8),fmtt(8)
    integer:: nsp(8),nspalt

    data nsp/8,8,17,30,5,5,4,8/
    data fmt /'8f8.4)',            '8F8.2)',   &         
              '17f4.1)',           '30f4.1)',  &         
              '5f12.8)',           '5f12.4)',  &        
              '4e13.6)',           '8e15.8)'/           
    data fmti/'(11x,8(i4,4x))',    '(11x,8(i4,4x))',    &
              '(11x,17(1x,i2,1x))','(11x,30(1x,i2,1x))',&
              '(11x,6(i4,8x))',    '(11x,10(i4,8x))',   &
              '(11x,5(i4,9x))',    '(11x,8(i4,11x))'/   
    data fmtt/'(11x,8a8)',         '(11x,8a8)',       &  
              '(11x,17a4)',        '(11x,30a4)',      &  
              '(11x,6(a8,2x))',    '(11x,10(4x,a8))', &  
              '(11x,5(a8,5x))',    '(11x,8(a8,7x))'/

    fmt1   = '(4x,i3,4x,' // fmt(ns)
    fmt2   = '(1x,a8,2x' // fmt(ns) 
    fti    = fmti(ns)               
    ftt    = fmtt(ns)               
    nspalt = nsp(ns)
    write(u_write,'(//,3x,a)') text
    ka = 1                           
    ke = nspalt  ! end coulum number to display in first Part                      
    nteil = n2/nspalt 
    if (nteil*nspalt.ne.n2) nteil = nteil + 1
    do 10  nt = 1,nteil                                            
        if (n2.gt.nspalt)  write(u_write,100)  nt                       
 100    format(//, 10x,'Part',i5,' of the Matrix',/)                
        if (nt.eq.nteil) ke = n2    ! in last Part, coulum end to n2                                  
        if (it.lt.3) then                                           
            write(u_write,fti) (k,k=ka,ke) ! number for coulums                               
        else                                                        
            write(u_write,ftt) (t2(k),k=ka,ke)  ! text for coulums                          
        endif  
                                                    
        do 20  i=1,n1                                                  
            kee=ke                                                      
            if (is.eq.2.and.ke.gt.i) kee=i                              
            if (ka.gt.kee) goto 20                                      
            if (is.eq.3) then                                           
                if (it.eq.1) then                                        
                    write(u_write,fmt1) i,(a(i+(k-1)*(n1+n1-k)/2),k=ka,kee)
                else                                                     
                    write(u_write,fmt2) t1(i),(a(i+(k-1)*(n1+n1-k)/2),k=ka,kee)
                endif                                                    
            else                                                        
                if (it.eq.1) then                                        
                    write(u_write,fmt1) i,(a(i+(k-1)*ma),k=ka,kee)             
                else                                                     
                    write(u_write,fmt2) t1(i),(a(i+(k-1)*ma),k=ka,kee)         
                endif                                                    
            endif                                                       
    20  continue                                                       

    ka=ka+nspalt                                                   
    ke=ke+nspalt                                                   
 10 continue
end subroutine aprint

subroutine calculate_meson_propagators(ifPrint) 
    !-------------------------------------------------------------------
    ! Calculates the meson-propagators GG.                              
    ! DD =/nabla^2 is the Laplace operator in oscillator space                   
    ! GG = m**2/(-DD+m**2)
    ! Actually, what is calculated is (-DD+m**2)/m**2 in Oscillator Basis
    !--------------------------------------------------------------------
    logical,intent(in),optional :: ifPrint
    integer :: i2,nz2,nr2,i1,nz1,nr1,no,imes,k,i,ifl
    real(r64),dimension(4) :: mass
    real(r64),dimension(ntb_max*ntb_max) ::dd,aa
    real(r64),dimension(ntb_max*ntb_max,4) :: gg
    real(r64) :: az2,ap2,t,m2i,d

    mass(1) = force%masses%amsig
    mass(2) = force%masses%amome
    mass(3) = force%masses%amdel
    mass(4) = force%masses%amrho

    az2 = one/(HO%b0*HO%bz)**2
    ap2 = one/(HO%b0*HO%bp)**2

    no = HO%ntb
    do i2 = 1,HO%ntb
        nz2 = HO%nzb(i2)
        nr2 = HO%nrb(i2)
        dd(i2+(i2-1)*no) = -(az2*(nz2+half)+ap2*(nr2+nr2+1))
        do i1 = 1, i2-1
            nz1 = HO%nzb(i1)
            nr1 = HO%nrb(i1)
            t   = zero
            if (nr1.eq.nr2) then
               if (nz2.eq.nz1+2) t = az2*half*gfv%sq(nz1+1)*gfv%sq(nz2)
               if (nz2.eq.nz1-2) t = az2*half*gfv%sq(nz1)*gfv%sq(nz2+1)
            endif
            if (nz1.eq.nz2) then
               if (nr2.eq.nr1+1) t = -ap2*nr2
               if (nr2.eq.nr1-1) t = -ap2*nr1
            endif
            dd(i1+(i2-1)*no) = t
            dd(i2+(i1-1)*no) = t
        enddo
    enddo
                                                    
    do imes = 1,4                                                     
        m2i = (one/(mass(imes)+1.d-10))**2                                     
        do k = 1,no                                                    
            do i = 1,no                                                 
                aa(i+(k-1)*no) = -dd(i+(k-1)*no)*m2i              
                gg(i+(k-1)*no,imes) = zero                               
            enddo                                                       
            aa(k+(k-1)*no) = aa(k+(k-1)*no) + one                  
            gg(k+(k-1)*no,imes) = one                                
        enddo                                                          
        call lingd(no,no,no,no,aa,gg(1,imes),d,ifl)
    enddo
    if(ifl == -1) stop '[Field]: error in subroutine lingd!'
    matrices%meson_propagators(:,:) = gg(:,:)

    if(ifPrint) call printMesonPropagator
    contains

    subroutine printMesonPropagator
        character(len=8) :: t(no)
        t(1:no) = ''
        open(u_config,file=outputfile%config,status='unknown',position='append')
        write(u_config,*) '*************************BEGIN calculate_meson_propagators ********************'
        do imes =1,4
            call aprint(u_config,1,1,1,no,no,no,gg(1,imes),t,t,'Meson-Propagator')
        enddo
        write(u_config,"(a,/)") '*************************BEGIN calculate_meson_propagators ********************'
        close(u_config)
    end subroutine printMesonPropagator

end subroutine calculate_meson_propagators

subroutine lingd(ma,mx,n,m,a,x,d,ifl)
    !----------------------------------------------------------------------
    ! Solves the system of linear equations A*X = B .
    ! At the beginning, the matrix B is stored in X. During the calculation 
    ! it will be overwritten. D is the determinant of A .
    !----------------------------------------------------------------------
    integer :: ma, mx, n, m, ifl
    real(r64),dimension(ma,m) :: a, x
    real(r64) :: d

    integer :: i,j,k,l,k1,n1
    real(r64) :: p,q,tol,cp,cq

    ! constant
    real(r64) :: tollim = 1.d-10, one = 1.d0, zero = 0.d0
    ifl = 1
    p = zero
    do i = 1,n
        q = zero
        do j = 1,n
        q = q + abs(a(i,j))
        enddo
        if (q.gt.p)   p=q
    enddo
    tol = tollim*p
    d   = one
    do k = 1,n
        p = zero
        do j = k,n
            q = abs(a(j,k))
            if (q.ge.p) then
                p = q
                i = j
            endif
        enddo
        if (p.le.tol) then
!             write (*,100) ('-',j=1,80),tol,i,k,a(i,k),('-',j=1,80)
! 100         format(1x,80a1,' *****  ERROR IN LINGD , TOLERANZ =', e10.4,' VALUE OF A(',i3,',',i3,') IS ',e10.4/1x,80a1) ! check it !
            ifl = -1
            return
        endif
        cp = one/a(i,k)
        if (i.ne.k) then
            d = -d
            do l = 1,m
                cq     = x(i,l)
                x(i,l) = x(k,l)
                x(k,l) = cq
            enddo
            do l = k,n
                cq     = a(i,l)
                a(i,l) = a(k,l)
                a(k,l) = cq
            enddo
        endif
        d = d*a(k,k)
        if (k.eq.n) goto 1
        k1 = k+1
        do i = k1,n
            cq = a(i,k)*cp
            do l = 1,m
                x(i,l) = x(i,l) - cq*x(k,l)
            enddo
            do l = k1,n
                a(i,l) = a(i,l) - cq*a(k,l)
            enddo
        enddo
    enddo
  1 do l = 1,m
        x(n,l) = x(n,l)*cp
    enddo
    if (n.gt.1) then
        n1 = n-1
        do k=1,n1
            cp = one/a(n-k,n-k)
            do l=1,m
                cq = x(n-k,l)
                do i = 1,k
                    cq = cq - a(n-k,n+1-i)*x(n+1-i,l)
                enddo
                x(n-k,l) = cq*cp
            enddo
        enddo
    endif
    return 
end subroutine lingd

END MODULE Matrix