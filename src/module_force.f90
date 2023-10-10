!==============================================================================!
! MODULE Force                                                                 !
!                                                                              !
! This module sets the selected force parameters                               !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_force_parameter                                             !
!==============================================================================!
MODULE Forces

use Constants, only: r64, zero, one,two, hbc, u_config
use Globals, only: force, outputfile, pairing
implicit none

contains
subroutine set_force_parameters(ifPrint)
    logical,intent(in),optional :: ifPrint
    real(r64), dimension(5) :: temp_array
    real(r64) :: alfs, bets, gams, dels, alfv, gamv,delv, alftv, deltv, alfts, delts
    !Statement Function
    ! set_a_s(c_s,d_s,b_s,d_s)  = (one + c_s*(one+d_s)**2)/(one + b_s*(one+d_s)**2)

    !===============================================================
    !     HS:  Horowitz and Serot:
    !---------------------------------------------------------------
    if (force%parname.eq.'HS') then
        force%masses%amu    = 939.0d0
        force%masses%amsig  = 520.0d0
        force%masses%amome  = 783.0d0
        force%masses%amrho  = 770.0d0
        force%couplm%gsig   = 10.47d0
        force%couplm%gome   = 13.80d0
        force%couplm%grho   = 8.07d0
        force%nonlin%g2     = zero
        force%nonlin%g3     = zero
        force%option%inl    = 0
        force%option%ipc    = 0
        force%option%idd    = 0
    endif   ! HS
    !=============================================================== 
    !     NL3:
    !---------------------------------------------------------------
    if (force%parname.eq.'NL3') then
        force%masses%amu    = 939.d0
        force%masses%amsig  = 508.194d0
        force%masses%amome  = 782.501d0
        force%masses%amrho  = 763.0d0 
        force%couplm%gsig   = 10.217d0
        force%couplm%gome   = 12.868d0
        force%couplm%grho   = 4.474d0 
        force%nonlin%g2     = -10.431d0
        force%nonlin%g3     = -28.885d0
        force%option%inl    = 1
        force%option%ipc    = 0
        force%option%idd    = 0
    endif   ! NL3
    !===============================================================
    if (force%parname.eq.'NL1') then
        force%masses%amu    = 938.d0
        force%masses%amsig  = 492.25d0
        force%masses%amome  = 795.359d0
        force%masses%amrho  = 763.0d0
        force%couplm%gsig   = 10.138d0
        force%couplm%gome   = 13.285d0
        force%couplm%grho   = 4.975d0
        force%nonlin%g2     = -12.172d0  
        force%nonlin%g3     = -36.265d0
        force%option%inl    = 1
        force%option%ipc    = 0
        force%option%idd    = 0
    endif   ! NL1
    !=============================================================== 
    if (force%parname.eq.'DD-ME1') then
    !---------------------------------------------------------------
        force%masses%amu    =  939.d0                   ! MeV
        force%masses%amsig  = 549.525547d0              ! MeV
        force%masses%amome  = 783.d0                    ! MeV
        force%masses%amrho  = 763.d0                    ! MeV
        force%couplm%gsig   = 10.443356d0
        force%couplm%gome   = 12.8939219d0
        force%couplm%grho   =  3.8053377d0
        force%dtypel%b_s    =  0.9780577d0
        force%dtypel%c_s    =  1.5342d0
        force%dtypel%c_v    =  1.3566d0
        force%dtypel%a_tv   =  0.500768d0
        force%dtypel%dsat   =  0.152d0
        temp_array= typel(force%dtypel%b_s,force%dtypel%c_s, force%dtypel%c_v, force%dtypel%a_tv)
        force%dtypel%a_s    = temp_array(1)
        force%dtypel%d_s    = temp_array(2)
        force%dtypel%a_v    = temp_array(3)
        force%dtypel%b_v    = temp_array(4)
        force%dtypel%d_v    = temp_array(5)
        force%rosat         = force%dtypel%dsat 
        force%option%idd    =  1
        force%option%ipc    =  0
        force%option%inl    =  0
        
    endif   ! DD-ME1
    !=============================================================== 
    if (force%parname.eq.'DD-ME2') then
    !---------------------------------------------------------------
        force%masses%amu    = 939.d0                   ! MeV
        force%masses%amsig  = 550.123780d0             ! MeV
        force%masses%amome  = 783.d0                   ! MeV
        force%masses%amrho  = 763.d0                   ! MeV
        force%masses%amdel  = 1000.d0
        force%couplm%gsig   = 10.53958d0
        force%couplm%gome   = 13.018913d0
        force%couplm%grho   = 3.683618d0
        force%couplm%gdel   = 0.d0
        force%dtypel%b_s    = 1.094325d0
        force%dtypel%c_s    = 1.705658d0
        force%dtypel%c_v    = 1.461999d0
        force%dtypel%a_tv   = 0.564729d0
        force%dtypel%dsat   = 0.152d0
        temp_array= typel(force%dtypel%b_s,force%dtypel%c_s, force%dtypel%c_v, force%dtypel%a_tv)
        force%dtypel%a_s    = temp_array(1)
        force%dtypel%d_s    = temp_array(2)
        force%dtypel%a_v    = temp_array(3)
        force%dtypel%b_v    = temp_array(4)
        force%dtypel%d_v    = temp_array(5)
        force%rosat         = force%dtypel%dsat
        force%option%idd    =  1    
        force%option%ipc    =  0
        force%option%inl    =  0     
    endif   ! DD-ME2
    !=============================================================== 
    if (force%parname.eq.'PC-F1') then
    !---------------------------------------------------------------
        force%masses%amu   =  939.d0                   ! MeV
        alfs   =  -3.83577d-4
        bets   =  +7.68567d-11
        gams   =  -2.90443d-17
        dels   =  -4.18530d-10
        alfv   =  +2.59333d-4
        gamv   =  -3.87900d-18
        delv   =  -1.19210d-10
        alftv  =  +3.46770d-5
        deltv  =  -4.20000d-11
        alfts  =   0.d0
        delts  =   0.d0

        !--------in unit of fm^2
        !------- quadratic terms
        force%couplg%ggsig  = alfs*hbc**2
        force%couplg%ggome  = alfv*hbc**2
        force%couplg%ggdel  = alfts*hbc**2 
        force%couplg%ggrho  = alftv*hbc**2
        !------- derivative terms         
        force%coupld%ddsig  = dels*hbc**4
        force%coupld%ddome  = delv*hbc**4
        force%coupld%ddrho  = deltv*hbc**4
        force%coupld%dddel  = delts*hbc**4
        !------- non-linear terms
        force%coupnl%ggbet  = bets*hbc**5 
        force%coupnl%gggams = gams*hbc**8
        force%coupnl%gggamv = gamv*hbc**8         
        force%dtypel%dsat   = 0.152d0
        force%rosat         = force%dtypel%dsat
        force%option%idd    =  0
        force%option%ipc    =  1
        force%option%inl    =  1
        !------- new combinations of parameter
        force%dpolyn%alf(1) =  2 * force%coupnl%ggbet / (3 * force%couplg%ggsig)
        force%dpolyn%alf(2) =  zero
        force%dpolyn%alf(3) =  zero
        force%dpolyn%alf(4) =  zero

        force%dpolyn%bet(1) =  force%coupnl%gggams / (2 * force%couplg%ggsig)
        force%dpolyn%bet(2) =  force%coupnl%gggamv/(2 * force%couplg%ggome)
        force%dpolyn%bet(3) =  zero
        force%dpolyn%bet(4) =  zero         
      endif   !  PC-F1
    !===============================================================
    if (force%parname.eq.'PC-PK1') then
    !---------------------------------------------------------------
        force%masses%amu    =  939.d0                   ! MeV
        alfs   =  -3.96291d-4
        bets   =  +8.6653d-11
        gams   =  -3.80724d-17
        dels   =  -1.09108d-10
        alfv   =  +2.69040d-4
        gamv   =  -3.64219d-18
        delv   =  -4.32619d-10
        alfts  =   0.d0
        delts  =   0.d0
        alftv  =  +2.95018d-5
        deltv  =  -4.11112d-10
    
        !--------in unit of fm^2
        !------- quadratic terms
        force%couplg%ggsig  = alfs*hbc**2
        force%couplg%ggome  = alfv*hbc**2
        force%couplg%ggdel  = alfts*hbc**2
        force%couplg%ggrho  = alftv*hbc**2
        !------- derivative terms
        force%coupld%ddsig  = dels*hbc**4
        force%coupld%ddome  = delv*hbc**4
        force%coupld%ddrho  = deltv*hbc**4
        force%coupld%dddel  = delts*hbc**4
        !------- non-linear terms
        force%coupnl%ggbet  = bets*hbc**5
        force%coupnl%gggams = gams*hbc**8
        force%coupnl%gggamv = gamv*hbc**8
        force%dtypel%dsat   = 0.152d0
        force%rosat         = force%dtypel%dsat
        force%option%idd    =  0
        force%option%ipc    =  1
        force%option%inl    =  1
        !------- new combinations of parameter
        force%dpolyn%alf(1) = 2 * force%coupnl%ggbet / (3 * force%couplg%ggsig)
        force%dpolyn%alf(2) = zero
        force%dpolyn%alf(3) = zero
        force%dpolyn%alf(4) = zero

        force%dpolyn%bet(1) = force%coupnl%gggams / (2 * force%couplg%ggsig)
        force%dpolyn%bet(2) = force%coupnl%gggamv/(2 * force%couplg%ggome)
        force%dpolyn%bet(3) = zero
        force%dpolyn%bet(4) = zero
    endif   !  PC-PK1
    !===============================================================       
    if(force%parname .eq. 'DD-PC1') then
    !---------------------------------------------------------------      
        force%masses%amu    = 939.d0
        force%couplm%gsig   = -10.04616d0
        force%couplm%gome   = 5.91946d0
        force%couplm%grho   = 0.d0
        
        force%dtypel%a_s    = -9.15042d0
        force%dtypel%b_s    = -6.42729d0
        force%dtypel%c_s    = 0.d0
        force%dtypel%d_s    = 1.37235d0
            
        force%dtypel%a_v    = 8.86370d0
        force%dtypel%b_v    = 0.d0
        force%dtypel%c_v    = 0.d0
        force%dtypel%d_v    = 0.65835d0
            
        force%dtypel%a_tv   = 1.84210d0
        force%dtypel%b_tv   = 0.d0
        force%dtypel%c_tv   = 0.d0
        force%dtypel%d_tv   = 0.64352d0
            
        dels  =-0.8149d0
        delv  = 0.d0
        deltv = 0.d0
        delts = 0.d0
            
        force%couplg%ggsig  = 1.d0
        force%couplg%ggome  = 1.d0
        force%couplg%ggrho  = 1.d0
        force%coupld%ddsig  = dels
        force%coupld%ddome  = delv
        force%coupld%ddrho  = deltv
        force%coupld%dddel  = delts

        force%dtypel%dsat   =    0.152d0
        force%rosat         =    force%dtypel%dsat
            
        force%option%idd = 2
        force%option%ipc = 1
        force%option%inl = 0
    endif
    !----------------------------------------------------------------------
    !------Delta force  
    !       vpair(1) = -326.d0
    !       vpair(2) = -326.d0
    !       ecut     = 18.d0
    !       dcut     = 0.5d0
    !=============================================================== 
    
    !---- masses in units of fm**(-1)
    force%masses%amu      = force%masses%amu/hbc
    force%masses%amsig    = force%masses%amsig/hbc
    force%masses%amome    = force%masses%amome/hbc
    force%masses%amdel    = force%masses%amdel/hbc
    force%masses%amrho    = force%masses%amrho/hbc
    force%masses%ampi     = force%masses%ampi/hbc
    
    if (force%option%ipc.eq.0) then
    force%couplg%ggsig  = -(force%couplm%gsig/(force%masses%amsig+1.d-10))**2
    force%couplg%ggome  = +(force%couplm%gome/(force%masses%amome+1.d-10))**2
    force%couplg%ggdel  = -(force%couplm%gdel/(force%masses%amdel+1.d-10))**2
    force%couplg%ggrho  = +(force%couplm%grho/(force%masses%amrho+1.d-10))**2
    if (abs(force%couplm%gsig).gt.1.d-5) force%couplg%lmes(1) = 1
    if (abs(force%couplm%gome).gt.1.d-5) force%couplg%lmes(2) = 1
    if (abs(force%couplm%gdel).gt.1.d-5) force%couplg%lmes(3) = 1
    if (abs(force%couplm%grho).gt.1.d-5) force%couplg%lmes(4) = 1
    endif

    if(ifPrint) call printForces
    contains
    function typel(b_s, c_s, c_v, a_tv)
        !prepares parameters of the density dependence
        real(r64), intent(in) :: b_s, c_s, c_v, a_tv
        real(r64) :: a_s, d_s, a_v, b_v, d_v, x, facs, faco, fac1, fac2
        real(r64), dimension(5) :: typel
        a_s  = (one + c_s*(one+d_s)**2)/(one + b_s*(one+d_s)**2)
        d_s  = one/sqrt(3.d0*c_s)
        d_v  = one/sqrt(3.d0*c_v)
        facs = two*a_s*(b_s - c_s)*(one - 3.d0*c_s*(one+d_s)**2)/(one + c_s*(1+d_s)**2)**3
        faco = (one - 3.d0*c_v*(one + d_v)**2)/(one + c_v*(1 + d_v)**2)**3
        x    = facs/(two*faco)
        fac1 = x + c_v*(one + c_v*(one + d_v)**2)
        fac2 = one + c_v*(one+d_v)**2-x*(one+d_v)**2
        b_v  = fac1/fac2
        a_v  =(one + c_v*(one + d_v)**2)/(one + b_v*(one + d_v)**2)
        typel = [a_s, d_s, a_v, b_v, d_v]
    end function typel

    subroutine printForces
        character(len=*), parameter :: format1 = "(a10,f10.6)", &
                                       format2 = "(3(a10,f10.6,5x))", &
                                       format3 = "(4(a10,f10.6,5x))"

        open(u_config, file=outputfile%config, status='unknown', position='append')
        write(u_config,*) '*************************BEGIN set_force_parameter********************'
        write(u_config,*) 'Force Name: ',force%parname
        ! non-linear meson exchange
        if(force%option%ipc.eq.0 .and. force%option%idd.eq.0 .and. force%option%inl.eq.1) then       
            write(u_config,format1) 'mu     = ',force%masses%amu*hbc
            write(u_config,format2) 'msig   = ',force%masses%amsig*hbc, &
                                    'mome   = ',force%masses%amome*hbc, &
                                    'mrho   = ',force%masses%amrho*hbc
            write(u_config,format2) 'gsig   = ',force%couplm%gsig, &
                                    'gome   = ',force%couplm%gome, &
                                    'grho   = ',force%couplm%grho
            write(u_config,format2) 'g2     = ',force%nonlin%g2, &
                                    'g3     = ',force%nonlin%g3, &
                                    'c3     = ',force%nonlin%c3
        ! typel-wolter, meson exchange
        else if(force%option%ipc.eq.0 .and. force%option%idd.eq.1 .and. force%option%inl.eq.1) then
            write(u_config,format1) 'mu     = ',force%masses%amu*hbc
            write(u_config,format2) 'msig   = ',force%masses%amsig*hbc, &
                                    'mome   = ',force%masses%amome*hbc, &
                                    'mrho   = ',force%masses%amrho*hbc
            write(u_config,format2) 'gsig   = ',sqrt(-force%couplg%ggsig)*force%masses%amsig, &
                                    'gome   = ',sqrt(-force%couplg%ggome)*force%masses%amome, &
                                    'grho   = ',sqrt(-force%couplg%ggrho)*force%masses%amrho
            write(u_config,format3) 'a_s    = ',force%dtypel%a_s, &
                                    'b_s    = ',force%dtypel%b_s, &
                                    'c_s    = ',force%dtypel%c_s, &
                                    'd_s    = ',force%dtypel%d_s
            write(u_config,format3) 'a_v    = ',force%dtypel%a_v, &
                                    'b_v    = ',force%dtypel%b_v, &
                                    'c_v    = ',force%dtypel%c_v, &
                                    'd_v    = ',force%dtypel%d_v
            write(u_config,format3) 'a_tv   = ',force%dtypel%a_tv, &
                                    'b_tv   = ',force%dtypel%b_tv, &
                                    'c_tv   = ',force%dtypel%c_tv, &
                                    'd_tv   = ',force%dtypel%d_tv
            write(u_config,format1) 'rosat  = ',force%rosat
        ! Buervenich
        else if(force%option%ipc.eq.1 .and. force%option%idd.eq.0 .and. force%option%inl.eq.1) then
            write(u_config,format1) 'mu     = ',force%masses%amu*hbc
            write(u_config,format3) 'alps   = ',force%couplg%ggsig, &
                                    'alpv   = ',force%couplg%ggome, &
                                    'alptv  = ',force%couplg%ggdel, &
                                    'alpts  = ',force%couplg%ggrho
            write(u_config,format3) 'dels   = ',force%coupld%ddsig, &
                                    'delv   = ',force%coupld%ddome, &
                                    'delts  = ',force%coupld%dddel, &
                                    'deltv  = ',force%coupld%ddrho
            write(u_config,format2) 'bets   = ',force%coupnl%ggbet, &
                                    'gams   = ',force%coupnl%gggams,&
                                    'gamv   = ',force%coupnl%gggamv
        ! point-coupling, density-dependent
        else if(force%option%ipc.eq.1 .and. force%option%idd.eq.2 .and. force%option%inl.eq.0) then
            write(u_config,format1) 'mu     = ',force%masses%amu*hbc
            write(u_config,format2) 'gsig   = ',force%couplm%gsig, &
                                    'gome   = ',force%couplm%gome, &
                                    'grho   = ',force%couplm%grho
            write(u_config,format3) 'a_s    = ',force%dtypel%a_s, &
                                    'b_s    = ',force%dtypel%b_s, &
                                    'c_s    = ',force%dtypel%c_s, &
                                    'd_s    = ',force%dtypel%d_s
            write(u_config,format3) 'a_v    = ',force%dtypel%a_v, &
                                    'b_v    = ',force%dtypel%b_v, &
                                    'c_v    = ',force%dtypel%c_v, &
                                    'd_v    = ',force%dtypel%d_v
            write(u_config,format3) 'a_tv   = ',force%dtypel%a_tv, &
                                    'b_tv   = ',force%dtypel%b_tv, &
                                    'c_tv   = ',force%dtypel%c_tv, &
                                    'd_tv   = ',force%dtypel%d_tv
            write(u_config,format1) 'rosat  = ',force%rosat
        ! stop, wrong combination of parameters
        else 
            write(u_config,*)'Stop, wrong combination of parameters!'
             stop '[Force]set_force_parameter: ipc not properly defined'
        endif
        write(u_config,*) 'Pairing strength(Vpair): ',pairing%vpair(1),pairing%vpair(2)
        write(u_config,*) 'Cut-off parameters: '
        write(u_config,"(a,/)") '*************************END set_force_param**************************'
    end subroutine printForces
end subroutine set_force_parameters

END MODULE Forces