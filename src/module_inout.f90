!==============================================================================!
! MODULE Inout                                                                 !
!                                                                              !
! This module contains functions and routines for reading and writing files.   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE Inout
use Constants, only: i8,i16,i64,r64,u_start,u_config,pi,ngh,ngl,nghl,itx
use Globals,only : constraint,HO,woodssaxon,iteration,inin,&
                   nucleus_attributes,pairing,force,outputfile,erwar,fermi,fields
implicit none
integer, private :: u_b23 = u_start + 1  ! the unit of b23.dat
integer, private :: u_dio = u_start + 2  ! the unit of dio.dat
integer, private :: u_wfs = u_start + 3  ! the unit of  ???
contains

subroutine read_file_b23
    character(len=*), parameter :: file_path_b23 = './input/b23.dat' ! the path of b23.dat
    integer(i64) :: fileb23_lines_max = 100   ! maximum number of lines in file b23.dat

    integer :: length_count = 0 ! count the exact number of data rows in the 'b23.dat' file
    integer :: iostat ! store the io status of the 'read' function 
    integer(i64) :: index

    allocate(constraint%betac(fileb23_lines_max))
    allocate(constraint%bet3c(fileb23_lines_max))
    allocate(constraint%clam2(fileb23_lines_max))

    open(u_b23, file=file_path_b23, status='old')
    do index = 1, fileb23_lines_max+1
        if(length_count >= fileb23_lines_max) then
            write(*,*) "[Initialization]: 'fileb23_lines_max' too small"
            stop
        end if 
        read(u_b23,*,iostat=iostat) constraint%betac(index),constraint%bet3c(index),constraint%clam2(index)
        if (iostat /= 0) then
            if (iostat < 0) then
                ! read the end
                exit
            else
                ! error
                write(*, *) "Error reading from b23.dat!","iostat:",iostat
                stop
            end if
        end if
        length_count = length_count + 1
    end do
    constraint%length = length_count
    close(u_b23)
end subroutine read_file_b23

subroutine read_file_dio(ifPrint)
    logical,intent(in),optional :: ifPrint
    character(len=*), parameter :: file_path_dio = './input/dio.dat'

    character :: first_character
    real(r64) :: tmp
    character(len=*), parameter ::  format1= "(10x, i5)", &
                                    format2= "(10x, 2i3)", &
                                    format3= "(10x, f10.6)", &
                                    format4= "(10x, 2f10.6)", &
                                    format5= "(12x, a10)", &
                                    format6= "(a1, 9x, 2f10.3)", &
                                    format7= "(a2, i4)"

    open(u_dio, file=file_path_dio, status='old')
    read(u_dio, format2) HO%n0f, HO%n0b
    read(u_dio, format3) HO%b0
    read(u_dio, format6) first_character, tmp
    if(first_character == 'q') then
        HO%q = tmp
        HO%beta0 = dlog(HO%q)/(3*sqrt(5/(16*pi)))
    else
        HO%beta0 = tmp
        HO%q = exp(HO%beta0 * 3*sqrt(5/(16*pi)))
    endif
    read(u_dio, format6) first_character, tmp
    if(first_character == 'q') then
        woodssaxon%qs = tmp
        woodssaxon%beta2 = dlog(woodssaxon%qs)/(3*sqrt(5/(16*pi)))
    else
        woodssaxon%beta2 = tmp
        woodssaxon%qs = exp(woodssaxon%beta2 * 3*sqrt(5/(16*pi)))
    endif
    read(u_dio, format3) woodssaxon%beta3
    read(u_dio, format1) iteration%iteration_max
    read(u_dio, format3) iteration%xmix
    read(u_dio, format1) inin
    read(u_dio, format7) nucleus_attributes%name, nucleus_attributes%mass_number_int
    read(u_dio, format2) pairing%ide
    read(u_dio, format4) pairing%dec
    read(u_dio, format4) pairing%ga
    read(u_dio, format4) pairing%del
    read(u_dio, format5) force%parname
    read(u_dio, format4) pairing%vpair
    read(u_dio, format1) constraint%icstr
    read(u_dio, format3) constraint%cspr
    read(u_dio, format3) constraint%cmax
    close(u_dio)
    if(ifPrint) call printInputConfig
    contains
    subroutine printInputConfig
        character(len=*), parameter :: format11 = "(a,2i5)", &
                                       format12 = "(a,2f12.6)", &
                                       format13 = "(a,f8.5,5x,a,f8.5)", &
                                       format14 = "(a,i6)", &
                                       format15 = "(a,20x,a,a2,i4)"

        open(u_config,file=outputfile%config,status='unknown')
        write(u_config,*)'*************************BEGIN read_file_dio*************************'
        write(u_config,format11) ' Number of oscillator shells : ',HO%n0f,HO%n0b
        write(u_config,format12) ' Oscillator length b0 (fm)   : ',HO%b0  
        write(u_config,format13) ' Basis deformation           :  beta0 =',HO%beta0,'q =',HO%q  
        write(u_config,format13) ' Initial deformation         :  betas =',woodssaxon%beta2,'qs =',woodssaxon%qs       
        ! write(u_config,format12) ' Q1/Q2=                      : ',delta3
        write(u_config,format14) ' Maximal number of iterations: ',iteration%iteration_max         
        write(u_config,format12) ' Mixing parameter            : ',iteration%xmix             
        if (inin.eq.0) then                                                    
            write(u_config,*) 'Initial wavefunctions       :  will read from tape'  
        else if (inin.eq.1) then                                                    
            write(u_config,*) 'Initial wavefunctions       :  Saxon-Woods'  
        endif        
        write(u_config,format15) ' Nucleus ',':     ',nucleus_attributes%name,nucleus_attributes%mass_number_int      
        write(u_config,format11) ' Pairing control number      : ',pairing%ide
        write(u_config,format12) ' Frozen Gap Parameters       : ',pairing%dec               
        write(u_config,format12) ' Pairing-Constants           : ',pairing%ga                
        write(u_config,format12) ' Initial values for Gap      : ',pairing%del               
        write(u_config,format12) ' Vpair for delta force       : ',pairing%vpair
        write(u_config,format14) ' Calculation with constraint : ',constraint%icstr
        ! write(u_config,format12) ' Constraint beta2-values     : ',constraint%betac
        ! write(u_config,format12) ' Constraint beta3-values     : ',constraint%bet3c
        write(u_config,format12) ' Spring constant cspr        : ',constraint%cspr
        write(u_config,format12) ' cutoff for dE/db            : ',constraint%cmax
        write(u_config,"(a,/)")'*************************END read_file_dio****************************'
        close(u_config)
    end subroutine printInputConfig
end subroutine read_file_dio

subroutine set_output_filename(constraint_beta2,constraint_beta3)
    real(r64) :: constraint_beta2,constraint_beta3
    character :: sign_beta2, sign_beta3
    real(r64) :: abs2c, abs3c
    integer(i16), dimension(6) :: name
    integer(i16) :: name_nf1,name_nf2
    if (constraint_beta2 >= 0.d0) then
        sign_beta2 = '+'
    else
        sign_beta2 = '-'
    end if
    if (constraint_beta3 >= 0.d0) then
        sign_beta3 = '+'
    else
        sign_beta3 = '-'
    end if
    abs2c = abs(constraint_beta2)
    abs3c = abs(constraint_beta3)
    name(1) = abs2c + 48 !In ASCII, character '0' start from 48. 
    name(2) = mod(abs2c*10,10.d0)+48
    name(3) = mod(abs2c*100,10.d0)+48
    name(4) = abs3c+48
    name(5) = mod(abs3c*10,10.d0)+48
    name(6) = mod(abs3c*100,10.d0)+48
    name_nf1 = mod(HO%n0f/10,10) + 48
    name_nf2 = mod(HO%n0f,10) + 48
    ! the structure of ouput filenames are `dio_eMax`//HO%n0f//constraint_beta2*100//constraint_beta3*100//type
    ! like 'dio_eMax08+140+080.out' means HO%n0f=08, constraint_beta2= +1.40, constraint_beta3= +0.80, type is '.out'
    outputfile%outputf='dio'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //sign_beta2//char(name(1))//char(name(2))//char(name(3)) &
                        //sign_beta3//char(name(4))//char(name(5))//char(name(6))//'.out'
    outputfile%outputw='dio'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //sign_beta2//char(name(1))//char(name(2))//char(name(3)) &
                        //sign_beta3//char(name(4))//char(name(5))//char(name(6))//'.wel'
     outputfile%outdel='dio'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //sign_beta2//char(name(1))//char(name(2))//char(name(3)) &
                        //sign_beta3//char(name(4))//char(name(5))//char(name(6))//'.del'
     outputfile%outputwf='dio'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //sign_beta2//char(name(1))//char(name(2))//char(name(3)) &
                        //sign_beta3//char(name(4))//char(name(5))//char(name(6))//'.wf'
end subroutine set_output_filename

!--------------------------------------------------------------------------------------!
!subroutine read_fields:                                                               ! 
!reading of meson fields from the specified file                                       !
!--------------------------------------------------------------------------------------!
subroutine read_fields
    character(len=500) :: WFS_DIR
    character(len=*),parameter :: format1 = "(1x,a2,8i4)", &
                                  format2 = "(2(6x,f12.6),6x,f12.8)",&
                                  format3 = "(2(6x,f12.6),6x,f20.12)",&
                                  format4 = "(10x,5f10.4)", &
                                  format5 = "(10x,5f12.6)", &
                                  format6 = "(4e20.12)"
    character(len=2) :: nucleus_name
    real(r64):: nucleus_mass_number,nucleus_proton_number,gauss_nh,basis_b0,basis_beta0,&
                iteration_si, force_amsig, force_amome,force_amrho,force_gsig,force_gome,&
                force_grho,force_g2,force_g3,tmp_vcn
    integer :: basis_NF,basis_NB, basis_block_number, basis_levels_number_fermions,&
                   index,it
    real(r64),dimension(nghl,itx) :: tmp_vps,tmp_vms,tmp_vpstot,tmp_vmstot


    WFS_DIR = find_file("GCM_FILES_DIR",'HFB.wfs')
    open(u_wfs,file=trim(WFS_DIR)//outputfile%outputw,status='old')
    read(u_wfs,format1)  nucleus_name,nucleus_mass_number,nucleus_proton_number,gauss_nh,&
                         basis_NF,basis_NB,basis_block_number,basis_levels_number_fermions
    read(u_wfs,format2) basis_b0,basis_beta0,iteration_si
    read(u_wfs,format3) erwar%ea, erwar%rms, erwar%qp
    read(u_wfs,format4) force_amsig, force_amome, force_amrho
    read(u_wfs,format4) force_gsig, force_gome, force_grho, force_g2, force_g3
    read(u_wfs,format5) pairing%ga,pairing%gg,pairing%pwi
    read(u_wfs,format5) pairing%del,pairing%dec
    read(u_wfs,format5) pairing%spk
    read(u_wfs,format5) fermi%ala,fermi%tz
    read(u_wfs,format5) constraint%c1x, constraint%c2x, constraint%c3x
    read(u_wfs,format5) constraint%c1xold, constraint%c2xold, constraint%c3xold
    read(u_wfs,format5) constraint%calq1, constraint%calq2, constraint%calq3
    read(u_wfs,format6) fields%vps
    read(u_wfs,format6) fields%vms
    if (pairing%ide(1)==4) then
        read(u_wfs,*)
        read(u_wfs,format6) pairing%delq
    endif
    close(u_wfs)
    tmp_vps = reshape(fields%vps,[nghl,itx])
    tmp_vms = reshape(fields%vms,[nghl,itx])
    do index=1,nghl
        tmp_vcn = constraint%c1x * constraint%vc(index,1) +&
                  constraint%c2x * constraint%vc(index,2) + &
                  constraint%c3x * constraint%vc(index,3)
        do it=1,itx
            tmp_vpstot(index,it) = tmp_vps(index,it) + tmp_vcn
            tmp_vmstot(index,it) = tmp_vms(index,it) + tmp_vcn
        enddo
    enddo
    fields%vpstot = reshape(tmp_vpstot,[ngh,ngl,itx])
    fields%vmstot = reshape(tmp_vmstot,[ngh,ngl,itx])
end subroutine 


!--------------------------------------------------------------------------------------!
!function find_file:                                                                   ! 
!returns directory that find file_name from input environment variable(directory_name).!
!--------------------------------------------------------------------------------------!
function find_file(environment_name,file_name)
    character(500) :: find_file
    character(*) :: file_name,environment_name
    character(1000) :: path
    integer :: istart,iend,i
    logical :: isthere,last_chance
    call GETENV(trim(environment_name),path)
        path=adjustl(path)
        istart = 1
        do while (.true.)
            i = istart
            do while (path(i:i).ne.':')
                iend = i
                i = i+ 1
                if (path(i:i) == ' ') then
                    i = 0
                exit
                end if
            end do
            inquire(file=path(istart:iend)//'/'//trim(file_name),exist=isthere)
            if (isthere) then
                find_file = path(istart:iend)//'/'
                return
            else if (i == 0) then
                inquire(file='./'//file_name,exist=last_chance)
                if (last_chance) then
                    find_file = './'
                    return
                else
                    print*,'FILE NOT FOUND: ',trim(adjustl(file_name))
                    print*,'CHECK ENVIRONMENT VARIABLE: ',trim(adjustl(environment_name))
                    stop
                end if
            end if
            istart = iend+2
        end do
end function find_file

END MODULE Inout