
!! Interface between Gaussian09 CNS and DL-FIND
!! COPYRIGHT
!!
!!  Copyright 2007 Johannes Kaestner (kaestner@theochem.uni-stuttgart.de),
!!  Tom Keal (thomas.keal@stfc.ac.uk)
!!
!!  This file is part of DL-FIND.
!!
!!  DL-FIND is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as 
!!  published by the Free Software Foundation, either version 3 of the 
!!  License, or (at your option) any later version.
!!
!!  DL-FIND is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public 
!!  License along with DL-FIND.  If not, see 
!!  <http://www.gnu.org/licenses/>.
!!

module GauCNS_parameter_module
  use dlf_parameter_module
  !method number : 1 for G09 external call CNS
  !                2 for G09 and CNS separated
  integer :: method_num
  integer :: H_natom !number with H atoms
  integer :: line_len ! 58 for simple computation, 61 for freeze info, >63 for Oniom
  real(rk), PARAMETER :: bohr2A = 0.52917721, h2kcal = 627.509
  integer, dimension(:), allocatable :: gjf_pdb_map  ! atom map between gjf (with H) and pdb (no H)
  !Gaussian-CNS computation parameters definition

end module GauCNS_parameter_module

program main
  !use dlfind_module, only: dlf_nframe, dlf_nweight, dlf_nmass, dlf_nz, &
  !                  & dlf_ncons, dlf_nconn, dlf_coords, dlf_n_po_scaling
  !use basato_module ! for xa
  use dlf_parameter_module, only: rk
  implicit none
  integer :: ivar
  !add variables
  integer :: narg
  character(len=25) :: inputf

  integer :: natoms, iostatus, H_natoms
  integer :: nvar, nvar2, nspec, master
  logical :: b, alive
  character(len=50)          :: c, gjffile

  narg = iargc()

  if (narg == 1) then
    call getarg(1, inputf)

    call system('cp '//trim(inputf)//' parafile')

    INQUIRE(FILE = trim(inputf), EXIST = alive)
    IF (.not. alive) THEN
      write(*,*) 'file '//trim(inputf)//' do not exit!'
      stop
    ENDIF

    open(21, file = trim(inputf), action = 'read')
    b = .false.
  !  do while(.not. b)
  !    read(21, *) c 
  !    b = (trim(c) =='G09FILE')
  !  enddo
  !  read(21, *) gjffile
!
  !  b = .false.
  !  do while(.not. b)
  !    read(21, *) c 
  !    b = (trim(c) =='NAT')
  !  enddo
  !  read(21, *) natoms
  !  close(21)

    H_natoms = 1
    do while(.not. b)
      read(21, '(A)', iostat=iostatus) c
      if (iostatus < 0) then
        exit
      endif
  
      select case (trim(c))
        case('GAUFILE')
          read(21, *) gjffile
        case('NAT')
          read(21, *) natoms
        case('HNAT')
          read(21, *) H_natoms
      end select
  
    enddo
    close(21)


    call dlf_mpi_initialize() ! only necessary for a parallel build; 
                              ! can be present for a serial build
    call dlf_output(6,0)

    call dlf_init()

    call dl_find(3*max(natoms, H_natoms), 0, 2*max(natoms, H_natoms), 1)

    call dlf_mpi_finalize() ! only necessary for a parallel build;
                            ! can be present for a serial build

    call system('rm parafile *.temp force2 mmen2')

  else
    write(*,*) 'Just indicate the parameter filename!'
  endif
end program main


subroutine dlf_init()
  use GauCNS_parameter_module
  implicit none

  character(len=80) :: c, gjffile, pdbfile, method, mapfile
  logical :: b, alive
  integer :: num_atom, iostatus, mapping, i

  b = .false.
  open(21, file = 'parafile', action = 'read')

  H_natom = 0
  mapping = 0
  

  do while(.not. b)
    read(21, '(A)', iostat=iostatus) c
    if (iostatus < 0) then
      exit
    endif

    select case (trim(c))
      case('GAUFILE')
        read(21, *) gjffile
      case('NAT')
        read(21, *) num_atom
      case('MAPPING')
        read(21, *) mapping
      case('MAPFILE')
        read(21, *) mapfile
      case('HNAT')
        read(21, *) H_natom
      case('METHOD')
        read(21, *) method_num
    end select

  enddo
  close(21)

  allocate(gjf_pdb_map(num_atom))

  if (method_num == 0) then
    write(*,*) 'using external QM program'
    write(*,*) 'For external QM program, coord.xyz file with all atomic coordinates,'
      write(*,*) 'relaxed.dat file with all relaxed atom numbers in the coord.xyz files'
      write(*,*) 'and the script external.sh should be completed for calling external programs and transfer coordinates info'
  elseif (method_num == 1) then
    write(*,*) 'using only Gaussian '
    open(19, file = trim(gjffile), action = 'read')

    b = .false.
    do while(.not. b)
      read(19,'(A)') c 
      b = (trim(c) == '')
    enddo
  
    b = .false.
    do while(.not. b)
      read(19,'(A)') c 
      b = (trim(c) == '')
    enddo
  
    read(19,'(A)') c 
    read(19,'(A)') c 
    line_len = len_trim(c)
    close(19)
    !write(*,*) line_len
  elseif (method_num == 2) then
    write(*,*) 'using Gaussian + CNS'
    open(19, file = trim(gjffile), action = 'read')

    b = .false.
    do while(.not. b)
      read(19,'(A)') c 
      b = (trim(c) == '')
    enddo
  
    b = .false.
    do while(.not. b)
      read(19,'(A)') c 
      b = (trim(c) == '')
    enddo
  
    read(19,'(A)') c 
    read(19,'(A)') c 
    line_len = len_trim(c)
    close(19)
    !write(*,*) line_len
  else
    write(*,*) 'wrong input for method, only GAUSSIAN and GAUCNS are supported for now!'
    stop
  endif

  if (mapping == 0 ) then
    do i = 1, num_atom
      gjf_pdb_map(i) = i
    enddo
  elseif (mapping == 1) then
    INQUIRE(FILE = trim(mapfile), EXIST = alive)
    IF (.not. alive) THEN
      write(*,*) trim(mapfile)//' does not exit, or no input for mapfile'
      stop
    ENDIF
    open (19, file = trim(mapfile), action = 'read')
      do i = 1, num_atom
        read(19, *) gjf_pdb_map(i)
      enddo
    close(19)
  elseif (mapping == 2) then
    if (method_num == 0) then
      call map_xyz_cns('coord.xyz', H_natom, 'mm3.pdb', num_atom, gjf_pdb_map)
    else
      call map_g09_cns(trim(gjffile), H_natom, 'mm3.pdb', num_atom, gjf_pdb_map)
    endif

    open (19, file = 'mapfile.dat', status = 'replace')
      do i = 1, num_atom
        write(19, '(I6)') gjf_pdb_map(i)
      enddo
    close(19)
  else
    write(*,*) 'Only 0 (no map) 1 (from map file) 2 (gjf ~ pdb) availble for MAPPING option'
    stop
  endif
  INQUIRE(FILE = 'result.temp', EXIST = alive)
  if (.not. alive) then
    call lanch_external(trim(gjffile), num_atom)
  endif

end subroutine dlf_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!external subroutine (gjf pdb map)
subroutine map_g09_cns(gjffile, num_gjf, pdbfile, num_pdb, map_list)
  use dlf_parameter_module, only: rk
  use GauCNS_parameter_module
  implicit none
  character(len=*), intent(in) :: gjffile, pdbfile
  integer, intent(in) :: num_gjf, num_pdb
  
  integer, dimension(num_pdb) :: map_list

  real(rk), dimension(3, num_pdb) :: coord_cns
  real(rk), dimension(3, num_gjf) :: coord_gjf
  character(len=2), dimension(num_gjf) :: ele_gjf

  real(rk), dimension(3) :: coord_temp

  character(len = 150) :: line
  character(len = 50) :: c1, c2, c3

  logical :: b 
  integer :: i, j, pos_temp, number_0

  !read gjf file
  open(51, file = trim(gjffile), action = 'read')
  b = .false.
  do while(.not. b)
      read(51,'(A)') c1
      b = (trim(c1) == '')
  enddo
  b = .false.
  do while(.not. b)
      read(51,'(A)') c1
      b = (trim(c1) == '')
  enddo
  read(51,'(A)') c1

  if (line_len < 60) then
    do i = 1, num_gjf
      read(51, '(A)') line
      read(line(17:58),'(3F14.8)') coord_gjf(1:3,i)
      read(line(1:10), *) ele_gjf(i)
  enddo
  else
    do i = 1, num_gjf
        read(51, '(A)') line
        read(line(20:61),'(3F14.8)') coord_gjf(1:3,i)
        read(line(1:10), *) ele_gjf(i)
    enddo
  endif
  close(51)

  !read pdb file
  open(52, file = trim(pdbfile), action = 'read')
  b = .false.
  do while(.not. b)
      read(52, '(A)') c1
      b = (c1(1:4) == 'ATOM')
  enddo
  backspace(52)
  do i = 1, num_pdb 
      read(52,'(A26,4X,3F8.3,A26)') c1, coord_cns(1:3,i), c2
  enddo
  close(52)

  number_0 = 0
  do i = 1, num_pdb
      coord_temp(1:3) = coord_cns(1:3,i)
      call find_pos(coord_temp,  coord_gjf, ele_gjf, num_gjf, pos_temp)
      if (pos_temp == 0) then
          write(*,'(I5)',advance='no') i 
          number_0 = number_0 + 1
      endif
      map_list(i) = pos_temp
  enddo

  if (number_0 > 0) then
      write(*, *) number_0,' atoms are not mapped in the Gaussian gjf file, possible reason is the conformer.'
      stop
  endif
endsubroutine

!external subroutine (xyz pdb map)
subroutine map_xyz_cns(xyzfile, num_xyz, pdbfile, num_pdb, map_list)
  use dlf_parameter_module, only: rk
  use GauCNS_parameter_module
  implicit none
  character(len=*), intent(in) :: xyzfile, pdbfile
  integer, intent(in) :: num_xyz, num_pdb
  
  integer, dimension(num_pdb) :: map_list

  real(rk), dimension(3, num_pdb) :: coord_cns
  real(rk), dimension(3, num_xyz) :: coord_xyz
  character(len=2), dimension(num_xyz) :: ele_xyz

  real(rk), dimension(3) :: coord_temp

  character(len = 150) :: line
  character(len = 50) :: c1, c2, c3

  logical :: b 
  integer :: i, j, pos_temp, number_0, num_temp

  !read xyz file
  open(51, file = trim(xyzfile), action = 'read')
  read(51, *) num_temp
  if (num_temp /= num_xyz) then
    write(*,*) 'num_xyz is different with number in the .xyz file'
    stop
  endif

  read(51,'(A)') c1

  do i = 1, num_xyz
    read(51, *) ele_xyz(i), coord_xyz(1:3,i)
  enddo
  
  close(51)

  !read pdb file
  open(52, file = trim(pdbfile), action = 'read')
  b = .false.
  do while(.not. b)
      read(52, '(A)') c1
      b = (c1(1:4) == 'ATOM')
  enddo
  backspace(52)
  do i = 1, num_pdb 
      read(52,'(A26,4X,3F8.3,A26)') c1, coord_cns(1:3,i), c2
  enddo
  close(52)

  number_0 = 0
  do i = 1, num_pdb
      coord_temp(1:3) = coord_cns(1:3,i)
      call find_pos(coord_temp,  coord_xyz, ele_xyz, num_xyz, pos_temp)
      if (pos_temp == 0) then
          write(*,'(I5)',advance='no') i 
          number_0 = number_0 + 1
      endif
      map_list(i) = pos_temp
  enddo

  if (number_0 > 0) then
      write(*, *) number_0,' atoms are not mapped in the Gaussian gjf file, possible reason is the conformer.'
      stop
  endif
endsubroutine

subroutine find_pos(pos, list, ele_list, dim_list, findpos)
  use dlf_parameter_module, only: rk
  implicit none
  real(rk), dimension(3), intent(in) :: pos 
  integer, intent(in) :: dim_list
  real(rk), dimension(3, dim_list), intent(in) :: list
  character(len=2), dimension(dim_list), intent(in) :: ele_list
  integer, intent(out) :: findpos

  integer :: i 
  real(rk) :: diff

  findpos = 0
  do i = 1, dim_list
    if ( ele_list(i) == 'H-' .or. trim(ele_list(i)) == 'H') then
      cycle
    endif
    diff = ((pos(1) - list(1,i))**2 + (pos(2) - list(2,i))**2 + (pos(3) - list(3,i))**2 )**0.5
    if (diff < 0.0001 ) then
        findpos = i 
        exit
    endif
  enddo
endsubroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine lanch_external(gjffile, num_atom)
  use GauCNS_parameter_module
  use dlf_parameter_module, only: rk
  implicit none
  character(len=*), intent(in) :: gjffile
  integer, intent(in) :: num_atom
  
  real(rk), dimension(3, max(num_atom, H_natom)) :: grad_g09
  real(rk), dimension(3,num_atom) :: grad_cns
  real(rk) :: ene_g09, ene_cns, r_factor, rfree
  real(rk) :: ene_temp, x, y, z
  character(len = 80) :: c, c1
  character(len=20) :: str_num
  integer :: i
  logical :: alive

  !if ( H_natom > num_atom) then
  !  write(str_num, *) H_natom+1
  !else
  !  write(str_num, *) num_atom+1
  !endif

  write(str_num, *) max(num_atom, H_natom)+1
  !write(*,*) adjustl(trim(str_num))
  !write(*,*) 'tail -n +'//adjustl(trim(str_num))//' fort.7 >> result.temp'
  
    

  if ( method_num == 1) then
    call system('$Gauexe < '//trim(gjffile)//' > gaussian.out')
    
    !gradient to grad.temp file
    call system("grep 'extrapolated energy' gaussian.out | awk '{print $5}' > result.temp ")
    open(16, file='result.temp',action='read')
    if(EOF(16)) then
      close(16)
      call system("grep 'SCF Done:' gaussian.out | awk '{print $5}' > result.temp ")
    else
      close(16)
    endif
    
    !if (line_len < 63 ) then
    !  call system("grep 'SCF Done:' gaussian.out | awk '{print $5}' > result.temp ")
    !else
    !  call system("grep 'extrapolated energy' gaussian.out | awk '{print $5}' > result.temp ")
    !endif

    !call system('tail -n +'//adjustl(trim(str_num))//' fort.7 >> result.temp')
    !using reading file
    open(20, file='result.temp', action='read')
    read(20,*) ene_temp
    close(20)

    open(21,file='fort.7', action='read')
    open(22,file='result.temp',status='replace')
    write(22,'(F16.8)') ene_temp

    do i= 1, max(num_atom, H_natom)
      read(21,*)
    enddo

    do i=1, max(num_atom, H_natom)
      read(21,*) x, y, z
      write(22,'(3F16.8)') x, y, z
    enddo
    close(22)
    close(21)

  elseif (method_num == 0 .or. method_num == 2) then
    if ( method_num == 0) then
      !QM part
      call system('sh external.sh')
      !output a file QM_grad.dat
      INQUIRE(FILE = "QM_grad.dat", EXIST = alive)
      IF (.not. alive) THEN
          write(*,*) 'No QM_grad.dat output file from the external QM program'
          stop
      ENDIF
      call system("sed -n '1p' QM_grad.dat > energy.temp")
      call system("sed -n '2,$'p QM_grad.dat > grad.temp")

    else
      !Gaussian part
      call system('$Gauexe < '//trim(gjffile)//' > gaussian.out')

      !call system('tail -n +'//adjustl(trim(str_num))//' fort.7 > grad.temp')
      open(21,file='fort.7', action='read')
      open(22,file='grad.temp',status='replace')

      do i= 1, max(num_atom, H_natom)
        read(21,*)
      enddo

      do i= 1, max(num_atom, H_natom)
        read(21,*) x, y, z
        write(22,'(3F16.8)') x, y, z
      enddo
      close(22)
      close(21)

      !call system('mv grad.temp fort.7')
      call system("grep 'extrapolated energy' gaussian.out | awk '{print $5}' > energy.temp ")
      open(16, file='energy.temp',action='read')
      if(EOF(16)) then
        close(16)
        call system("grep 'SCF Done:' gaussian.out | awk '{print $5}' > energy.temp ")
      else
        close(16)
      endif
    
      !if (line_len < 63 ) then
      !  call system("grep 'SCF Done:' gaussian.out | awk '{print $5}' > energy.temp ")
      !else
      !  call system("grep 'extrapolated energy' gaussian.out | awk '{print $5}' > energy.temp ")
      !endif

    endif

    open(16, file = 'grad.temp', action = 'read')
    do i = 1, max(num_atom, H_natom)
      read(16, *) grad_g09(1:3, i)
    enddo
    close(16)

    open (14, file = 'energy.temp', action='read')
    read(14,*) ene_g09
    close(14)
    write(*,*) 'QM energy (a.u.): ', ene_g09
    !CNS part
    INQUIRE(FILE = "force2", EXIST = alive)
    IF (alive) THEN
        CALL SYSTEM("rm force2")
    ENDIF

    INQUIRE(FILE = "mmen2", EXIST = alive)
    IF (alive) THEN
        CALL SYSTEM("rm mmen2")
    ENDIF

    INQUIRE(FILE = "minimize.pdb", EXIST = alive)
    IF (alive) THEN
          CALL SYSTEM("rm minimize.pdb")
    ENDIF

    INQUIRE(FILE = "minimize.pdb1", EXIST = alive)
    IF (alive) THEN
          CALL SYSTEM("rm minimize.pdb1")
    ENDIF

    call system('cns < minimize.inp > minimize.out')

    open(18, file = 'mmen2', action = 'read')
    read(18, '(A)') c1
    read(18, *) ene_cns
    read(18, '(A)') c1
    read(18, '(A)') c1
    read(18, '(A)') c1
    read(18, *) r_factor, rfree
    close(18)

    write(*,*) 'XREF energy (a.u.): ', ene_cns/h2kcal
    write(*,*) 'R(factor) = ', r_factor, '   R(free) = ', rfree

    open(17, file = 'force2', action = 'read')
    read(17, '(A)') c1
    read(17, '(A)') c1
    do i = 1, num_atom
      read(17, *) grad_cns(1:3,i)
      !write(*,*) gjf_pdb_map(i)
      !grad_cns(1:3,i) = grad_cns(1:3,i)* bohr2A / h2kcal
      grad_g09(1:3, gjf_pdb_map(i)) = grad_g09(1:3, gjf_pdb_map(i)) + grad_cns(1:3,i)* bohr2A / h2kcal
    enddo
    close(17)

    !write the final gradient and energy
    open (32, file = 'result.temp', status = 'replace')
    write(32, '(F16.8)') ene_g09 + ene_cns/h2kcal

    do i = 1, max(num_atom, H_natom)
      write(32, '(3F16.8)') grad_g09(:,i)
    enddo
    write(32, *)
    close(32)

    call system('cns < bindividual.inp > bindividual.out')
    INQUIRE(FILE = "bindividual.pdb1", EXIST = alive)
    IF (alive) THEN
      call system('mv bindividual.pdb mm3.pdb')
      call system('mv bindividual.pdb1 mm3.pdb1')
    else
      call system('rm -f bindividual.pdb')
    ENDIF
    

    !do i = 1, num_atom
    !  write(32, '(3F16.8)') grad_cns(:,i) + grad_g09(:,i)
    !enddo
    !if (num_atom < H_natom) then
    !  do i = num_atom+1, H_natom
    !    write(32, '(3F16.8)') grad_g09(:,i)
    !  enddo
    !endif
    !write(32, *)
    !close(32)

  else
    write(*,*) 'wrong input for method, 0 (external+CNS), 1(Gaussian), 2(Gaussian + CNS) are supported for now!'
    stop
  endif
end subroutine
! **********************************************************************
! subroutines that have to be provided to dl_find from outside
! **********************************************************************

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_params(nvar,nvar2,nspec,coords,coords2,spec,ierr, &
    tolerance,printl,maxcycle,maxene,tatoms,icoord, &
    iopt,iline,maxstep,scalestep,lbfgs_mem,nimage,nebk,dump,restart,&
    nz,ncons,nconn,update,maxupd,delta,soft,inithessian,carthessian,tsrel, &
    maxrot,tolrot,nframe,nmass,nweight,timestep,fric0,fricfac,fricp, &
    imultistate, state_i,state_j,pf_c1,pf_c2,gp_c3,gp_c4,ln_t1,ln_t2, &
    printf,tolerance_e,distort,massweight,minstep,maxdump,task,temperature, &
    po_pop_size,po_radius,po_contraction,po_tolerance_r,po_tolerance_g, &
    po_distribution,po_maxcycle,po_init_pop_size,po_reset,po_mutation_rate, &
    po_death_rate,po_scalefac,po_nsave,ntasks,tdlf_farm,n_po_scaling, &
    neb_climb_test, neb_freeze_test, &
    nzero, coupled_states, qtsflag,&
    imicroiter, maxmicrocycle, micro_esp_fit)
  use dlf_parameter_module, only: rk
  use GauCNS_parameter_module
  !use vib_pot
  implicit none
  integer   ,intent(in)      :: nvar 
  integer   ,intent(in)      :: nvar2
  integer   ,intent(in)      :: nspec
  real(rk)  ,intent(inout)   :: coords(nvar) ! start coordinates
  real(rk)  ,intent(inout)   :: coords2(max(nvar2,1)) ! a real array that can be used
                                               ! depending on the calculation
                                               ! e.g. a second set of coordinates
  integer   ,intent(inout)   :: spec(nspec)  ! specifications like fragment or frozen
  integer   ,intent(out)     :: ierr
  real(rk)  ,intent(inout)   :: tolerance
  real(rk)  ,intent(inout)   :: tolerance_e
  integer   ,intent(inout)   :: printl
  integer   ,intent(inout)   :: maxcycle
  integer   ,intent(inout)   :: maxene
  integer   ,intent(inout)   :: tatoms
  integer   ,intent(inout)   :: icoord
  integer   ,intent(inout)   :: iopt
  integer   ,intent(inout)   :: iline
  real(rk)  ,intent(inout)   :: maxstep
  real(rk)  ,intent(inout)   :: scalestep
  integer   ,intent(inout)   :: lbfgs_mem
  integer   ,intent(inout)   :: nimage
  real(rk)  ,intent(inout)   :: nebk
  integer   ,intent(inout)   :: dump
  integer   ,intent(inout)   :: restart
  integer   ,intent(inout)   :: nz
  integer   ,intent(inout)   :: ncons
  integer   ,intent(inout)   :: nconn
  integer   ,intent(inout)   :: update
  integer   ,intent(inout)   :: maxupd
  real(rk)  ,intent(inout)   :: delta
  real(rk)  ,intent(inout)   :: soft
  integer   ,intent(inout)   :: inithessian
  integer   ,intent(inout)   :: carthessian
  integer   ,intent(inout)   :: tsrel
  integer   ,intent(inout)   :: maxrot
  real(rk)  ,intent(inout)   :: tolrot
  integer   ,intent(inout)   :: nframe
  integer   ,intent(inout)   :: nmass
  integer   ,intent(inout)   :: nweight
  real(rk)  ,intent(inout)   :: timestep
  real(rk)  ,intent(inout)   :: fric0
  real(rk)  ,intent(inout)   :: fricfac
  real(rk)  ,intent(inout)   :: fricp
  integer   ,intent(inout)   :: imultistate
  integer   ,intent(inout)   :: state_i
  integer   ,intent(inout)   :: state_j
  real(rk)  ,intent(inout)   :: pf_c1  
  real(rk)  ,intent(inout)   :: pf_c2  
  real(rk)  ,intent(inout)   :: gp_c3  
  real(rk)  ,intent(inout)   :: gp_c4
  real(rk)  ,intent(inout)   :: ln_t1  
  real(rk)  ,intent(inout)   :: ln_t2  
  integer   ,intent(inout)   :: printf
  real(rk)  ,intent(inout)   :: distort
  integer   ,intent(inout)   :: massweight
  real(rk)  ,intent(inout)   :: minstep
  integer   ,intent(inout)   :: maxdump
  integer   ,intent(inout)   :: task
  real(rk)  ,intent(inout)   :: temperature
  integer   ,intent(inout)   :: po_pop_size
  real(rk)  ,intent(inout)   :: po_radius
  real(rk)  ,intent(inout)   :: po_contraction
  real(rk)  ,intent(inout)   :: po_tolerance_r
  real(rk)  ,intent(inout)   :: po_tolerance_g
  integer   ,intent(inout)   :: po_distribution
  integer   ,intent(inout)   :: po_maxcycle
  integer   ,intent(inout)   :: po_init_pop_size
  integer   ,intent(inout)   :: po_reset
  real(rk)  ,intent(inout)   :: po_mutation_rate
  real(rk)  ,intent(inout)   :: po_death_rate
  real(rk)  ,intent(inout)   :: po_scalefac
  integer   ,intent(inout)   :: po_nsave
  integer   ,intent(inout)   :: ntasks
  integer   ,intent(inout)   :: tdlf_farm
  integer   ,intent(inout)   :: n_po_scaling
  real(rk)  ,intent(inout)   :: neb_climb_test
  real(rk)  ,intent(inout)   :: neb_freeze_test
  integer   ,intent(inout)   :: nzero
  integer   ,intent(inout)   :: coupled_states
  integer   ,intent(inout)   :: qtsflag
  integer   ,intent(inout)   :: imicroiter
  integer   ,intent(inout)   :: maxmicrocycle
  integer   ,intent(inout)   :: micro_esp_fit
  ! local variables
  real(rk)                   :: svar
  integer                    :: iat,jat, i, atomic, j, num_relaxed
  logical                    :: b, alive
  integer                    :: iostatus, num_atom, tolereance_int,resnum, coordfile, num_temp
  character(len=80)          :: c, c1, c2, c3, c4, giffile,resfile
  character(len=1)           :: layer
  integer, dimension(:), allocatable :: res_number
  real(rk), dimension(:,:), allocatable :: coord_pdb
  interface
    subroutine read_rand(arr)
      use dlf_parameter_module, only: rk
      use dlf_global, only : glob
      real(rk)  :: arr(:)
    end subroutine read_rand
  end interface
! **********************************************************************
!read the parameters from para.dat file (defin also the default values)


!default values
  ierr=0
  tsrel=1

  tolerance=4.5D-5 ! negative: default settings
  printl=4
  printf=4
  maxcycle=100 !200
  maxene=100000

  tolrot=1.D2
!  tolrot=0.1D0
!  maxrot=100 !was 100

  task=0 !1011

  distort=0.D0 !0.4 
  tatoms=0
  icoord=0 !0 cartesian coord !210 Dimer !190 qts search !120 NEB frozen endpoint
  massweight=0

! TO DO (urgent): better dtau for endpoints (interpolate energy) when reading dist

  imicroiter=1 ! means: use microiterative optimization

  iopt=3 ! 20 !3 or 20 later change to 12
  temperature = 0.D0! K ! T_c for the MB-potential for hydrogen is 2206.668 K

  iline=0
  maxstep=0.1D0
  scalestep=1.0D0
  lbfgs_mem=100
  nimage=39 !*k-k+1
  nebk=1.0D0 ! for QTS calculations with variable tau, nebk transports the parameter alpha (0=equidist. in tau)
  qtsflag=0 ! 1: Tunneling splittings, 11/10: read image hessians

  !"accurate" etunnel for 0.4 Tc: 0.3117823831234522

  ! Hessian
  delta=1.D-2
  soft=-6.D-4
  update=2
  maxupd=0
  inithessian=0
  minstep=1.D0 **2 ! 1.D-5
  minstep=1.D-5

!!$  ! variables for exchange by sed
!!$  iopt=IOPT_VAR
!!$  nimage=NIMAGE_VAR
!!$  temperature = TFAC ! K ! T_c for the MB-potential for hydrogen is 2206.668 K
!!$  nebk= NEBK_VAR


  ! damped dynamics
  fric0=0.1D0
  fricfac=1.0D0
  fricp=0.1D0

  dump=0
  restart=0

  ! Parallel optimization

  po_pop_size=25
  po_radius=0.5D0
  po_init_pop_size=50
  po_contraction=0.95D0
  po_tolerance_r=1.0D-8
  po_tolerance_g=1.0D-6
  po_distribution=3
  po_maxcycle=100000
  po_reset=500
  po_mutation_rate=0.15D0
  po_death_rate=0.5D0
  po_scalefac=10.0D0
  po_nsave=10
  n_po_scaling=0 ! meaning that the base radii values for the sampling and tolerance 
                 ! are used for all components of the working coordinate vector.
                 ! Remember to change the second arg in the call to dl_find if a 
                 ! non-zero value for n_po_scaling is desired, and also to add the 
                 ! necessary values to the coords2 array...

  ! Taskfarming 
  ! (the two lines below, with any values assigned, 
  ! may be safely left in place for a serial build) 
  ntasks = 1
  tdlf_farm = 1

  tatoms=1

  nmass = 0
  nframe = 0
  nweight = 0

!user-defined from input file

  resnum = 0
  coordfile = 0
  imicroiter = 0
  !H_natom = 0
  b = .false.
  open(22, file= 'parafile', action='read')
  do while(.not. b)
    read(22, *, iostat=iostatus) c
    if (iostatus < 0) then
      exit
    endif

    select case (trim(c))
    case('GAUFILE')
      read(22, *) giffile
    case('NAT')
      read(22, *) num_atom
    !case('HNAT')
    !  read(22, *) H_natom
    case('RESNUM')
      read(22, *) resnum
    case('RESFILE')
      read(22, *) resfile
    case('TOLERANCE')
      read(22, *) tolerance, tolerance_e
    case('MAXCYCLE')
      read(22, *) maxcycle
    case('NCON')
      read(22, *) ncons, nconn
    case('PRINT')
      read(22, *) printl, printf
    case('ICOORD')
      read(22, *) icoord
    case('IMULTISTATE')
      read(22, *) imultistate
    case('IOPT')
      read(22, *) iopt
    case('ILINE')
      read(22, *) iline
    case('INITHESSIAN')
      read(22, *) inithessian
    case('UPDATE')
      read(22, *) update
    case('IMICROITER')
      read(22, *) imicroiter
    case('LBFGS_MEM')
      read(22, *) lbfgs_mem
    case('MAXENE')
      read(22, *) maxene
    case('COORDFILE')
      read(22, *) coordfile
    case('MAXSTEP')
      read(22, *) maxstep
    end select

  enddo
  close(22)
 
  print*,"External sizes:",nvar,nvar2,nspec

! res number
  allocate(res_number(max(num_atom, H_natom)))

  if (resnum == 0 ) then
    do i = 1, max(num_atom, H_natom)
      res_number(i) = 1
    enddo
  elseif(resnum == 1) then
    INQUIRE(FILE = trim(resfile), EXIST = alive)
    if (.not. alive) then
      write(*,*) trim(resfile)//' does not exit, or no input for resfile'
      stop
    endif

    open(23, file = trim(resfile), action = 'read')
    do i = 1, max(num_atom, H_natom)
      read(23, *) res_number(i)
    enddo
    close(23)
  endif
!*********************************************************
!!!!!!!!!!!!!!!!!!!!!
! read coordinates
if (method_num == 0) then
! other QM program
  ! read from coord.xyz file
  INQUIRE(FILE = 'coord.xyz', EXIST = alive)
  if (.not. alive) then
    write(*,*) 'coord.xyz file does not exit!'
    stop
  endif
  open(23, file = 'coord.xyz', action = 'read')
  read(23, *) num_temp
  if (3*num_temp /= nvar) then
    write(*,*) 'Number of atoms in coord.xyz is different with parameter definition'
    close(23)
    stop
  endif

  read(23,'(A)') c 
  do i = 1, max(num_atom, H_natom)
    read(23, *) c, coords((i-1)*3+1 : 3*i)
  enddo

  close(23)

  !read relaxed atom info
  spec(:) = -1
  INQUIRE(FILE = 'relaxed.dat', EXIST = alive)
  if (.not. alive) then
    write(*,*) 'relaxed.dat file does not exit!'
    stop
  endif
  !num_relaxed = 0
  open(24, file = 'relaxed.dat', action = 'read')
  read(24, *) num_relaxed
  do i = 1, num_relaxed
    read(24, *) num_temp
    spec(num_temp) = res_number(num_temp)
  enddo
 ! do while(.not. b)
 !   read(24,'(A)') c
 !   b = (trim(c) == '')
 !   if (.not. b ) then
 !     num_relaxed = num_relaxed +1 
 !     read(c, *) num_temp
 !     spec(num_temp) = res_number(num_temp)
 !   endif
 ! enddo
  write(*,*) 'There are ', num_relaxed, ' atoms relaxed during the refinement!'
  close(24)
!**********************************************************

elseif (method_num == 1 .or. method_num == 2) then
! Gaussian program
  open(23, file = trim(giffile), action = 'read')
  b = .false.
  do while(.not. b)
    read(23,'(A)') c 
    b = (trim(c) == '')
  enddo

  b = .false.
  do while(.not. b)
    read(23,'(A)') c 
    b = (trim(c) == '')
  enddo
  read(23,'(A)') c 

  spec(:) = -1
  if (line_len < 60 ) then
    do i = 1, max(num_atom, H_natom)
      read(23,'(A)') c 
      read(c(17:58),'(3F14.8)') coords((i-1)*3+1 : 3*i)
    enddo
    
  elseif (line_len < 63 ) then
    do i = 1, max(num_atom, H_natom)
      read(23,'(A)') c 
      !residue number treatement for PDB case
      read(c(18:61),'(I2,3F14.8)') spec(i), coords((i-1)*3+1 : 3*i)
      !if (spec(i) == 0 .and. i .LE. num_atom) then
      if (spec(i) == 0 ) then
        spec(i) = res_number(i)
      !elseif (spec(i) == 0 .and. i .GT. num_atom) then
      !  spec(i) = 999
      endif
    enddo
  else
    if (imicroiter ==0) then
      do i = 1, max(num_atom, H_natom)
        read(23,'(A)') c 
        !residue number treatement for PDB case
        read(c(18:61),'(I2,3F14.8)') spec(i), coords((i-1)*3+1 : 3*i)
        if (spec(i) == 0 ) then
          spec(i) = res_number(i)
        endif
      enddo
    else
      do i = 1, max(num_atom, H_natom)
        read(23,'(A)') c 
        !residue number treatement for PDB case
        read(c(18:63),'(I2,3F14.8,1X,A1)') spec(i), coords((i-1)*3+1 : 3*i), layer
        if (spec(i) == 0 ) then
          spec(i) = res_number(i)
        endif
        if (layer == 'M') then
          spec(i+ max(num_atom, H_natom)) = 1
        else
          spec(i+ max(num_atom, H_natom)) = 0
        endif
      enddo
    endif 
  endif
  close(23)

  if (coordfile == 1) then
    allocate(coord_pdb(3,num_atom))

    open(28, file = 'mm3.pdb', action = 'read')
    b = .false.
    do while(.not. b)
      READ(28, '(A)') c
      b = ( c(1:4) == 'ATOM')
    enddo
    BACKSPACE(28)

    do i = 1, num_atom
      read(28,'(A26,4X,3F8.3,A26)') c, coord_pdb(1:3,i), c1
      j = gjf_pdb_map(i)
      coords((j-1)*3+1 : 3*j) = coord_pdb(1:3,i)
    enddo
    close(28)
    deallocate(coord_pdb)
  endif
endif
  coords2(:)=0.D0
  !open(23, file='fort.7', action='read')
  !  do i = 1, num_atom
  !    read(23, *) atomic, coords((i-1)*3+1 : 3*i)
  !  enddo
  !close(23)
!
  !coords2(:)=0.D0
  !spec(:) = 1
  deallocate(res_number)

end subroutine dlf_get_params

subroutine test_update
  use dlf_parameter_module, only: rk
  use dlf_allocate, only: allocate,deallocate
  use dlf_global, only: glob,printl
  implicit none
  integer(4) :: varperimage,nimage,iimage,ivar4
  real(rk), allocatable:: coords(:,:),grad(:,:) ! varperimage,nimage
  real(rk), allocatable:: hess(:,:,:) ! varperimage,varperimage,nimage
  real(rk), allocatable:: fhess(:,:) ! 2*varperimage*nimage,2*varperimage*nimage
  real(rk), allocatable:: eigval(:),eigvec(:,:)
  real(rk), allocatable:: vec0(:)
  real(rk), allocatable:: determinant(:)
  real(rk), allocatable:: tmphess(:,:)
  real(rk) :: svar
  integer :: ivar,jvar,vp8,target_image,step,lastimage,jimage,turnimg
  logical :: havehessian,fracrecalc

  open(unit=102,file="grad_coor.bin",form="unformatted")
  read(102) varperimage,nimage
  print*,"varperimage,nimage",varperimage,nimage
  vp8=varperimage

!!$  call allocate(coords,int(varperimage,kind=8),int(nimage,kind=8))
!!$  call allocate(grad,int(varperimage,kind=8),int(nimage,kind=8))
!!$  call allocate(determinant,int(nimage,kind=8))

  do iimage=1,nimage
    read(102) ivar4
    print*,"Reading coords/grad of image",ivar4
    read(102) grad(:,iimage)
    read(102) coords(:,iimage)
  end do
  close(102)
  print*,"Coords sucessfully read"

  ! print coords and grad
  iimage=2
  print*,"Coords and grad for image ",iimage
  do ivar=1,vp8
    write(6,"(i6,1x,2es18.9)") &
        ivar,coords(ivar,iimage),grad(ivar,iimage)
  end do

  open(unit=101,file="hessian.bin",form="unformatted")
  ! Hessian in mass-weighted coordinates on the diagnal blocks - everything else should be zero
  read(101) iimage,ivar4!neb%varperimage,neb%nimage
  if(iimage/=varperimage.or.ivar4/=nimage) then
    print*,"Dimensions read",iimage,ivar4
    call dlf_fail("ERROR: wrong dimensions in hessian.bin!")
  end if
  ivar=2*varperimage*nimage
  print*,"File Hessian size",ivar
  call allocate(fhess,ivar,ivar)
  read(101) fhess
  close(101)
  print*,"Hessian sucessfully read"

  ! map hessian to different array:
!!$  call allocate(hess,int(varperimage,kind=8),int(varperimage,kind=8),int(nimage,kind=8))
  do iimage=1,nimage
    print*,"Image",iimage,"Hessian positions",(iimage-1)*varperimage+1,iimage*varperimage
    hess(:,:,iimage)=fhess((iimage-1)*varperimage+1:iimage*varperimage,(iimage-1)*varperimage+1:iimage*varperimage)
  end do
  call deallocate(fhess)

  call allocate(eigval,vp8)
  call allocate(eigvec,vp8,vp8)
  ! now we have all we need

  print*,"# Distance from previous image"
  do iimage=2,nimage
    print*,iimage,sqrt(sum( (coords(:,iimage)-coords(:,iimage-1))**2))
  end do

  print*,"# Determinant of Hessian"
  do iimage=1,nimage
    do ivar=1,vp8
      do jvar=ivar+1,vp8
        if(abs(hess(ivar,jvar,iimage)-hess(jvar,ivar,iimage))>1.D-20) &
            print*,"Unsymmetric:",ivar,jvar,iimage,hess(ivar,jvar,iimage),hess(jvar,ivar,iimage)
      end do
    end do
    call dlf_matrix_diagonalise(vp8,hess(:,:,iimage),eigval,eigvec)
    !write(6,"(i6,1x,9es12.3)") iimage,product(eigval(8:vp8)),eigval(1:8)
    do ivar=1,6
      eigval(minloc(abs(eigval)))=1.D0
    end do
    determinant(iimage)=product(eigval)
    write(6,"(i6,1x,9es12.3)") iimage,product(eigval),eigval(1:8)
  end do
  
  print*,"maxval(hess(:,:,nimage-1)-hess(:,:,nimage))",maxval(hess(:,:,nimage-1)-hess(:,:,nimage))

!!$  ! Richtungsableitung
!!$  call allocate(vec0,vp8)!int(varperimage,kind=8))
!!$  do iimage=2,nimage-1
!!$    ! eigval is Vector along which the derivative is taken
!!$    eigval=coords(:,iimage+1)-coords(:,iimage-1)
!!$    svar=sqrt(sum(eigval**2))
!!$    eigval=eigval/svar
!!$    vec0=matmul(hess(:,:,iimage),eigval)
!!$    do ivar=1,vp8
!!$      write(6,"(2i6,1x,2es18.9,1x,f10.5)") &
!!$          iimage,ivar,(grad(ivar,iimage+1)-grad(ivar,iimage-1))/svar,vec0(ivar),&
!!$          vec0(ivar)/((grad(ivar,iimage+1)-grad(ivar,iimage-1))/svar)
!!$    end do
!!$  end do

  !
  ! now test updates
  !
  call allocate(tmphess,vp8,vp8)
  havehessian=.true.
  fracrecalc=.false.
  printl=2
  glob%maxupd=30000

  target_image=1 !nimage
  ! update hessians to the one of the first image
  print*,"Updating Hessians to that of image",target_image
  print*,"Sum-of-squares difference"
  do iimage=1,nimage
    tmphess(:,:)=hess(:,:,iimage)
    call dlf_hessian_update(vp8, &
        coords(:,target_image),coords(:,iimage),&
        grad(:,target_image),grad(:,iimage), &
        tmphess, havehessian, fracrecalc)
    if(.not.havehessian) then
      print*,"Problem with hessian update, image",iimage
      havehessian=.true.
    end if
    print*,iimage,sum( (tmphess-hess(:,:,target_image))**2),&
        sum( (hess(:,:,iimage)-hess(:,:,target_image))**2)
  end do

  print*,"Minstep",glob%minstep
  print*,"Updating Hessians to that of image",target_image
  print*,"Determinant"
  open(file="determinant",unit=10)
  do iimage=1,nimage
    tmphess(:,:)=hess(:,:,iimage)
    call dlf_hessian_update(vp8, &
        coords(:,target_image),coords(:,iimage),&
        grad(:,target_image),grad(:,iimage), &
        tmphess, havehessian, fracrecalc)
    if(.not.havehessian) then
      print*,"Problem with hessian update, image",iimage
      havehessian=.true.
    end if
    call dlf_matrix_diagonalise(vp8,tmphess,eigval,eigvec)
    !write(6,"(i6,1x,9es12.3)") iimage,product(eigval(8:vp8)),eigval(1:8)
    do ivar=1,6
      eigval(minloc(abs(eigval)))=1.D0
    end do

    print*,iimage,product(eigval),determinant(iimage),determinant(target_image)
    write(10,*) iimage,product(eigval),determinant(iimage),determinant(target_image)
  end do
  close(10)

  print*,"Updating Hessians to that of image",target_image
  print*,"Determinant - incremental"
  open(file="determinant_incr",unit=10)
  do iimage=1,nimage
    tmphess(:,:)=hess(:,:,iimage)
    step=1
    if(iimage>target_image) step=-1
    lastimage=iimage
    do jimage=iimage+step,target_image,step
      !print*,"updating",lastimage," to ",jimage
      call dlf_hessian_update(vp8, &
          coords(:,jimage),coords(:,lastimage),&
          grad(:,jimage),grad(:,lastimage), &
          tmphess, havehessian, fracrecalc)

      if(.not.havehessian) then
        print*,"Problem with hessian update, image",iimage
        havehessian=.true.
      end if
      lastimage=jimage
    end do
    call dlf_matrix_diagonalise(vp8,tmphess,eigval,eigvec)
    !write(6,"(i6,1x,9es12.3)") iimage,product(eigval(8:vp8)),eigval(1:8)
    do ivar=1,6
      eigval(minloc(abs(eigval)))=1.D0
    end do

    print*,iimage,product(eigval),determinant(iimage),determinant(target_image)
    write(10,*) iimage,product(eigval),determinant(iimage),determinant(target_image)
  end do
  close(10)

  print*,"Updating Hessians to that of image",target_image
  print*,"Determinant - incremental turning around"
  turnimg=20
  open(file="determinant_turn",unit=10)
  do iimage=1,nimage
    tmphess(:,:)=hess(:,:,iimage)
    lastimage=iimage
    ! first upwards to turnimg
    do jimage=iimage+1,turnimg
      print*,"updating",lastimage," to ",jimage
      call dlf_hessian_update(vp8, &
          coords(:,jimage),coords(:,lastimage),&
          grad(:,jimage),grad(:,lastimage), &
          tmphess, havehessian, fracrecalc)

      if(.not.havehessian) then
        print*,"Problem with hessian update, image",iimage
        havehessian=.true.
      end if
      lastimage=jimage
    end do
    step=1
    if(lastimage>target_image) step=-1
    do jimage=lastimage+step,target_image,step
      print*,"updating",lastimage," to ",jimage
      call dlf_hessian_update(vp8, &
          coords(:,jimage),coords(:,lastimage),&
          grad(:,jimage),grad(:,lastimage), &
          tmphess, havehessian, fracrecalc)

      if(.not.havehessian) then
        print*,"Problem with hessian update, image",iimage
        havehessian=.true.
      end if
      lastimage=jimage
    end do
    call dlf_matrix_diagonalise(vp8,tmphess,eigval,eigvec)
    !write(6,"(i6,1x,9es12.3)") iimage,product(eigval(8:vp8)),eigval(1:8)
    do ivar=1,6
      eigval(minloc(abs(eigval)))=1.D0
    end do

    print*,iimage,product(eigval),determinant(iimage),determinant(target_image)
    write(10,*) iimage,product(eigval),determinant(iimage),determinant(target_image)
  end do
  close(10)

!!$  do ivar=1,vp8
!!$    WRITE(*,'(33f10.5)') hess(ivar,:,1)*1.D6
!!$  end do
!!$
!!$  do ivar=1,vp8
!!$    write(6,"(i6,1x,es18.9)") &
!!$            ivar,eigval(ivar)
!!$  end do

  call deallocate(coords)
  call deallocate(grad)
  
  call dlf_fail("stop in test_update")
end subroutine test_update

! just a test routine
subroutine test_ene
  use dlf_parameter_module, only: rk
  implicit none
  integer :: ivar,ivar2,status
  integer :: samples
  real(rk) :: halfSamples
  real(rk) :: coords(3),grad(3),hess(3,3),ene
  coords(:)=0.D0
!  open(file="energy",unit=13)
  open(file="energy-2d.dat",unit=13)
  samples = 100
  halfSamples = dble(samples) * 0.25D0
  do ivar2=1,samples
    do ivar=1,samples
      coords(1)=dble(ivar-samples/2)/halfSamples
      coords(2)=dble(ivar2-samples/2)/halfSamples
      call dlf_get_gradient(3,coords,ene,grad,1,-1,status)
      call dlf_get_hessian(3,coords,hess,status)
     write(13,*) coords(1),coords(2),ene,grad(1),grad(2),hess(1,1),hess(2,2),hess(2,1)
!     write(13,*) coords(1),coords(2),ene,grad(1)
    end do
    write(13,*) ""
  end do
  close(13)
  call dlf_fail("stop in test_ene")
end subroutine test_ene

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_gradient(nvar,coords,energy,gradient,iimage,kiter,status)
  use GauCNS_parameter_module
  use dlf_parameter_module, only: rk
  !use driver_module
  !use driver_parameter_module
  !use vib_pot
  implicit none
  integer   ,intent(in)    :: nvar
  real(rk)  ,intent(in)    :: coords(nvar)
  real(rk)  ,intent(out)   :: energy
  real(rk)  ,intent(out)   :: gradient(nvar)
  integer   ,intent(in)    :: iimage
  integer   ,intent(in)    :: kiter
  integer   ,intent(out)   :: status
  !
  integer  :: i
  character(len=80) :: c

  real(rk) :: acoords(3,nvar/3)
  real(rk) :: agrad(3,nvar/3)
  real(rk) :: r
  integer  :: nat,iat,jat
  ! additional variables non-cont. diff MEP potential
  real(rk) :: t ! variation of the depth of the flater saddle point

! **********************************************************************
!  get the gradient
  status=1
  open(23, file = 'result.temp', action = 'read')
  read(23, *) energy
  do i = 1, nvar/3
    read(23, *) gradient(3*(i-1)+1 : 3*i)
  enddo

  close(23)
!  open(23, file='fort.7', action='read')
!    do i = 1, nvar/3
!      read(23, '(A)') c
!    enddo
!
!    do i = 1, nvar/3
!      read(23, *) gradient(3*(i-1)+1 : 3*i)
!    enddo
!  close(23)
!
!! get the energy
!  call system("grep 'SCF Done:' gaussian.out | awk '{print $5}' > energy.temp ")
!  open (27, file = 'energy.temp', action='read')
!    read(27,*) energy
!  close(27)

  status=0
end subroutine dlf_get_gradient

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_hessian(nvar,coords,hessian,status)
  !  get the hessian at a given geometry
  use dlf_parameter_module, only: rk
  !use dlf_parameter_module
  !use driver_module
  !use driver_parameter_module
  !use vib_pot
  implicit none
  integer   ,intent(in)    :: nvar
  real(rk)  ,intent(in)    :: coords(nvar)
  real(rk)  ,intent(out)   :: hessian(nvar,nvar)
  integer   ,intent(out)   :: status
  real(rk) :: acoords(3,nvar/3),r,svar,svar2
  integer  :: posi,posj,iat,jat,m,n
  ! variables for Mueller-Brown potential
  real(rk) :: x,y
  integer  :: icount
  ! variables non-cont. diff. potential
  real(rk) :: t
! **********************************************************************
  hessian(:,:)=0.D0
  status=1

end subroutine dlf_get_hessian

! initialize parameters for the test potentials
!subroutine driver_init
!  use driver_parameter_module
!  implicit none
!  real(rk) :: ebarr_
!  integer :: icount
!  ! assign parameters for MB potential
!  ebarr_=0.5D0
!  acappar(1)=-200.D0*ebarr_/106.D0
!  acappar(2)=-100.D0*ebarr_/106.D0
!  acappar(3)=-170.D0*ebarr_/106.D0
!  acappar(4)=  15.D0*ebarr_/106.D0
!  apar(1)=-1.D0
!  apar(2)=-1.D0
!  apar(3)=-6.5D0
!  apar(4)=0.7D0
!  bpar(1)=0.D0
!  bpar(2)=0.D0
!  bpar(3)=11.D0
!  bpar(4)=0.6D0
!  cpar(1)=-10.D0
!  cpar(2)=-10.D0
!  cpar(3)=-6.5D0
!  cpar(4)=0.7D0
!  x0par(1)=1.D0
!  x0par(2)=0.D0
!  x0par(3)=-0.5D0
!  x0par(4)=-1.D0
!  y0par(1)=0.D0
!  y0par(2)=0.5D0
!  y0par(3)=1.5D0
!  y0par(4)=1.D0
!  ! parameters for polymul
!  dpar(1)=0.D0
!  do icount=2,num_dim
!    dpar(icount)=1.D-4+(dble(icount-2)/dble(num_dim-2))**2*0.415D0
!  end do
!end subroutine driver_init

! classical flux of the Eckart potential 
! returns log(flux)
!subroutine driver_clrate(beta_hbar,flux)
!  use driver_parameter_module
!  implicit none
!  real(rk),intent(in) :: beta_hbar
!  real(rk),intent(out):: flux
!  real(rk) :: vmax,pi
!  pi=4.D0*atan(1.D0)
!  vmax=(Va+Vb)**2/(4.D0*Vb)
!  print*,"V_max=",vmax
!  flux=-beta_hbar*vmax-log(2.D0*pi*beta_hbar)
!  
!end subroutine driver_clrate

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_put_coords(nvar,mode,energy,coords,iam)
  !use dlf_parameter_module
  use GauCNS_parameter_module
  use dlf_parameter_module, only: rk
  implicit none
  integer   ,intent(in)    :: nvar
  integer   ,intent(in)    :: mode
  integer   ,intent(in)    :: iam
  real(rk)  ,intent(in)    :: energy
  real(rk)  ,intent(in)    :: coords(nvar)
  integer                  :: iat, iostatus, natoms
! **********************************************************************

  logical :: b
  character(len=80) :: c, giffile
! output the 
! Only do this writing of files if I am the rank-zero processor
  if (iam /= 0) return

  b = .false.
  open(21, file= 'parafile', action='read')
  do while(.not. b)
    read(21, *, iostat=iostatus) c
    if (iostatus < 0) then
      exit
    endif

    select case (trim(c))
      case('GAUFILE')
        read(21, *) giffile
      case('NAT')
        read(21, *) natoms
    end select

  enddo
  close(21)


  !update the coordinates of the Gaussian gjf file
  if(mod(nvar,3)==0) then
    !assume coords are atoms
    if ( method_num == 1) then 
      call update_gjf(nvar/3, trim(giffile), coords)
    elseif (method_num == 2) then
      call update_gjf(nvar/3, trim(giffile), coords)
      call update_pdb(natoms, 'mm3.pdb', coords)
    elseif (method_num == 0) then
      call update_pdb(natoms, 'mm3.pdb', coords)
      call update_xyz(nvar/3, 'coord.xyz', coords)
    endif
  else
    !print*,"Coords in put_coords: ",coords
  end if

  call lanch_external(trim(giffile),natoms)
end subroutine dlf_put_coords

!called by dlf_put_coords
subroutine update_gjf(num_atom, gjffile, coords)
  !use driver_parameter_module
  use GauCNS_parameter_module
  use dlf_parameter_module, only: rk
  implicit none
  integer, intent(in) :: num_atom
  character(len=*), intent(in) :: gjffile
  real(rk), intent(in) ::  coords(3*num_atom)

  CHARACTER(len=120) :: line, chkfile
  CHARACTER(len=120) :: c1
  CHARACTER(len=18), DIMENSION(num_atom):: atomes_gjf
  CHARACTER(len=15), DIMENSION(num_atom):: atomes_gjf15
  CHARACTER(len=150), DIMENSION(num_atom*2):: gjfend
  CHARACTER(len=1), DIMENSION(num_atom) :: layer
  CHARACTER(len=150), DIMENSION(10) :: gjfheader
  CHARACTER(len=15), DIMENSION(num_atom) :: linend
  real(rk), DIMENSION(3,num_atom):: coord_gau
  integer :: num_head, i, num_end, isend
  logical :: b, isguess

  isguess = .false.
  !read old gjffile
  num_head = 0d0
  OPEN(24, FILE = TRIM(gjffile), ACTION = 'read')
  b= .FALSE.
  DO WHILE(.NOT. b)
        num_head = num_head +1d0
        READ(24, '(A)') gjfheader(num_head)
        if(index(gjfheader(num_head), 'guess=read') /= 0) then
          isguess = .true.
        elseif (index(gjfheader(num_head),'%chk=') /= 0 ) then
          c1 = gjfheader(num_head)
          chkfile = c1(6:len_trim(c1))
        endif
        b = (TRIM(gjfheader(num_head)) == '')
  ENDDO

  b= .FALSE.
  DO WHILE(.NOT. b)
        num_head = num_head +1d0
        READ(24, '(A)') gjfheader(num_head)
        b = (TRIM(gjfheader(num_head)) == '')
  ENDDO
  num_head = num_head +1d0
  READ(24, '(A)') gjfheader(num_head)

  !read atome block (with no PDB information)
  !READ(24,'(A)') line
  !BACKTRACE(24)
  !line_len = len_trim(line) 
  if (line_len < 60 ) then
    DO i=1, num_atom
        READ(24,'(A)') line
        READ(line, '(1X,A15,3F14.8)') atomes_gjf15(i), coord_gau(1:3,i)
        linend(i) = ''
    ENDDO
  elseif (line_len < 63 ) then
    DO i=1, num_atom
        READ(24,'(A)') line
        READ(line, '(1X,A18,3F14.8)') atomes_gjf(i), coord_gau(1:3,i)
        linend(i) = ''
    ENDDO
  else
    DO i=1, num_atom
      READ(24,'(A)') line
      READ(line, '(1X,A18,3F14.8,1X,A1)') atomes_gjf(i), coord_gau(1:3,i), layer(i)
      if (LEN_TRIM(line) > 63) then
        linend(i) = line(65:LEN_TRIM(line))
      else
        linend(i) = ''
      endif
    ENDDO
  endif

  !DO i=1, num_atom
  !      READ(24,'(A)') line
  !      READ(line, '(1X,A15,3F14.8)') atomes_gjf(i), coord_gau(1:3,i,)
  !      !if the line contains connections information between layers
  !      IF (LEN_TRIM(line) > 59) THEN
  !        linend(i) = line(60:LEN_TRIM(line))
  !      ELSE
  !        linend(i) = ''
  !      ENDIF
  !ENDDO

  !read to the end
  !READ(24,'(A)') c1
  !DO i=1, num_atom
  !    READ(24, '(A)') gjfend(i)
  !ENDDO
  num_end = 0
  READ(24,'(A)') c1
  do while(.true.)
    read(24, '(A)', IOSTAT = isend) gjfend(num_end+1)
    if (isend < 0) then
      EXIT
    else
      num_end = num_end +1
    endif
  enddo

  CLOSE(24)


  !write the new gjffile
  open(26, file = TRIM(gjffile), status='replace')
  do i = 1, num_head
    write(26,'(A)') trim(gjfheader(i))
  enddo

  if (line_len <60 ) then
    DO i=1, num_atom
        write(26, '(1X,A15,3F14.8)') atomes_gjf15(i), coords(3*(i-1)+1 : 3*i)
    ENDDO
  elseif (line_len < 63 ) then
    DO i=1, num_atom
      write(26, '(1X,A18,3F14.8)') atomes_gjf(i), coords(3*(i-1)+1 : 3*i)
    ENDDO
  else
    DO i=1, num_atom
      write(26, '(1X,A18,3F14.8,1X,A1)',advance='no') atomes_gjf(i), coords(3*(i-1)+1 : 3*i), trim(layer(i))
      write(26,*) trim(linend(i))  
    ENDDO
  endif
  write(26,'(A)') ''
  do i = 1, num_end
    write(26,'(A)') trim(gjfend(i))
  enddo
  !write(26,'(A)') ''
  !write(26,'(A)') ''
  close(26)

  INQUIRE(FILE = trim(chkfile), EXIST = b)

  !there is problem using external 
  !if ((.not. isguess) .and. b) then
  !  call system("sed -i 's/force/force guess=read/1' "//TRIM(gjffile))
  !elseif (isguess .and. ( .not. b)) then
  !  call system("sed -i 's/guess=read//1' "//TRIM(gjffile))
  !endif


 !coords = coords*0.52917721
  !open(26, file = TRIM(gjffile), status='replace')
  !do i = 1, num_head
  !  write(26,'(A)') trim(gjfheader(i))
  !enddo
!
  !do i = 1, num_atom
  !  if (len_trim(linend(i)) >0) then
  !    write(26,'(1X,A15,3F14.8)',advance='no') atomes_gjf(i), coords(3*(i-1)+1 : 3*i)*0.52917721
  !    write(26,*) trim(linend(i))
  !  else
  !    write(26,'(1X,A15,3F14.8)') atomes_gjf(i), coords(3*(i-1)+1 : 3*i)*0.52917721
  !  endif
  !enddo
  !write(26,'(A)') ''
  !do i = 1, num_atom
  !  write(26,'(A)') trim(gjfend(i))
  !enddo
  !write(26,'(A)') ''
  !write(26,'(A)') ''
  !close(26)

  !call system('$Gauexe < '//trim(gjffile)//' > gaussian.out')

endsubroutine

subroutine update_pdb(num_atom, pdbfile, coords)
  use dlf_parameter_module, only: rk
  use GauCNS_parameter_module
  implicit none
  integer, intent(in) :: num_atom
  character(len=*), intent(in) :: pdbfile
  real(rk), intent(in) ::  coords(3*max(num_atom, H_natom)) 

  CHARACTER(len=150), DIMENSION(100) :: pdbheader
  CHARACTER(len=26), DIMENSION(2,num_atom) :: pdblines
  real(rk), DIMENSION(3,num_atom) :: coord_pdb

  logical :: b
  integer :: num_head, i, j
  CHARACTER(len=150) :: c1
  real(rk), dimension(3) :: xyz, xyz1
  real(rk) :: x, y, z


  !test 
  !open(30, file = 'testlog.log', status = 'replace')

  !read old pdb file
  num_head = 0
  !write(30,*) num_head
  open(28, file = trim(pdbfile), action = 'read')
  b = .false.
  do while(.not. b)
    READ(28, '(A)') c1
    num_head = num_head + 1
    !write(30,*) num_head
    pdbheader(num_head) = trim(c1)
    b = ( c1(1:4) == 'ATOM')
  enddo
  num_head = num_head - 1
  BACKSPACE(28)

  !write(30,*) num_head
  !close(30)

  do i = 1, num_atom
    read(28,'(A26,4X,3F8.3,A26)') pdblines(1,i), coord_pdb(1:3,i), pdblines(2,i)
  enddo
  close(28)

  open(29, file = trim(pdbfile), status = 'replace')
  open(31, file = trim(pdbfile)//'1', status = 'replace')
  DO i=1,num_head
    WRITE(29, '(A)') TRIM(pdbheader(i))
    WRITE(31, '(A)') TRIM(pdbheader(i))
  ENDDO

  do i = 1, num_atom
    j = gjf_pdb_map(i)
    x = coords(3*(j-1)+1)
    y = coords(3*(j-1)+2)
    z = coords(3*j)

    !double
    !xyz(1) = IDNINT(x * 1000.0)/1000.0
    !xyz(2) = IDNINT(y * 1000.0)/1000.0
    !xyz(3) = IDNINT(z * 1000.0)/1000.0
!
    !xyz1(1) = (x - IDNINT(x * 1000.0)/1000.0)*100000.0
    !xyz1(2) = (y - IDNINT(y * 1000.0)/1000.0)*100000.0
    !xyz1(3) = (z - IDNINT(z * 1000.0)/1000.0)*100000.0

    !real
    xyz(1) = INT(x * 1000.0)/1000.0
    xyz(2) = INT(y * 1000.0)/1000.0
    xyz(3) = INT(z * 1000.0)/1000.0

    xyz1(1) = (x - INT(x * 1000.0)/1000.0)*100000.0
    xyz1(2) = (y - INT(y * 1000.0)/1000.0)*100000.0
    xyz1(3) = (z - INT(z * 1000.0)/1000.0)*100000.0

    write(29, '(A26,4X,3F8.3,A26)') pdblines(1,i), xyz(:), pdblines(2,i)
    write(31, '(A26,4X,3F8.3,A26)') pdblines(1,i), xyz1(:), pdblines(2,i)
  enddo
  write(29, *) 'END'
  write(29, *)
  write(31, *) 'END'
  write(31, *)
  close(29)
  close(31)

end subroutine

subroutine update_xyz(num_atom, xyzfile, coords)
  use GauCNS_parameter_module
  use dlf_parameter_module, only: rk
  implicit none
  integer, intent(in) :: num_atom
  character(len=*), intent(in) :: xyzfile
  real(rk), intent(in) ::  coords(3*num_atom)

  character(len = 4), dimension(num_atom) :: ele
  character(len=100) :: cline
  integer :: num ,i
  integer, dimension(3, num_atom) :: coord_temps
  !read old xyz file
  open(32, file = trim(xyzfile), action = 'read')
  read(32, *) num
  if (num /= num_atom) then
    write(*,*) 'Number of atoms in coord.xyz is different with parameter definition'
    stop
  endif

  read(32, '(A)') cline
  do i = 1, num_atom 
    read(32, *) ele(i), coord_temps(:,i)
  enddo
  close(32)

  !write new xyz file
  open(33, file = trim(xyzfile), status = 'replace')
  write(33, *) num_atom

  write(33, '(A)') trim(cline)
  do i = 1, num_atom 
    write(33, '(A4,3F14.6)') ele(i), coords(3*(i-1)+1:3*i)
  enddo
  write(33,*) ''
  close(33)

endsubroutine

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_error()
  implicit none
! **********************************************************************
  call dlf_mpi_abort() ! only necessary for a parallel build;
                       ! can be present for a serial build
  call exit(1)
end subroutine dlf_error

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_update()
  implicit none
! **********************************************************************
  ! only a dummy routine here.
end subroutine dlf_update


subroutine dlf_get_multistate_gradients(nvar,coords,energy,gradient,iimage,status)
  ! only a dummy routine up to now
  ! for conical intersection search
  use dlf_parameter_module
  implicit none
  integer   ,intent(in)    :: nvar
  integer   ,intent(in)    :: coords(nvar)
  real(rk)  ,intent(in)    :: energy(2)
  real(rk)  ,intent(in)    :: gradient(nvar,2)
  integer   ,intent(in)    :: iimage
  integer   ,intent(in)    :: status
end subroutine dlf_get_multistate_gradients


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_put_procinfo(dlf_nprocs, dlf_iam, dlf_global_comm)

  implicit none

  integer, intent(in) :: dlf_nprocs ! total number of processors
  integer, intent(in) :: dlf_iam ! my rank, from 0, in mpi_comm_world
  integer, intent(in) :: dlf_global_comm ! world-wide communicator
! **********************************************************************

!!! variable in the calling program = corresponding dummy argument

end subroutine dlf_put_procinfo


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_procinfo(dlf_nprocs, dlf_iam, dlf_global_comm)

  implicit none

  integer :: dlf_nprocs ! total number of processors
  integer :: dlf_iam ! my rank, from 0, in mpi_comm_world
  integer :: dlf_global_comm ! world-wide communicator
! **********************************************************************

!!! dummy argument = corresponding variable in the calling program

end subroutine dlf_get_procinfo


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_put_taskfarm(dlf_ntasks, dlf_nprocs_per_task, dlf_iam_in_task, &
                        dlf_mytask, dlf_task_comm, dlf_ax_tasks_comm)

  implicit none

  integer, intent(in) :: dlf_ntasks          ! number of taskfarms
  integer, intent(in) :: dlf_nprocs_per_task ! no of procs per farm
  integer, intent(in) :: dlf_iam_in_task     ! my rank, from 0, in my farm
  integer, intent(in) :: dlf_mytask          ! rank of my farm, from 0
  integer, intent(in) :: dlf_task_comm       ! communicator within each farm
  integer, intent(in) :: dlf_ax_tasks_comm   ! communicator involving the 
                                             ! i-th proc from each farm
! **********************************************************************

!!! variable in the calling program = corresponding dummy argument 

end subroutine dlf_put_taskfarm


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_taskfarm(dlf_ntasks, dlf_nprocs_per_task, dlf_iam_in_task, &
                        dlf_mytask, dlf_task_comm, dlf_ax_tasks_comm)

  implicit none

  integer :: dlf_ntasks          ! number of taskfarms
  integer :: dlf_nprocs_per_task ! no of procs per farm
  integer :: dlf_iam_in_task     ! my rank, from 0, in my farm
  integer :: dlf_mytask          ! rank of my farm, from 0
  integer :: dlf_task_comm       ! communicator within each farm
  integer :: dlf_ax_tasks_comm   ! communicator involving the
                                 ! i-th proc from each farm
! **********************************************************************

!!! dummy argument = corresponding variable in the calling program

end subroutine dlf_get_taskfarm


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_output(dum_stdout, dum_stderr)
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr,stdout,keep_alloutput
  implicit none
  integer :: dum_stdout
  integer :: dum_stderr
  integer :: ierr
  logical :: topened
  character(len=10) :: suffix

! sort out output units; particularly important on multiple processors
 
! set unit numbers for main output and error messages
  if (dum_stdout >= 0) stdout = dum_stdout 
  if (dum_stderr >= 0) stderr = dum_stderr

  if (glob%iam /= 0) then
     inquire(unit=stdout, opened=topened, iostat=ierr)
     if (topened .and. ierr == 0) close(stdout)
     if (keep_alloutput) then ! hardwired in dlf_global_module.f90
        write(suffix,'(i10)') glob%iam
        open(unit=stdout,file='output.proc'//trim(adjustl(suffix)))
     else
        open(unit=stdout,file='/dev/null')
     end if
  endif

  if (glob%nprocs > 1) then
     ! write some info on the parallelization
     write(stdout,'(1x,a,i10,a)')"I have rank ",glob%iam," in mpi_comm_world"
     write(stdout,'(1x,a,i10)')"Total number of processors = ",glob%nprocs
     if (keep_alloutput) then
        write(stdout,'(1x,a)')"Keeping output from all processors"
     else
        write(stdout,'(1x,a)')"Not keeping output from processors /= 0"
     end if
  end if

end subroutine dlf_output


! **********************************************************************
! **********************************************************************
! The following routine either writes random numbers to a file, or reads
! them. This is to have equal starting conditions for different compilers
subroutine read_rand(arr)
  use dlf_parameter_module, only: rk
  use dlf_global, only : glob
  real(rk)  :: arr(:)
  integer, parameter :: si1=12
  integer, parameter :: si2=3000
  logical,parameter :: readf=.true.
  real(rk) :: ar(si2)
  integer :: l(1),length
  l=ubound(arr)
  length=l(1)
  if(readf) then
    if(length<=si1) then
      open(unit=201,file="random1.bin",form="unformatted")
    else if(length<=si2) then
      open(unit=201,file="random2.bin",form="unformatted")
    else
      call dlf_mpi_finalize() ! only necessary for a parallel build;
                              ! can be present for a serial build
      stop "Too many coordinates to be read from random.bin file"
    end if
    read(201) ar(1:length)
    close(201)
    arr=ar(1:length)
  else
    if (glob%iam == 0) then
       call random_number(ar)
       open(unit=201,file="random1.bin",form="unformatted")
       write(201) ar(1:si1)
       close(201)
       open(unit=201,file="random2.bin",form="unformatted")
       write(201) ar(1:si2)
       close(201)
    end if
  end if
end subroutine read_rand

subroutine driver_eck(vb_,alpha_)
  use dlf_parameter_module
  !use driver_parameter_module
  real(rk),intent(out) :: vb_,alpha_
  vb_=vb
  alpha_=alpha
end subroutine driver_eck
