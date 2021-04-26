Module para_ele
    DOUBLEPRECISION, PARAMETER :: b2a=0.52917721067121D0
    CHARACTER(len=2), PARAMETER :: element(0:118)=(/'Bq','H ','He',& !0 is ghost atom, 1~2
        'Li','Be','B ','C ','N ','O ','F ','Ne',& !3~10
        'Na','Mg','Al','Si','P ','S ','Cl','Ar',& !11~18
        'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',& !19~36
        'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe',& !37~54
        'Cc','Ba',& !55~56
        'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',& !La, 57~71
        'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',& !72~86
        'Fr','Ra',& !87~88
        'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',& !Ac, 89~103
        'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/) !104~118
endmodule para_ele


program Gau2Centers
    implicit none

    integer :: narg, nb_center, i, j, num, headline, line_len, num_line
    character(len=120) :: parafile, InputFile, OutputFile, MsgFile, FChkFile, MatElFile
    character(len=120) :: gjffile, headf, tailf
    character(len=120) :: line, c1, c2, c3
    character(len=30) :: num_str

    

    integer, dimension(:), allocatable :: list_center, input_atomic, nb_gjf
    integer, dimension(:,:), allocatable :: pos_gjf

    !doubleprecision, dimension(:), allocatable :: r_center
    doubleprecision, dimension(:,:), allocatable :: coord_center
    character(len=120), dimension(:,:), allocatable :: headtailfile
    
    doubleprecision, dimension(:,:), allocatable :: coord_input
    
    logical :: alive, b

    narg=iargc()

    if (narg == 7) then
        call getarg(1,parafile)
        call getarg(3,InputFile)
        call getarg(4,OutputFile)
        call getarg(5,MsgFile)

        inquire(file=trim(parafile), exist=alive)
        if (.not. alive) then
            write(*,*) trim(parafile)//' does not exist'
            stop 
        endif 

        open(10,file=trim(parafile), action='read')
        read(10,*) gjffile
        read(10,*) nb_center
        !allocate(list_center(nb_center), r_center(nb_center), coord_center(3,nb_center), headtailfile(2,nb_center))
        allocate(list_center(nb_center), coord_center(3,nb_center), headtailfile(2,nb_center))
        do i=1,nb_center
            read(10,*) list_center(i)! r_center(i)
        enddo

        !head and tail file for Gaussian computation
        do i=1,nb_center
            read(10,*) headtailfile(:,i)
        enddo
        close(10)

        !gjf2xyz method
        !use sed n line info to get the center atom lines
        open(11,file=trim(gjffile),action='read')

        headline=0d0
        b = .false.
        do while(.not. b)
            read(11,'(A)') line
            headline=headline+1d0
            b = (trim(line) == '')
        enddo

        b = .false.
        do while(.not. b)
            read(11,'(A)') line
            headline=headline+1d0
            b = (trim(line) == '')
        enddo

        headline=headline+1d0

        read(11,'(A)') line
        read(11,'(A)') line
        !line_len = len_trim(line)
        call getncol(trim(line), line_len)
        close(11)

        !sed the line of center atoms
        !write(*,*) list_center(1)
        num_line=list_center(1)+headline
        write(num_str,*) num_line
        call system("sed -n "//trim(adjustl(num_str))//"p "//trim(gjffile)//" > centers.dat")
        do i=2,nb_center
            !write(*,*) list_center(i)
            num_line=list_center(i)+headline
            write(num_str,*) num_line
            call system("sed -n "//trim(adjustl(num_str))//"p "//trim(gjffile)//" >> centers.dat")
        enddo

        open(12,file="centers.dat",action='read')

        if (line_len == 4) then
            do i= 1, nb_center
                read(12,*) c1, coord_center(1:3,i)
            enddo
        else
            do i= 1, nb_center
                read(12,*) c1, c2, coord_center(1:3,i)
            enddo
        endif

        close(12)
        
        !if(nb_center==2) then
        if(nb_center > 1) then
            !call twocenter(coord_center, headtailfile, trim(InputFile), trim(OutputFile))
            call multicenter(nb_center, coord_center, headtailfile, trim(InputFile), trim(OutputFile))
        else 
            write(*,*) 'Only one center is not the case.'
        endif
        
        !1001 FORMAT(4I10)
        !1002 FORMAT(I10, 3F20.12)
        !open(12,file=trim(InputFile),action='read')
        !read(12,1001) atoms, deriva, charge, spin
        !allocate(coord_input(3,atoms), input_atomic(atoms), pos_gjf(atoms,nb_center))
        !do i=1,atoms
        !    read(11, 1002) input_atomic(i), coord_input(1,i), coord_input(2,i), coord_input(3,i)
        !enddo
        !close(12)

    else
        write(*,*) "The input is: parfile layer InputFile OutputFile MsgFile FChkFile MatElFile"
        write(*,*) "See Gaussian External keyword description for more information"
    endif

endprogram


subroutine getncol(strs, ncol)
    implicit none
    character(len=*), intent(in) :: strs 
    integer, intent(out) :: ncol
    character(len=len(strs)) :: tmpstrs, tmpstr
    integer :: pos

    tmpstrs = trim(adjustl(strs))
    if(len_trim(tmpstrs) > 0) then
        pos=index(trim(tmpstrs), " ")
        if (pos==0) then
            ncol=1
        else
            ncol=0
            do while (len_trim(tmpstrs) > 0)
                tmpstr=tmpstrs(1:index(tmpstrs, " ")-1)
                ncol = ncol + 1
                tmpstrs = trim(adjustl(tmpstrs(index(tmpstrs, " "):)))
            enddo
        endif
    else
        ncol=0
    endif
endsubroutine


!subroutine multcenter(coord_center, r_center, headtails, inputfile, outputtfile)
subroutine multicenter(numcenter, coord_center, headtails, inputfile, outputtfile)
    use para_ele
    implicit none
    integer, intent(in) :: numcenter
    doubleprecision, dimension(3,numcenter), intent(in) :: coord_center
    !doubleprecision, dimension(2), intent(in) :: r_center
    character(len=120), dimension(2,numcenter), intent(in) :: headtails
    character(len=*), intent(in) :: inputfile, outputtfile

    integer, dimension(numcenter) :: nb_gjf 
    integer, dimension(:,:), allocatable :: list_gjf
    doubleprecision, dimension(:,:,:), allocatable :: grad
    character(len=120) :: line, str_num

    integer, dimension(:), allocatable :: atomic
    doubleprecision, dimension(:,:), allocatable :: atomiccoord

    doubleprecision :: mindis
    doubleprecision, dimension(numcenter) :: dis, ene 
    integer :: minpos
    integer, dimension(numcenter) :: nb

    integer :: i, isend, j
    integer :: atoms, deriva, charge, spin

    logical :: b
    
    nb=0d0

    1001 FORMAT(4I10)
    1002 FORMAT(I10, 3F20.12)

   
    !separate the atoms from external input file to fragments

    !open InputFile files
    open(12,file=trim(InputFile),action='read')
    read(12,1001) atoms, deriva, charge, spin
    allocate(list_gjf(atoms,2), atomic(atoms), atomiccoord(3,atoms))
    
    do i=1,atoms
       
        read(12, 1002) atomic(i), atomiccoord(:,i)
        !gjf unit Angstrom
        !input unit a.u
        !we use Angstrom in the parameter file
        atomiccoord(:,i) = atomiccoord(:,i) * b2a
  
        do j=1,numcenter
            dis(j)=((atomiccoord(1,i)-coord_center(1,j))**2 + (atomiccoord(2,i)-coord_center(2,j))**2 + (atomiccoord(3,i)-coord_center(3,j))**2)**0.5
        enddo
    
        mindis=minval(dis)
        minpos=findloc(dis,mindis,1)

        nb(minpos)=nb(minpos)+1
        !list_gjf(i,1) the fragment of atom i
        list_gjf(i,1)=minpos
        !list_gjf(i,2) the index of atom i in the fragment
        list_gjf(i,2)=nb(minpos)
        
    enddo
    close(12)

    !save the gjfs of framentations in a backup file for checking
    call system(':> center_gjf.bak')

    allocate(grad(3,maxval(nb),numcenter))
    do j=1,numcenter 
        inquire(file=trim(headtails(1,j)), exist=b)
        if (.not. b) then
            write(*,*) trim(headtails(1,j))//' does not exist'
            stop 
        endif 
       
        !read the head  file
        open(23,file=trim(headtails(1,j)),action='read')
        !write gjf of fragment 
        open(13,file='Center_temp.gjf',status='replace')

        b = .false.
        do while(.not. b)
            read(23,'(A)') line
            write(13,*) trim(line)
            b = (trim(line) == '')
        enddo

        b = .false.
        do while(.not. b)
            read(23,'(A)') line
            write(13,*) trim(line)
            b = (trim(line) == '')
        enddo

        read(23,'(A)') line
        write(13,*) trim(line)
        close(23)

        !write atom block
        do i=1,atoms 
            if ((list_gjf(i,1)) == j ) then
                write(13,'(1X,A2,13X,3F14.8)') element(atomic(i)),atomiccoord(:,i)
            endif
        enddo

        !read the tail
        inquire(file=trim(headtails(2,j)), exist=b)
        !withno extral information
        if (.not. b) then
            write(13,*) ' '
            write(13,*) ' '
            write(13,*) ' '
        else
            open(33,file=trim(headtails(2,j)),action='read')
            do while(.true.)
                read(33, '(A)', IOSTAT = isend) line
                if (isend < 0) then
                    EXIT
                else
                    write(13,*) trim(line)
                endif
            enddo
            close(33)
        endif
        write(13,*) ' '
        close(13)

        call system('cat Center_temp.gjf >> center_gjf.bak')
        call system('echo >> center_gjf.bak')
        !run gaussian job
        call system('$Gauexe < Center_temp.gjf > Center_temp.out')
        !extract energy and force 
        call system("grep 'extrapolated energy' Center_temp.out | awk '{print $5}' > result.temp ")
        open(16, file='result.temp',action='read')
        if(EOF(16)) then
        close(16)
        call system("grep 'SCF Done:' Center_temp.out | awk '{print $5}' > result.temp ")
        else
            close(16)
        endif

        !get energy
        open(20, file='result.temp', action='read')
        read(20,*) ene(j)
        close(20)

        open(21,file='fort.7', action='read')

        do i= 1, nb(j)
            read(21,*)
        enddo

        do i= 1, nb(j)
            read(21,*) grad(:,i,j)
        enddo
        close(21)

    enddo
    
    !write output file
    write(*,*) 'Energy of fragments:'
    do i=1,numcenter
        write(*,*) 'Fragment ', i, ': ', ene(i), ' Bohrs'
    enddo

    open(25, file=trim(outputtfile),status='replace')
    write(25,"(4D20.12)") sum(ene), 0D0,0D0,0D0
    do i=1, atoms
        write(25,"(3D20.12)") grad(:,list_gjf(i,2),list_gjf(i,1))
    enddo
    close(25)
    call system('rm -f fort.7 reslt.tmp centers.dat')
    deallocate(list_gjf,atomic,atomiccoord,grad)
  
endsubroutine





