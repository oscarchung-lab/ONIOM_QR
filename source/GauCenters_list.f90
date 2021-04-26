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
    use para_ele
    implicit none
    integer :: narg, nb_center, i, j, num, headline, line_len, num_line, isend
    character(len=120) :: parafile, InputFile, OutputFile, MsgFile, FChkFile, MatElFile
    character(len=120) :: gjffile, headf, tailf
    character(len=120) :: line, c1, c2, c3
    character(len=30) :: num_str

    
    integer :: atoms, deriva, charge, spin

    integer, dimension(:), allocatable :: atomic, nb 
    doubleprecision, dimension(:,:), allocatable :: atomiccoord
    integer, dimension(:,:), allocatable :: list_gjf

    character(len=120), dimension(:,:), allocatable :: headtailfile
    character(len=150), dimension(:), allocatable :: liststr
    
    doubleprecision, dimension(:,:,:), allocatable :: grad
    doubleprecision, dimension(:), allocatable :: ene

    logical :: alive, b

    narg=iargc()

    1001 FORMAT(4I10)
    1002 FORMAT(I10, 3F20.12)

    !test case for checking the atomic label in external atom list
    if (narg == 6) then
        call getarg(2, InputFile) 
        !read gaussian input
        open(11,file = trim(InputFile), action='read')
        open(13,file = 'testtmp.gjf', status='replace')
        write(13,'(A)') '%chk=testtmp.chk'
        write(13,'(A)') '# hf/sto-3g '
        write(13,'(A)') ' '
        write(13,'(A)') 'gjf for external atom list'
        write(13,'(A)') ' '

        allocate( atomic(atoms), atomiccoord(3,atoms))
        read(11,1001) atoms, deriva, charge, spin

        write(13,"(2I4)") charge, spin

        do i = 1, atoms
            read(11, 1002) atomic(i), atomiccoord(:,i)
            atomiccoord(:,i) = atomiccoord(:,i) * b2a
            write(13, "(A2,1X,3F20.12)") element(atomic(i)), atomiccoord(:,i)
        enddo
        write(13, *) 
        write(13, *) 
        close(13)
        close(11)
        deallocate(atomic,atomiccoord)
    elseif (narg == 7) then
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
    
        read(10,*) nb_center
        
        
        allocate(liststr(nb_center),headtailfile(2,nb_center),nb(nb_center), ene(nb_center))
        
        do i=1,nb_center
            read(10,"(A)") liststr(i)
            
        enddo

        !head and tail file for Gaussian computation
        do i=1,nb_center
            read(10,*) headtailfile(:,i)
            
        enddo
        close(10)

        
        if(nb_center > 1) then
            
            !read gaussian input
            open(11,file = trim(InputFile), action='read')
         
            read(11,1001) atoms, deriva, charge, spin
            allocate(list_gjf(atoms,2), atomic(atoms), atomiccoord(3,atoms))
         
            do i = 1, atoms
                read(11, 1002) atomic(i), atomiccoord(:,i) 
                atomiccoord(:,i) = atomiccoord(:,i) * b2a
            enddo
            close(11)

           

            do j = 1,nb_center
                call getatomfrag(liststr(j), atoms, list_gjf, nb(j), j)
            enddo
            
            !save the gjfs of framentations in a backup file for checking
            call system(':> center_gjf.bak')

            allocate(grad(3,maxval(nb),nb_center))
            do j=1,nb_center 
                inquire(file=trim(headtailfile(1,j)), exist=b)
                if (.not. b) then
                    write(*,*) trim(headtailfile(1,j))//' does not exist'
                    stop 
                endif 
               
                !read the head  file
                open(23,file=trim(headtailfile(1,j)),action='read')
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
                inquire(file=trim(headtailfile(2,j)), exist=b)
                !withno extral information
                if (.not. b) then
                    write(13,*) ' '
                    write(13,*) ' '
                    write(13,*) ' '
                else
                    open(33,file=trim(headtailfile(2,j)),action='read')
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
            do i=1,nb_center
                write(*,*) 'Fragment ', i, ': ', ene(i), ' Bohrs'
            enddo
        
            open(25, file=trim(OutputFile),status='replace')
            write(25,"(4D20.12)") sum(ene), 0D0,0D0,0D0
            do i=1, atoms
                write(25,"(3D20.12)") grad(:,list_gjf(i,2),list_gjf(i,1))
            enddo
            close(25)
            call system('rm -f fort.7 reslt.tmp centers.dat')
            deallocate(liststr,headtailfile,nb,ene)
            deallocate(list_gjf,atomic,atomiccoord)

        else 
            write(*,*) 'Only one center is not the case.'
        endif
    else
        write(*,*) "The input is: parfile layer InputFile OutputFile MsgFile FChkFile MatElFile"
        write(*,*) "See Gaussian External keyword description for more information"
    endif

endprogram


!str2list
subroutine getatomfrag(strinfo, numatoms, lisgjf, numfrag, numcenter)
    implicit none
    character(len=*), intent(in) :: strinfo
    integer, intent(in) :: numatoms, numcenter
    integer, dimension(numatoms,2) :: lisgjf
    integer, intent(out) :: numfrag

    character(len=len(strinfo)) :: tmpstrs, tmpstr, str
    character(len=20) :: strp
    integer :: num_temp, num1, num2, i, numtmp

    !open(10, file=trim(listname),status='replace')

    tmpstrs=trim(adjustl(strinfo))
    !write(*,*) tmpstrs
    num_temp=0

    do while (len_trim(tmpstrs) > 0 .and. index(tmpstrs, ",") >0 )
        !write(*,*) index(tmpstrs, ",")
        !write(*,*) tmpstrs
        tmpstr=tmpstrs(1:index(tmpstrs, ",")-1)
        !write(*,*) tmpstr

        if (index(tmpstr, "-") == 0) then !only one atom for this part e.g.: ,num,
            num_temp=num_temp+1
            !write(10, *) trim(tmpstr)
            read(tmpstr, *) numtmp
            lisgjf(numtmp,1) = numcenter
            lisgjf(numtmp,2) = num_temp
        else !a list of atoms for this part e.g.: ,num1-num2, 
            read(tmpstr(1:index(tmpstr, "-")-1), *) num1
            read(tmpstr(index(tmpstr, "-")+1:), *) num2
            do i=num1,num2
                num_temp=num_temp+1
                lisgjf(i,1) = numcenter
                lisgjf(i,2) = num_temp
                !write(strp, *) i
                !write(10, *) adjustl(trim(strp))
            enddo
        endif
        tmpstrs = trim(adjustl(tmpstrs(index(tmpstrs, ",")+1:)))
    enddo
    tmpstr=tmpstrs
    if (index(tmpstr, "-") == 0) then !only one atom for this part e.g.: ,num,
        num_temp=num_temp+1
        !write(10, *) trim(tmpstr)
        read(tmpstr, *) numtmp
        lisgjf(numtmp,1) = numcenter
        lisgjf(numtmp,2) = num_temp
    else !a list of atoms for this part e.g.: ,num1-num2, 
        read(tmpstr(1:index(tmpstr, "-")-1), *) num1
        read(tmpstr(index(tmpstr, "-")+1:), *) num2
        do i=num1,num2
            num_temp=num_temp+1
            lisgjf(i,1) = numcenter
            lisgjf(i,2) = num_temp
            !write(strp, *) i
            !write(10, *) adjustl(trim(strp))
        enddo
    endif
    numfrag=num_temp
    !close(10)
endsubroutine






