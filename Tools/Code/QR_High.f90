program QR_Hign
    implicit none

    integer :: narg, numH, i, isok, frez
    character(len=100) :: gjff, xyzf
    character(len=1) :: layer

    doubleprecision, dimension(:,:), allocatable :: coordH
    doubleprecision :: x, y, z 
    character(len=150) :: line, c1, c2
    character(len=20) :: ele
    logical :: b, ornot


    narg=iargc()
    if (narg == 3) then
        call getarg(1, gjff)
        call getarg(2, xyzf)
        call getarg(3, layer)

        !read xyz file
        open(10, file=trim(xyzf), action='read')
        read(10, *) numH
        read(10, *)
        allocate(coordH(3, numH))
        do i=1,numH
            read(10, *) ele, coordH(1, i), coordH(2, i), coordH(3, i)
        enddo
        close(10)

        !read gjf file
        open(11, file=trim(gjff),action='read')
        open(12, file='new_'//trim(gjff),status='replace')
        b=.false.
        do while(.not. b)
            read(11, '(A)') line
            write(12, '(A)') trim(line)
            b=(trim(line) == '')
        enddo

        b=.false.
        do while(.not. b)
            read(11, '(A)') line
            write(12, '(A)') trim(line)
            b=(trim(line) == '')
        enddo

        !charge, spin
        read(11, '(A)') line
        write(12, '(A)') trim(line)

        !coordinate
        b=.false.
        do while(.not. b)
            read(11, '(A)') line
            b=(trim(line) == '')
            if (.not. b) then
                read(line, *) c1, frez, x, y, z, c2
                call test_pos((/x, y, z/), coordH, numH, ornot)
            endif
           
            if ( ornot .and. (.not. b)) then
                write(12, '(1X,A16,I2,3F14.8,1X,A1)') c1(1:16), frez, x, y, z, layer
            else
                write(12, '(A)') trim(line)
            endif                        
        enddo
        !connectivity and MM parameters
        do while (.true.)
            read(11, '(A)', iostat=isok) line
            if (isok ==0) then
                write(12, '(A)') trim(line)
            else
                exit      
            endif

        enddo

        close(11)
        close(12)

    else
        write(*,*) 'This program is to modifiy the layer atom of gjf file'
        write(*,*) 'Input: mol.gjf, High.xyz H (H or M)'
    endif 

endprogram

subroutine test_pos(pos, mat, dim, isin)
    implicit none
    logical, intent(out) :: isin
    integer, intent(in) :: dim 
    doubleprecision, dimension(3) :: pos 
    doubleprecision, dimension(3,dim) :: mat 
    doubleprecision :: diff

    integer :: i 

    isin = .false.
    do i=1, dim 
        diff = ((pos(1)-mat(1,i))**2 + (pos(2)-mat(2,i))**2 +(pos(3)-mat(3,i))**2)**0.5
        if (diff < 0.01) then
            isin = .true.
            exit
        endif 
    enddo
endsubroutine