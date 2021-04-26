program main
    integer :: narg, num_atom
    logical :: alive

    character(len=120) :: gjffile_ori

    narg = iargc()

    if (narg == 1) then
        call getarg(1, gjffile_ori)

        INQUIRE(file = trim(gjffile_ori), exist= alive )
        if (.not. alive) then
            write(*,*) 'File '//trim(gjffile_ori)//' does not exist'
            stop
        endif

        call get_num(gjffile_ori, num_atom)

        call check_gjflink(gjffile_ori, num_atom)

    endif

endprogram

subroutine check_gjflink(gjfin, num_atom)
    implicit none
    character(len=*), intent(in) :: gjfin
    integer, intent(in) :: num_atom

    integer :: i, num_head, num_end, isend, line_len, pos
    logical :: b 
    DOUBLEPRECISION, DIMENSION(3,num_atom):: coord_gau
    CHARACTER(len=16), DIMENSION(num_atom):: atomes_gjf
    !CHARACTER(len=15), DIMENSION(num_atom):: atomes_gjf15
    CHARACTER(len=150), DIMENSION(num_atom*2):: gjfend
    CHARACTER(len=1), DIMENSION(num_atom) :: layer
    CHARACTER(len=150), DIMENSION(10) :: gjfheader
    CHARACTER(len=15), DIMENSION(num_atom) :: linend
    integer, dimension(num_atom) :: freeze, link_info
    CHARACTER(len=15) :: c1

    CHARACTER(len=120) :: line

    !read gjf information from gjfin

    num_head = 0
    OPEN(24, FILE = TRIM(gjfin), ACTION = 'read')
    b= .FALSE.
    DO WHILE(.NOT. b)
            num_head = num_head +1
            READ(24, '(A)') gjfheader(num_head)
            b = (TRIM(gjfheader(num_head)) == '')
    ENDDO

    b= .FALSE.
    DO WHILE(.NOT. b)
            num_head = num_head +1
            READ(24, '(A)') gjfheader(num_head)
            b = (TRIM(gjfheader(num_head)) == '')
    ENDDO
    num_head = num_head +1
    READ(24, '(A)') gjfheader(num_head)

!read atome block (with no PDB information)
    read(24, '(A)') line
    line_len = len_trim(line)
    backspace(24)
    if (line_len < 63 ) then
        write(*,*) 'there is no freeze and Oniom information, no need to check'
        stop
    else

        DO i=1, num_atom

          READ(24,'(A)') line

          READ(line, '(1X,A16,I2,3F14.8,1X,A1)') atomes_gjf(i), freeze(i), coord_gau(1:3,i), layer(i)

          if (LEN_TRIM(line) > 63) then
            linend(i) = line(65:LEN_TRIM(line))
            READ(linend(i), *) c1, link_info(i)
          else
            linend(i) = ''
            link_info(i) = -1
          endif
        ENDDO
    endif

    num_end = 0
    READ(24,'(A)') line
    do while(.true.)
      read(24, '(A)', IOSTAT = isend) gjfend(num_end+1)
      if (isend < 0) then
        EXIT
      else
        num_end = num_end +1
      endif
    enddo
    CLOSE(24)

!check the link infomation
    do i = 1, num_atom
        if (link_info(i) > 0) then
            if (freeze(i) /= freeze(link_info(i))) then
                write(*,*) 'Link atoms', i, ' and', link_info(i), ' have different freeze infomation!'
                freeze(i) = 0
                freeze(link_info(i)) = 0
            endif
        endif

    enddo

!write to the new gjf file
    pos = index(TRIM(gjfin), '.')-1
    open(26, file = gjfin(1:pos)//'_new'//'.gjf', status='replace')
    do i = 1, num_head
      write(26,'(A)') trim(gjfheader(i))
    enddo
  
    DO i=1, num_atom
        write(26, '(1X,A16,I2,3F14.8,1X,A1)',advance='no') atomes_gjf(i), freeze(i), coord_gau(1:3,i), trim(layer(i))
        write(26,*) trim(linend(i))  
    ENDDO

    write(26,'(A)') ''
    do i = 1, num_end
      write(26,'(A)') trim(gjfend(i))
    enddo
    !write(26,'(A)') ''
    !write(26,'(A)') ''
    close(26)
    write(*,*) 'New gjf with link atoms relaxed has been generated!'

endsubroutine


subroutine get_num(gjfin, num_atom)
    implicit none
    character(len=*), intent(in) :: gjfin
    integer, intent(out) :: num_atom

    logical :: b 
    character(len=100) :: c


    open(21, file = trim(gjfin), action = 'read')

    b = .false.
    DO WHILE(.NOT. b)
        read (21, '(A)') c
        b = (TRIM(c) == '')
    ENDDO

    b = .false.
    DO WHILE(.NOT. b)
        read (21, '(A)') c
        b = (TRIM(c) == '')
    ENDDO

    read (21, '(A)') c

! atom block
    b = .false.
    num_atom = 0
    DO WHILE(.NOT. b)
        num_atom = num_atom + 1
        read (21, '(A)') c
        b = (TRIM(c) == '')
    ENDDO
    num_atom = num_atom - 1

    close(21)

endsubroutine

