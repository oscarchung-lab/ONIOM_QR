program main
    integer :: narg, num_atom
    logical :: alive

    character(len=120) :: gjffile_ori, gjffile_tar

    narg = iargc()

    if (narg == 2) then
        call getarg(1, gjffile_ori)
        call getarg(2, gjffile_tar)

        INQUIRE(file = trim(gjffile_ori), exist= alive )
        if (.not. alive) then
            write(*,*) 'File '//trim(gjffile_ori)//' does not exist'
            stop
        endif

        INQUIRE(file = trim(gjffile_tar), exist= alive )
        if (.not. alive) then
            write(*,*) 'File '//trim(gjffile_tar)//' does not exist'
            stop
        endif

        call get_num(gjffile_ori, num_atom)

        call update_gjf(trim(gjffile_ori), trim(gjffile_tar), num_atom)

    else
        write(*,*) 'Wrong input'
        write(*,*) 'give me the coodinate gjf and target gjf file'
    endif

end program

subroutine update_gjf(gjfin, gjfout, num_atom)
    implicit none
    character(len=*), intent(in) :: gjfin, gjfout
    integer, intent(in) :: num_atom
    integer :: i, num_head, num_end, isend, line_len
    logical :: b 
    DOUBLEPRECISION, DIMENSION(3,num_atom):: coord_gau, coord_temp
    CHARACTER(len=18), DIMENSION(num_atom):: atomes_gjf
    CHARACTER(len=15), DIMENSION(num_atom):: atomes_gjf15
    CHARACTER(len=150), DIMENSION(num_atom*2):: gjfend
    CHARACTER(len=1), DIMENSION(num_atom) :: layer
    CHARACTER(len=150), DIMENSION(10) :: gjfheader
    CHARACTER(len=15), DIMENSION(num_atom) :: linend

    CHARACTER(len=120) :: line


    !read coord from gjfin
    open(20, file = trim(gjfin), action = 'read')
    b = .false.
    DO WHILE(.NOT. b)
        read (20, '(A)') line
        b = (TRIM(line) == '')
    ENDDO

    b = .false.
    DO WHILE(.NOT. b)
        read (20, '(A)') line
        b = (TRIM(line) == '')
    ENDDO

    read (20, '(A)') line
    

! atom block
    b = .false.
    read (20, '(A)') line
    line_len = len_trim(line)
    backspace(20)
    if (line_len < 60) then
        do i = 1, num_atom
            read(20, '(1X,A15,3F14.8)') atomes_gjf15(i), coord_gau(1:3,i)
            linend(i) = ''
        enddo
    elseif (line_len < 63) then
        do i = 1, num_atom
            read(20, '(1X,A18,3F14.8)') atomes_gjf(i), coord_gau(1:3,i)
            linend(i) = ''
        enddo
    else
        do i = 1, num_atom
            READ(20,'(1X,A18,3F14.8,1X,A1)') atomes_gjf(i), coord_gau(1:3,i), layer(i)
        enddo
    endif

    close(20)

    !read gjf information from gjfout
    num_head = 0
    OPEN(24, FILE = TRIM(gjfout), ACTION = 'read')
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
    if (line_len < 60 ) then
        DO i=1, num_atom
            READ(24,'(A)') line
            READ(line, '(1X,A15,3F14.8)') atomes_gjf15(i), coord_temp(1:3,i)
            linend(i) = ''
        ENDDO
      elseif (line_len < 63 ) then
        DO i=1, num_atom
            READ(24,'(A)') line
            READ(line, '(1X,A18,3F14.8)') atomes_gjf(i), coord_temp(1:3,i)
            linend(i) = ''
        ENDDO
      else
        DO i=1, num_atom
          READ(24,'(A)') line
          READ(line, '(1X,A18,3F14.8,1X,A1)') atomes_gjf(i), coord_temp(1:3,i), layer(i)
          if (LEN_TRIM(line) > 63) then
            linend(i) = line(65:LEN_TRIM(line))
          else
            linend(i) = ''
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
!write the new gjffile
    open(26, file = TRIM(gjfout), status='replace')
    do i = 1, num_head
      write(26,'(A)') trim(gjfheader(i))
    enddo
  
    if (line_len <60 ) then
      DO i=1, num_atom
          write(26, '(1X,A15,3F14.8)') atomes_gjf15(i), coord_gau(1:3,i)
      ENDDO
    elseif (line_len < 63 ) then
      DO i=1, num_atom
        write(26, '(1X,A18,3F14.8)') atomes_gjf(i), coord_gau(1:3,i)
      ENDDO
    else
      DO i=1, num_atom
        write(26, '(1X,A18,3F14.8,1X,A1)',advance='no') atomes_gjf(i), coord_gau(1:3,i), trim(layer(i))
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