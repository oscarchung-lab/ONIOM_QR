program main
    integer :: narg, num_atom, i
    logical :: alive
    character(len=120) :: file_obj
    character(len=20) :: num_str

    doubleprecision, dimension(:,:), allocatable :: mat

    narg = iargc()
    if (narg == 0) then
        write(*,*) 'transfer gradient and force (sign of values)'
        write(*,*) 'example: force_grad filename num_atom'
    elseif (narg == 2) then
        call getarg(1, file_obj)
        call getarg(2, num_str)
        read(num_str,*) num_atom
        allocate(mat(3,num_atom))
        open(20, file=trim(file_obj), action='read')
        do i =1, num_atom
            read(20, *) mat(:,i)
        enddo
        mat = mat * (-1.0)
        close(20)

        open(21, file=trim(file_obj), status='replace')
        do i =1, num_atom
            write(21, '(3D20.12)') mat(:,i)
        enddo
        close(21)
        deallocate(mat)
    endif
endprogram