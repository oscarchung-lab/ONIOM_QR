!this program is to generate A B conformer based on PDB file
program AB_Conformer
    implicit none
    integer :: narg, num, i, isend, pos
    character(len=120) :: filename, line
    character(len=30) :: word, res1, res2
    logical :: alive, b

    narg=iargc()
    if (narg == 1) then
        call getarg(1, filename)

        inquire(file=trim(filename),exist=alive)
        if (.not. alive) then
            write(*,*) trim(filename)//' does not exist'
            stop
        endif

        pos=index(trim(filename), '.')

        open(10, file=trim(filename),action='read')
        open(11, file=filename(1:pos-1)//'_a.pdb',status='replace')
        open(12, file=filename(1:pos-1)//'_b.pdb',status='replace')

        read(10,'(A)') line
        read(line, *) word
        b = .false.
        do while (.not. b)
            write(11,'(A)') trim(line)
            write(12,'(A)') trim(line)

            read(10,'(A)') line
            read(line, *) word
            b=(trim(word) == 'ATOM')
        enddo

        backspace(10)
        num=0
        res1=line(23:26)
        do while (.true.)
            read(10,'(A)',iostat=isend) line
            if (isend<0) then
                exit
            endif
            read(line, *) word
            if (trim(word)=='ATOM' .or. trim(word)=='ANISOU' .or. trim(word)=='HETATM') then
                if (line(17:17) == 'A') then
                    line(17:17)=' '
                    line(57:60)='1.00'
                    write(11,'(A)') trim(line)
                    res2=line(23:26)
                    if(trim(res2)/=trim(res1)) then
                        num=num+1  
                        res1=trim(res2)              
                    endif
                elseif(line(17:17) == 'B') then
                    line(17:17)=' '
                    line(57:60)='1.00'
                    write(12,'(A)') trim(line)
                else
                    write(11,'(A)') trim(line)
                    write(12,'(A)') trim(line)
                endif
            elseif(trim(word)=='CONECT' .or. trim(word)=='MASTER') then
                cycle
            else
                write(11,'(A)') trim(line)
                write(12,'(A)') trim(line)
            endif
        enddo
        close(11)
        close(12)
        if (num==0) then
            write(*,*) 'this is the unique pdb structure'
            call system('rm -f '//filename(1:pos-1)//'_a.pdb '//filename(1:pos-1)//'_b.pdb')
        else
            write(*,*) 'there are ', num, 'confomer residues in the pdb'
        endif
    else
        write(*,*) "This program is to generate A B conformer based on PDB file"
        write(*,*) "AB_Conformer xxx.pdb"
    endif

endprogram