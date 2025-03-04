program fch2hess

    implicit none

    integer :: N_atom, i, j, mod_number
    character(len=80) :: char_tmp, trash, label_atom_number, fchk_file, hess_file
    character(len=30) :: label_coord_start, label_coord_end
    character(len=30) :: label_hess_start, label_hess_end
    character(len=30) :: label_mass_start, label_mass_end
    character(len=30) :: label_atomid_start, label_atomid_end
    integer, allocatable :: atom_number(:)
    character(len=30), allocatable :: atom_symbol(:)
    real*8, allocatable :: Hessian(:,:), atom_coord(:,:), atom_mass(:)

    ! N_atom = 18 ! 原子数，按实际情况修改
    fchk_file = "S0-freq.fchk"
    hess_file = "S0-free.hess"

    label_atom_number = "Number of atoms                            I"
    label_atomid_start = "Atomic numbers"
    label_atomid_end = "Nuclear charges"
    label_mass_start = "Real atomic weights"
    label_mass_end = "Atom fragment info"
    label_coord_start = "Current cartesian coordinates"
    label_coord_end = "Number of symbols in /Mol/"
    label_hess_start = "Cartesian Force Constants"
    label_hess_end = "Nonadiabatic coupling"

    open(unit=17, file=fchk_file, status="old", action="read")

    read(17, *) 
    read(17, *)
    read(17, '(A)') char_tmp
    read(char_tmp(46:), "(I3)") N_atom

    ! do i = 1, 20
    !     read(17, '(a)') char_tmp
    !     print *, char_tmp
    
    !     if (index(char_tmp, label_atom_number) .ne. 0) then
    !         print *, char_tmp
    !         read(char_tmp, "(A45, I3)") trash, N_atom
    !     else
    !         cycle
    !     end if
    ! end do

    ! print *, N_atom
    close(17)

    allocate(atom_number(N_atom))
    allocate(atom_symbol(N_atom))
    allocate(atom_mass(N_atom))
    allocate(Hessian(3*N_atom, 3*N_atom))
    allocate(atom_coord(N_atom, 3))

    open(unit=18, file=fchk_file, status="old", action="read")

    do 
        read(18, '(a)') char_tmp
        ! print *, char_tmp
    
        if (index(char_tmp, label_atomid_start) .ne. 0) then
            ! print *, char_tmp
            read(18, "(6(I12))") (atom_number(i), i = 1, N_atom)
        else if (index(char_tmp, label_atomid_end) .ne. 0) then   
            exit
        else
            cycle
        end if
    end do

    close(18)

    open(unit=19, file=fchk_file, status="old", action="read")

    do 
        read(19, '(a)') char_tmp
        ! print *, char_tmp
    
        if (index(char_tmp, label_coord_start) .ne. 0) then
            ! print *, char_tmp
            read(19, "(5(1PE16.8))") ((atom_coord(i, j), j = 1, 3), i = 1, N_atom)
        else if (index(char_tmp, label_coord_end) .ne. 0) then   
            exit
        else
            cycle
        end if
    end do

    close(19)

    open(unit=20, file=fchk_file, status="old", action="read")

    do 

        read(20, '(a)') char_tmp
        ! print *, char_tmp
    
        if (index(char_tmp, label_hess_start) .ne. 0) then
            ! print *, char_tmp
            read(20, "(5(1PE16.8))") ((Hessian(i, j), j = 1, i), i = 1, 3*N_atom)
        else if (index(char_tmp, label_hess_end) .ne. 0) then   
            exit
        else
            cycle
        end if
    end do

    close(20)

    open(unit=21, file=fchk_file, status="old", action="read")

    do 

        read(21, '(a)') char_tmp
        ! print *, char_tmp
    
        if (index(char_tmp, label_mass_start) .ne. 0) then
            ! print *, char_tmp
            read(21, "(5(1PE16.8))") (atom_mass(i), i = 1, N_atom)
        else if (index(char_tmp, label_mass_end) .ne. 0) then   
            exit
        else
            cycle
        end if
    end do

    close(21)

    do i = 1, N_atom

        if (atom_number(i) == 6) then 
            atom_symbol(i) = "C "
        else if (atom_number(i) == 7) then
            atom_symbol(i) = "N "
        else if (atom_number(i) == 8) then
            atom_symbol(i) = "O "
        else if (atom_number(i) == 1) then
            atom_symbol(i) = "H "
        else if (atom_number(i) == 79) then
            atom_symbol(i) = "Au"
        else if (atom_number(i) == 47) then
            atom_symbol(i) = "Ag"
        else if (atom_number(i) == 29) then
            atom_symbol(i) = "Cu"
        else if (atom_number(i) == 53) then
            atom_symbol(i) = "I "
        else if (atom_number(i) == 16) then
            atom_symbol(i) = "S "
        else if (atom_number(i) == 15) then
            atom_symbol(i) = "P "
        else if (atom_number(i) == 17) then
            atom_symbol(i) = "Cl"
        end if
        
    end do

    Hessian = Hessian + transpose(Hessian)

    do i = 1, 3*N_atom
        Hessian(i, i) = Hessian(i, i)/2.0d0
        ! print *, Hessian(i,1)
    end do

    mod_number = mod(3*N_atom, 5)

    open(unit=22, file=hess_file, status="new", action="write")

    write(22, *) 
    write(22, "(A18)") "$orca_hessian_file"
    write(22, *) 
    write(22, "(A9)") "$act_atom"
    write(22, "(A1)") "0"
    write(22, *)
    write(22, "(A10)") "$act_coord"
    write(22, *) "0"
    write(22, *)
    write(22, "(A11)") "$act_energy"
    write(22, *) "        0.000000"
    write(22, *)
    write(22, "(A8)") "$hessian"
    write(22, *) 3 * N_atom
    do i = 1, 3*N_atom/5
        write(22, '(3X, 5I19)') 5*(i-1), 5*(i-1)+1, 5*(i-1)+2, 5*(i-1)+3, 5*(i-1)+4
        do j = 1, 3*N_atom
            write(22, '(I6, 3X, 5E19.10)') j-1, hessian(j, 5*(i-1)+1), hessian(j, 5*(i-1)+2), hessian(j, 5*(i-1)+3), &
            hessian(j, 5*(i-1)+4), hessian(j, 5*(i-1)+5)
        end do
    end do

    write(22, '(3X, 4I19)') (3*N_atom-mod_number+i-1, i = 1, mod_number)
    do j = 1, 3*N_atom
        write(22, '(I6, 3X, 4E19.10)') j-1, (hessian(j, 3*N_atom-mod_number+i), i = 1, mod_number)
    end do

    write(22, *) 
    write(22, "(A6)") "$atoms"
    write(22, "(I3)") N_atom

    do i = 1, N_atom
        write(22, "(1X, A2, 3X, F9.5, 1X, 3F19.12)") atom_symbol(i), atom_mass(i), atom_coord(i, 1), &
        atom_coord(i, 2), atom_coord(i, 3)
    end do

    write(22, *)

    close(22)

    deallocate(atom_number)
    deallocate(atom_symbol)
    deallocate(atom_mass)
    deallocate(Hessian)
    deallocate(atom_coord)

end program fch2hess
