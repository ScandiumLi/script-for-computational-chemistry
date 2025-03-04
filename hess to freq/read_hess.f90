module read_hess

    implicit none

    contains


    subroutine get_Natom(Natom)
        integer :: Natom, i
        character(len=12) :: atom_start

        open(unit=18, file="Azulene-S0-opt.hess", status="old", action="read")

            do i = 1, 3000
            read (18, '(a6)') atom_start

            if (trim(atom_start) == "$atoms") then 
                ! print *, "yes"   
                read (18, *) Natom
                ! print * , Natom
                exit
            else 
                read (18, *) 
                cycle   
            end if
        end do

        close(18)
        
    end subroutine get_Natom

    subroutine read_orca_atom_mass(Natom, atom_mass)

        integer :: Natom, i, j
        real*8 :: atom_mass(Natom)
        character(len=12) :: atom_start
        character(len=3) :: atom_symbol(Natom)

        open(unit=19, file="Azulene-S0-opt.hess", status="old", action="read")

        do i = 1, 3000
            read (19, '(a6)') atom_start

            if (trim(atom_start) == "$atoms") then 
                ! print *, "yes"   
                read (19, *) 
                do j = 1, Natom
                    read (19, '(A3, f12.5)') atom_symbol(j), atom_mass(j)
                    ! print *, atom_mass(j)
                end do
                exit
            else 
                read (19, *) 
                cycle   
            end if
        end do

        close(19)
             
    end subroutine read_orca_atom_mass

    subroutine read_orca_hess(Natom, hessisian)

        integer :: Natom, i, j, trash, mod_number
        real*8 :: hessisian(3*Natom,3*Natom)
       
        mod_number = mod(3*Natom, 5)

        ! print *, mod_number
        
        open(unit=20, file="Azulene-S0-opt.hess", status="old", action="read")
        
        do i = 1, 15
            read(20, *) 
        end do

        do i = 1, 3*Natom/5

            do j = 1, 3*Natom
                read(20, '(I8, 5E19.10)') trash, hessisian(j, 5*(i-1)+1), hessisian(j, 5*(i-1)+2), hessisian(j, 5*(i-1)+3), &
                hessisian(j, 5*(i-1)+4), hessisian(j, 5*(i-1)+5)
            end do

            read(20, *) 
        
        end do
        
        do j = 1, 3*Natom
            read(20, *) trash, (hessisian(j, 3*Natom-mod_number+i), i = 1, mod_number)     
        end do

    
        close(20)
        
    end subroutine read_orca_hess
    

end module read_hess
