program hess2freq

    use read_hess
    ! use write_freq

    implicit none

    integer :: Natom, info, i, j
    real*8, allocatable :: hessisian(:,:), normal_mode(:,:)
    real*8, allocatable :: force_constant(:,:), force_constant_tmp(:,:)
    real*8, allocatable :: frequency(:), atom_mass(:)
    real*8, allocatable :: work(:), half_mass(:,:)
    

    call get_Natom(Natom)

    allocate(hessisian(3*Natom,3*Natom))
    allocate(normal_mode(3*Natom,3*Natom))
    allocate(frequency(3*Natom))
    allocate(atom_mass(Natom))
    allocate(half_mass(3*Natom,3*Natom))
    allocate(work(100000))

    call read_orca_atom_mass(Natom, atom_mass)

    ! do i = 1, Natom
    !     print *, atom_mass(i)
    ! end do

    do i = 1, 3*Natom
        do j = 1, 3*Natom
            half_mass(i, j) = 0.0d0                   
        end do
    end do 

    do i = 1, Natom
        do j = 1, 3
            half_mass(3*(i-1)+j, 3*(i-1)+j) =  half_mass(3*(i-1)+j, 3*(i-1)+j) + 1/sqrt(atom_mass(i))
        end do
    end do

    ! print * ,

    ! do i = 1, 3*Natom
    !     print *, half_mass(i, i)
    ! end do


    call  read_orca_hess(Natom, hessisian)

    allocate(force_constant(3*Natom,3*Natom))
    allocate(force_constant_tmp(3*Natom,3*Natom))


    ! call dgemm('N', 'N', 54, 54, 54, 1.0d0, half_mass, 54, hessisian, 54, 0.0d0, force_constant_tmp, 54)

    ! call dgemm('N', 'N', 54, 54, 54, 1.0d0, force_constant_tmp, 54, half_mass, 54, 0.0d0, force_constant, 54)

    do i = 1, 3*Natom
        do j = 1, 3*Natom
            force_constant(i, j) = hessisian(i, j)*half_mass(i, i)*half_mass(j, j)
        end do
    end do


    normal_mode = force_constant

    call dsyev('V', 'L', 3*Natom, normal_mode, 3*Natom, frequency, work, 100000, info)

    ! do i = 1, 3*Natom
    !     print *, sqrt(abs(frequency(i)*4.35974434d-18/(0.529177249*0.529177249d-20*1.660538d-27)))/(2*3.1415926*2.9979d10)
    ! end do

    do i = 1, 3*Natom
        print *, normal_mode(i,7), normal_mode(i,8), normal_mode(i,9)
    end do


    deallocate(work)
    deallocate(force_constant)
    deallocate(force_constant_tmp)
    deallocate(half_mass)
    deallocate(atom_mass)
    deallocate(hessisian)
    deallocate(frequency)
    deallocate(normal_mode)
    


end program hess2freq
