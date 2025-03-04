program xyz2arc

    implicit none

    integer, parameter :: natom = 61, nstep = 20000, n_skip = 100 
    ! natom为单个结构原子数，nstep为轨迹数, n_skip=100为每100步xyz输出1步arc
    integer :: i, j
    character(len=3) :: elem(nstep, natom)
	real*8 :: rx(nstep, natom),ry(nstep, natom), rz(nstep, natom) 

    open(unit=20, file="Au25SH18-MD-pos-1.xyz", status="old", action="read")
    do i = 1, nstep
        read(20, *)
        read(20, *)
        do j = 1, natom
            read(20, *) elem(i, j), rx(i, j), ry(i, j), rz(i, j)
        end do     
    end do

    close(20)

    open(21, file="Au25SH18-MD-pos-1.arc", status="new", action="write")
   
    write(21, '(a)')  '!BIOSYM archive 3' 
    write(21, '(a)')  'PBC=OFF'

    do i = 1, nstep
        if (mod(i, n_skip) == 0) then
            write(21, '(a)') '                       -11571.7839'
            write(21, '(a)')     '!DATE     Jun 07 18:06:48 2008'    
            do j = 1, natom
                write(21, '(A,3X,3f15.9,21X,A,f7.3)') elem(i, j), rx(i, j), ry(i, j), rz(i, j), elem(i, j), 0.000
            end do
            write(21,'(A)') 'end'
            write(21,'(A)') 'end'
        else
            cycle
        end if       
    end do

    close(21)

end program xyz2arc
