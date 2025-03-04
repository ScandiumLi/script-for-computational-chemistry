	program COVERT XYZ TO ARC
      implicit double precision (a-h, o-z) 
      character elem*2, line*100
	parameter(natom=3,nstep=24)
	dimension elem(nstep,natom), rx(nstep,natom),ry(nstep,natom) 
      dimension rz(nstep,natom)

	open(01,file = 'au3-pos-1.xyz')
	   do j = 1, nstep 
	   print*, j
	   read(01,*)
	   read(01,*) 
	   do i = 1, natom
	      read(01,*) elem(j,i), rx(j,i), ry(j,i), rz(j,i)
	print*, i
	   enddo
	   enddo
	close(01)

      open(02,file = 'au3.arc') 
        write(02,'(a)')  '!BIOSYM archive 3'                                                              
	  write(02,'(a)')  'PBC=OFF'   
        do nrec = 1, nstep                                                                    
           write(02,'(a)') '                       -11571.7839'
           write(02,'(a)')     '!DATE     Jun 07 18:06:48 2008'     
           do i = 1, natom 
              WRITE (02, '(A,3X,3f15.9,21X,A,f7.3)') elem(nrec,i),
     +      rx(nrec,i), ry(nrec,i),
     +      rz(nrec,i),elem(nrec,i),0.000  
           enddo


           write(02,'(A)') 'end'
           write(02,'(A)') 'end'
        enddo
      close(02)

      end


