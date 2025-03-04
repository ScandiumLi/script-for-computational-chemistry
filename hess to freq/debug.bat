gfortran.exe -c .\read_hess.f90
gfortran.exe -c .\hess2freq.f90
gfortran.exe .\read_hess.o .\hess2freq.o -o main.exe
.\main.exe