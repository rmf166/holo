# GNU Fortran Compiler
F90    = gfortran
LIBS   = -llapack
FFLAGS = -O0 -g -fconvert=big-endian -fno-automatic -fbacktrace \
         -ffpe-trap=invalid,zero,overflow -fPIC -fcheck=all -Wall -Wextra \
         -fimplicit-none -Wuninitialized -pedantic # -Warray-temporaries

## Intel Fortran Compiler
#F90    = ifort
#LIBS   =
#FFLAGS = -g -C -traceback -mkl

fa.exe : main.f90
	$(F90) $(FFLAGS) $(LIBS) -o fa.exe main.f90

clean:
	rm -f fa.exe *.mod *.o
