.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o
      
FC      = gfortran
LDR     = gfortran
OPT     = -O2 -ffree-line-length-none
LIB     = -L/usr/local/lib  -llapack -lblas


DEPS    = timef, default, rnd, lr, ptype, tt, ttind, ttio, dmrgg

OBJ    = $(DEPS:,=.o).o main.o

MAIN = main

main:       $(OBJ)
		$(LDR) $(OBJ) -o main $(LIB)

.f.o:
		$(FC) $(OPT) -c $<
.f90.o:
		$(FC) $(OPT) -c $<

clean:
		rm -f *.o *.obj *.mod main
