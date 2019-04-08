.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o
      
FC      = mpif90 -ffree-line-length-none
CC      = mpicc
LDR     = mpif90
OPT     = -O2 -fopenmp #-fdefault-real-8 # enable for quad prec
BLASLIB = -llapack -lblas # Set this to your BLAS & LAPACK distribution
MPFLIB  =    # Set this to your MPFR distribution,
             # for example -I/home/sd901/opt/mpfr/include -L/home/sd901/opt/mpfr/lib on my PC

MPD    = ./mpfun-mpfr-v08
SRC    = zero, nan, trans, default, timef, say, rnd, ptype, ort, lr, mat, quad, tt, ttaux, ttind, ttio, dmrgg, qmc, mc
MPF    = mpfuna, mpfunf, mpfung1, mpinterface, mpmodule, mpblas, ttmp, dmrggmp
OBJ    = $(SRC:,=.o).o
MPOBJ  = $(MPF:,=.o).o

MAIN = main

all: ising gauss

ising:   test_crs_ising.exe test_qmc_ising.exe test_mpf_ising.exe test_mc_ising.exe
gauss:   test_crs_gauss.exe test_qmc_gauss.exe 
main2014: main2014.exe

mpfuna.o:   $(MPD)/mpfuna.f90
		$(FC) $(OPT) -fno-underscoring -c $< $(MPFLIB)
mpfunf.o:   $(MPD)/mpfunf.f90
		$(FC) $(OPT) -fno-underscoring -c $< $(MPFLIB)
mpfung1.o:  $(MPD)/mpfung1.f90
		$(FC) $(OPT) -fno-underscoring -c $< $(MPFLIB)
mpmodule.o: $(MPD)/mpmodule.f90
		$(FC) $(OPT) -fno-underscoring -c $< $(MPFLIB)
mpinterface.o: $(MPD)/mpinterface.c
		$(FC) $(OPT) -c $< $(MPFLIB)
.f.o:
		$(FC) $(OPT) -c $<
.f90.o:
		$(FC) $(OPT) -c $<
.F90.o:
		$(FC) $(OPT) -c $<
.c.o:
		$(CC) $(OPT) -c $<

test_mpf_ising.exe: $(OBJ) $(MPOBJ) test_mpf_ising.o 
	$(LDR) $(OPT) $^ -o $@ -lmpfr $(BLASLIB) $(MPFLIB)
%.exe:  $(OBJ)  %.o
	$(LDR) $(OPT) $^ -o $@ $(BLASLIB)

clean:
		rm -f *.o *.obj *.mod *.exe



