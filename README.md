ttcross
=======

Cross interpolation of high-dimensional arrays in tensor train format.
This code implements the parallel version of the tensor cross interpolation algorithm, see
  * Dmitry Savostyanov, [Quasioptimality of maximum-volume cross interpolation of tensors](http://dx.doi.org/10.1016/j.laa.2014.06.006), Linear Algebra and its Applications, 2014.
  * Sergey Dolgov and Dmitry Savostyanov, [Parallel cross interpolation for high-precision calculation of high-dimensional integrals](http://arxiv.org/abs/1903.11554), arXiv: 1903.11554, 2019.

If you benefit from this code, please consider citing our papers in your work.

Basic usage
-----------

  0. Read the articles.
  1. Download the code.
  2. Edit the Makefile to match your system.
  **Note:** this software makes use of BLAS and LAPACK libraries. Please set the ***BLASLIB*** variable in the *Makefile* accordingly.
  3. Compile the code by `make` command.
  4. Run appropriate `./test_*` command to execute an experiment.
     Use `mpirun -np NN ./test_*` to run the experiment in parallel on NN processors.
     Use `OMP` system variables (e.g. `OMP_NUM_THREADS`) to control the number of OMP threads.

Basic tests
-----------
   
   * `test_mc_ising.exe`:  calculate Ising susceptibility integrals using Monte Carlo algorithm
   * `test_qmc_ising.exe`: calculate Ising susceptibility integrals using quasi Monte Carlo algorithm
   * `test_crs_ising.exe`: calculate Ising susceptibility integrals using tensor cross interpolation algorithm

Quadruple precision
-------------------
   Double precision calculations are accurate to approximately 15 decimal digits. To reach higher precision, you can compile the  code enabling `-fdefault-real-8` or `-fdefault-real-16` in the compiler.
   Note that you also will need to obtain and link quadruple-precision BLAS/LAPACK libraries.

Multiple precision
------------------

   For those who desires even higher precision, we provide a multi-precision version of the code. It relies on the MPFUN2015 library written by David H. Bailey. 
   The source code of MPFUN2015 is included in the `mpfun-mpfr-v08` directory. The code uses GNU MPFR Library [https://www.mpfr.org/](https://www.mpfr.org/). Please install it on your system and set the ***MPFLIB*** variable in the *Makefile* accordingly.
   Compile the multi-precision version with `make test_mpf_ising.exe`, noting the following:
   * BLAS/LAPACK libraries are *not* needed for this experiment --- all necessary subroutines are provided in the `mpblas.f90` file.
   * The `-fdefault-real-XX` flag has to be disabled.
   * The precision is adjusted via the parameter `mpipl` in `./mpfun-mpfr-v08/mpfunf.f90`.
   After compiling, run the `test_mpf_ising.exe` to calculate Ising susceptibility integrals using tensor cross interpolation algorithm in multiple precision.

