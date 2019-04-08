ttcross
=======

Cross interpolation of high-dimensional arrays in tensor train format.
This code implements the parallel version of the tensor cross interpolation algorithm, see
  * Dmitry Savostyanov, [Quasioptimality of maximum-volume cross interpolation of tensors](http://dx.doi.org/10.1016/j.laa.2014.06.006), Linear Algebra and its Applications, 2014.
  * Sergey Dolgov and Dmitry Savostyanov, [Parallel cross interpolation for high-precision calculation of high-dimensional integrals](http://arxiv.org/abs/1903.11554), arXiv: 1903.11554, 2019.


Usage
-----

  0. Read the articles.
  1. Download the code.
  2. Edit the Makefile to match your system.
  3. Compile the code by `make` command. Use `make ising` or `make gauss` to compile specific group of tests.
  **Note:** this software makes use of BLAS and LAPACK libraries. Please set the ***BLASLIB*** variable in the *Makefile* accordingly.
  4. Run appropriate `./test_*` command to execute an experiment.
     Use `mpirun -np NN ./test_*` to run the experiment in parallel on NN processors.
     Use `OMP` system variables (e.g. `OMP_NUM_THREADS`) to control the number of OMP threads.

Tests
-----

   * `test_mc_ising`:  calculate Ising susceptibility integrals using Monte Carlo algorithm
   * `test_qmc_ising`: calculate Ising susceptibility integrals using quasi Monte Carlo algorithm
   * `test_crs_ising`: calculate Ising susceptibility integrals using tensor cross interpolation algorithm
   * `test_mpf_ising`: calculate Ising susceptibility integrals using tensor cross interpolation algorithm in multiple precision
 **Note:** the mpf test uses the GNU MPFR Library [https://www.mpfr.org/](https://www.mpfr.org/) Please install it on your system and set the ***MPFLIB*** variable in the *Makefile* accordingly.
   * `test_qmc_gauss`: calculate generalised Gaussian integrals using quasi Monte Carlo algorithm
   * `test_crs_gauss`: calculate generalised Gaussian integrals using tensor cross interpolation algorithm

If you benefit from this code, please consider citing our papers in your work.
