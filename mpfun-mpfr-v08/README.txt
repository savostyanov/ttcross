MPFUN2015: A thread-safe arbitrary precision package

Revision date: 20 August 2018

AUTHOR:
David H. Bailey
Lawrence Berkeley National Lab (retired) and University of California, Davis
Email: dhbailey@lbl.gov
 
COPYRIGHT AND DISCLAIMER:
All software in this package (c) 2018 David H. Bailey. By downloading or using this software you agree to the copyright, disclaimer and license agreement in the accompanying file DISCLAIMER.txt.


I. PURPOSE OF PACKAGE:

This system permits one to perform floating-point computations (real and complex) to arbitrarily high numeric precision, by making only relatively minor changes to existing Fortran-90 programs. All basic arithmetic operations and transcendental functions are supported, together with several special functions.

In addition to fast execution times, one key feature of this package is a 100% THREAD-SAFE design, which means that user-level applications can be easily converted for parallel execution, say using a threaded parallel environment such as OpenMP. There are no global shared variables (except static compile-time data), and no initialization is necessary unless extremely high precision is required.

There are two versions of this package:

MPFUN-Fort: This is the all-Fortran version. It compiles in just a few seconds on any system with a Fortran-2003 compliant compiler.

MPFUN-MPFR: This is similar to MPFUN-Fort in its user interface, but it calls the MPFR package for all low-level functions and operations. The MPFUN-MPFR version is significantly faster on most operations than MPFUN-Fort, but installation is more complicated (because the GMP and MPFR packages must first be installed).

Most users may wish to first try the MPFUN-Fort package, and then switch to the MPFUN-MPFR package only if needed for best performance.


II. DOCUMENTATION:

A detailed description of this software, with instructions for writing Fortran code to use the package, is available in this technical paper:
 
David H. Bailey, "MPFUN2015: A thread-safe arbitrary precision package," 
http://www.davidhbailey.com/dhbpapers/mpfun2015.pdf.


III. INSTALLING COMPILERS:

Installation, compilation and linking is relatively straightforward, provided that you have a Unix-based system, such as Apple OSX or Linux, and a command-line interface such as the Terminal application of Apple OSX.

For Apple OSX systems, you first must install the "Command Line Tools" package, which is available (for free) from the Apple Developer website. Login (or register, if this your first access), using your Apple ID at: https://developer.apple.com/download/more. Then select "Command Line Tools" for your particular version of the MAC OSX operating system. Install the downloaded package on your system.

The gfortran compiler, which is highly recommended for this package, is available for a variety of systems at the website https://gcc.gnu.org/wiki/GFortranBinaries.

On a Mac, when you attempt to install gfortran, you may get the message "gfortran.pkg can't be opened because it is from an unidentified developer". If so, go to the Security page of your Mac's System Preferences, and click "Open anyway".

The packages should also work with IBM's xlf_r compiler, Intel's ifort compiler and Portland Group's pgf90 (in the case of MPFUN-MPFR, the corresponding C compilers are also required).


IV. DOWNLOADING MPFUN-Fort AND MPFUN-MPFR:

From the website http://www.davidhbailey.com/dhbsoftware, download the file "mpfun-fort-vnn.tar.gz" or "mpfun-mpfr-vnn.tar.gz, whichever version you wish you use (replace "vnn" by whatever is the current version on the MPFUN website, such as "v22"). If the file is not decompressed by your browser, use gunzip at the shell level to do this. Some browsers (such as the Apple OSX Safari browser) do not drop the ".gz" suffix after decompression; if so, remove this suffix manually at the shell level.

Then type
  tar xfv mpfun-fort-vnn.tar 
or
  tar xfv mpfun-mpfr-vnn.tar
  
whichever version you are using (where again "vnn" is replaced by whatever is the current version). This should create the directory and unpack all files.


V. INSTALLING MPFUN-MPFR (if you use MPFUN-Fort, skip to section VI below):

To use MPFUN-MPFR, you must first install the GMP and MPFR packages, as follows (this presumes you have a command-line interface and gcc installed):
  
1. Using a browser (e.g., Apple Safari or Firefox, download the file "gmp-6.0.0a.tar.xz" (or whatever is latest version) from the site https://gmplib.org. The .xz file is more likely to be decompressed by your browser than the .lz file. Move this file to a suitable spot on your system, typically to the Documents folder.

2. Open the Terminal application or an equivalent command-line interface. Type "which gcc" to see if /usr/local/bin is in your default search path. If not, create a file .bashrc or the equivalent in your home directory with the line "PATH=/usr/bin:/usr/local/bin:$PATH". Then type

  source .bashrc

3. In the Documents folder (or wherever the tar file was moved), type 

  tar -xf gmp-6.0.0a.tar.xz

(or whatever is the latest version of GMP). This should create the directory "gmp-6.0.0" (or a similar name).  

4. Change directory to this GMP directory, then type "./configure", followed by "make", then "make check". All tests should pass. Then type "make install". On Apple systems and some others, you may need to type instead "sudo make install",  which will request your computer system's admin password. This should place several files, including "libgmp.10.dylib", "libgma.a", "libgmp.dylib" and "libgmp.la" in /usr/local/lib.

5. Using a browser, download "mpfr-3.1.3.tar.xz" (or whatever is the latest version) from http://www.mpfr.org/mpfr-current/ and move it to a suitable spot in your Documents folder (or wherever else is convenient).

6. In the Documents folder (or wherever the tar file was moved), type

  tar -xf mpfr-3.1.3.tar.xz

(or similar name).  Then change directory to mpfr-3.1.3 (or similar name) and type

  ./configure

followed by "make", "make check" (see if all tests pass or skipped), and then either "make install" or "sudo make install", as appropriate for your system. This should place the files "libmpfr.4.dylib", "libmpfr.a", "libmpfr.dylib" and "libmpfr.la" in /usr/local/lib.  

7. To test the installation, place the C program "sample" (beginning with the line "#include <stdio.h>") from the URL http://www.mpfr.org/sample.html into a file "sample.c" (located anywhere within the Documents folder), then compile by typing

  gcc -o sample sample.c -lmpfr -lgmp

Then when you type "./sample", you should see the single line of output given at the bottom of the URL http://www.mpfr.org/sample.html.
  
Note that both GMP and MPFR require 5-10 minutes to install as described above.


VI. COMPILING MPFUN-Fort and/or MPFUN-MPFR:

For both MPFUN-Fort and MPFUN-MPFR, there are actually two variants of the software, both of which are included in the distribution file:

Variant 1: This is recommended for basic applications that do not dynamically change the precision level (or do so only rarely).

Variant 2: This is recommended for more sophisticated applications that dynamically change the precision level.

The two variants of the packages correspond to two variants of module MPFUNG, the high-level language interface module. Compile/link scripts are available in the fortran directory for the gfortran compiler and other supported compilers. For the gfortran and ifort compilers, which support the real*16 datatype, the respective compile-link scripts include the proper modules.

For example, to compile variant 1 of either the MPFUN-Fort or MPFUN-MPFR library with the GNU gfortran compiler, go to the fortran directory and type

  ./gnu-complib1.scr

Then to compile and link the application program prog.f90 for variant 1, producing the executable file prog, type

  ./gnu-complink1.scr prog

NOTE: The very first time you compile the library (using either the complib1.scr or complib2.scr scripts), you may see "fatal" errors, such as various modules not found. This is normal -- just repeat the library compile scripts. The library compile scripts invoke the compiler twice for this reason.

Several test programs are included in the fortran directory of the packages, together with output files -- see Section VIII below.


VII. BRIEF SUMMARY OF CODING INSTRUCTIONS AND USAGE:

Here is a brief summary of Fortran coding instructions. For full details, see:

David H. Bailey, "MPFUN2015: A thread-safe arbitrary precision package," 
http://www.davidhbailey.com/dhbpapers/mpfun2015.pdf

First set the parameter mpipl, the "default" precision level in digits, which is the maximum precision level to be used for subsequent computation, and is used to specify the amount of storage required for multiprecision data. mpipl is set in a parameter statement at the start of module MPFUNF, which is in file mpfunf.f90 in the fortran directory. In the code as distributed, mpipl is set to 1200 digits (sufficient to run each of the test programs), but it can be set to any level greater than or equal to 32 digits. mpipl is automatically converted to mantissa words by the formula 

  mpwds = int (mpipl / mpdpw + 2),

where mpdpw is a system parameter, and where int () means truncate to integer. For MPFUN-Fort, mpdpw is log_{10} (2^{48}) = 14.44943979187..., whereas for MPFUN-MPFR it is log_{10}(2^{64}) = 19.26591972249... (both values are double precision approximations). The resulting parameter mpwds is the internal default precision level, in words. All computations are performed to mpwds precision unless the user, within an application code, specifies a lower precision level.

After setting the value of mpipl in module MPFUNF, compile the appropriate version of the library, using one of the scripts mentioned above.

Next, place the following line in every subprogram of the user's application code that contains a multiprecision variable or array, at the beginning of the declaration section, before any implicit or type statements:

  use mpmodule

To designate a variable or array as multiprecision real (MPR) in your application code, use the Fortran-90 type statement with the type "mp_real", as in this example:

  type (mp_real) a, b(m), c(m,n)

Similarly, to designate a variable or array as multiprecision complex (MPC), use a type statement with "mp_complex".

Thereafter when one of these variables or arrays appears in code, e.g.,

  d = a + b(i) * sqrt(3.d0 - c(i,j))

the proper multiprecision routines are automatically called by the Fortran compiler.

Most common mixed-mode combinations (arithmetic operations, comparisons and assignments) involving MPR, MPC, double precision (DP) and integer arguments are supported. A complete list of supported mixed-mode operations is given in the documentation paper.

However, there are some hazards. For example, the code

  r1 = 3.14159d0

where r1 is MPR, does NOT produce the true multiprecision equivalent of 3.14159, unless the numerical value is a modest-sized whole number or exact binary fraction. In fact, the software will flag such usage with a run-time error. To obtain the full MP converted value, write this as

r1 = "3.14159d0"

instead. Similarly, the code

  r2 = r1 + 3.d0 * sqrt (2.d0)

where r1 and r2 are MPR, does NOT produce the true multiprecision value, since the expression 3.d0 * sqrt (2.d0) will be performed in double precision (according to standard Fortran-90 precedence rules). In fact, as above, the above line of code will result in a run-time error. To obtain the fully accurate r2 result, write this as

  r2 = r1 + 3.d0 * sqrt (mpreal (2.q0))

See documentation for details.

Input/output of MP variables or array elements is done using the subroutines mpread and mpwrite. See documentation for details.

Standard Fortran intrinsics are supported with MPR and MPC arguments, and they operate similarly to the standard double precision (DP) and double complex (DC) equivalents. A complete list of supported functions and subroutines is given in the documentation paper. 

Many user applications do not need to change the working precision from the initially-defined default level, whereas more sophisticated applications may need to change this frequently. As noted above, for both MPFUN-Fort and MPFUN-MPFR, there are two variants of the language interface module MPFUNG:

Variant 1: This is recommended for basic applications that do not dynamically change the precision level (or do so only rarely).

Variant 2: This is recommended for more sophisticated applications that dynamically change the precision level.

See documentation for full details on the differences between these two variants.


VIII. SAMPLE APPLICATION PROGRAMS:

The current release of the software includes these application programs in the fortran directory (exactly the same programs are in both the mpfun-fort/fortran directory and the mpfun-mpfr/fortran directory):

testmpfun.f90    Tests most arithmetic and transcendental functions.

tpphix3.f90      Performs a 3-level Poisson polynomial application.

tpslq1.f90       Performs the standard 1-level PSLQ integer relation algorithm.

tpslqm1.f90      Performs the 1-level multipair PSLQ integer relation algorithm.

tpslqm2.f90      Performs the 2-level multipair PSLQ integer relation algorithm.

tpslqm3.f90      Performs the 3-level multipair PSLQ integer relation algorithm.

tquad.f90        Evaluates integrals using tanh-sinh, exp-sinh and sinh-sinh algorithms.

tquadgs.f90      Evaluates integrals using Gaussian quadrature.

The fortran directory also includes the corresponding output files (e.g., tpphix3.ref.txt) for each of these programs. It also includes the script "mpfun-test.scr", which runs each of these sample programs above except tquadgs.f90, which takes considerably more run time. To run this script, first edit the file "mpfun-test.scr", replacing "complink2.scr" with, say, "gnu-complink2.scr" or whatever is one's preferred compiler script. Then simply type

./mpfun-tests.scr

If, after running this script (several minutes run time), the results in the output files match those in the reference output files (except for timings, etc.), then one can be fairly confident that the software is working properly. Full descriptions of these programs are included in the documentation paper mentioned above.


