program main
 use time_lib
 use quad_lib
 use mat_lib
 implicit none
 include "mpif.h"
 integer :: i,j,m,n,q,r,piv,decay,a,info,nproc,me,neval
 double precision :: f,bnd,t1,t2,tcrs, einf,efro,ainf,afro, mem,logdet, val,tru,stddev
 double precision,allocatable::par(:),mat(:,:)
 double precision,parameter :: tpi=6.28318530717958647692528676655900577d0
 integer,parameter :: NCHK=2**10
 integer,parameter :: fileid=2
 character(len=32) :: aa
 double precision,external :: dfunc_gauss_cont

 ! Read params
 call getarg(1,aa); read(aa,'(i10)')m; if(m.eq.0)m=8         ! dimension of the problem
 call getarg(2,aa); read(aa,'(i10)')q; if(q.eq.0)q=22        ! QMC level
 call getarg(3,aa); read(aa,'(i10)')r; if(r.eq.0)r=4         ! QMC repeats

 call mpi_init(info)
 if(info.ne.0)then;write(*,*)'mpi: init fail: ',info;stop;endif
 call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
 if(info.ne.0)then;write(*,*)'mpi: comm_size fail: ',info;stop;endif
 call mpi_comm_rank(MPI_COMM_WORLD,me,info)
 if(info.ne.0)then;write(*,*)'mpi: comm_rank fail: ',info;stop;endif
 !write(*,'(a,i3,a,i3)')'mpi: I am ',me,' of ',nproc
 
 if(me.eq.0)then
  write(*,'(a)') 'Hi, this is QMC integrating generalised Gaussian...'
  write(*,'(3x,a,i10)') 'dimension  :',m
  write(*,'(3x,a,i10)') 'QMC level  :',q
  write(*,'(3x,a,i10)') 'QMC repeats:',r
  write(*,'(3x,a,i10)') 'MPI procs  :',nproc
!!$OMP PARALLEL
!  if (omp_get_thread_num().eq.0) then
!   write(*,'(3x,a,i10)')'OMP threads:', omp_get_num_threads()
!  end if 
!!$OMP END PARALLEL
 end if

 allocate(par(2+m*m+m),mat(m,m), stat=info)
 if(info.ne.0)then;write(*,*)'cannot allocate par:',info;stop;endif
 par(1)=dble(121)
 bnd=10.d0
 par(2)=bnd
!  call eye(mat)
 !call dscal(m*m,2.d0,mat,1)
 !call laplace(mat); call matinv(mat)
 forall(i=1:m,j=1:m)mat(i,j)=dexp(-dble(abs(i-j)))
 call dcopy(m*m,mat,1,par(3),1)
 call dpotrf('l',m,mat,m,info)
 if(info.ne.0)then;write(*,*)'dportf: info: ',info; stop; endif
 logdet=0.d0
 do i=1,m
  logdet=logdet+2*dlog(dabs(mat(i,i)))
 end do

 t1=timef()
 call qmc(m, q,r, val, stddev, dfunc_gauss_cont, par)
 ! Scale the integrand for QMC
 val = val*((2*bnd)**m)
 t2=timef()
 tcrs=t2-t1
 
 if(me.eq.0)write(*,'(a,i20,a,e12.4,a)') '...with ',(2_8**q)*r,' evaluations completed in ',tcrs,' sec.'

 if(me.eq.0) then
  tru = dexp( dlog(tpi)*m/2 - logdet/2)
  !write(*,*) dexp(logdet)
  write(*,'(a,e30.22)') 'computed value    :', val
  write(*,'(a,e30.22)') 'standard deviation:', stddev
  write(*,'(a,e30.22)') 'analytic value    :', tru
  write(*,'(a,f7.2)')   'correct digits    :', -dlog(dabs(1.d0-val/tru))/dlog(10.d0)

  open(UNIT=fileid,FILE='qmcresults.dat',position="append")
  write(fileid,'(i5,i5,i5,i10,f25.20,2x,e10.3,2x,i5)') decay, m, nproc, neval, val, tcrs, r
  close(fileid)
  write(*,'(a)')'Good bye.'
 endif
 call mpi_finalize(info)
 if(info.ne.0)then;write(*,*)'mpi: finalize fail: ',info;stop;endif
end program

double precision function dfunc_gauss_cont(m,y,par) result(f)
 ! function to compute tensor entry
  implicit none
  integer,intent(in) :: m
  double precision,intent(in) :: y(m)
  double precision,intent(inout),optional :: par(*)
  double precision,parameter :: tpi=6.28318530717958647692528676655900577d0
  integer :: matrix,i
  double precision :: bnd
  double precision,allocatable :: ay(:)
  allocate(ay(m))
  matrix=3
  bnd=par(2)
  call dgemv('n',m,m,1.d0,par(matrix),m,y,1,0.d0,ay,1)
  f=dot_product(y,ay)
  f = f*(bnd**2)  ! scale [-1,1] -> [-bnd,bnd]
  f=dexp(-f/2)
  deallocate(ay)
end function
