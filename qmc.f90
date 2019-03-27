subroutine qmc(d, l,r, F, stddev, dfunc_c, par)
 implicit none
 include "mpif.h"
 interface
   function dfunc_c(m,y,par) result (func_result)
     integer,intent(in) :: m
     double precision,intent(in) :: y(m)
     double precision,intent(inout),optional :: par(*)
     double precision :: func_result
   end function dfunc_c
 end interface
integer, intent(in) :: d, l,r
double precision, intent(out) :: F, stddev
double precision,intent(inout),optional :: par(*)

integer(kind=8) :: i,N
integer(kind=8),parameter :: one=1_8,two=2_8
integer :: nproc,me,info,j,k
integer :: z(d) 
double precision :: y(d), delta(d,r), Fall(r)
integer, parameter :: fileid=2

call mpi_comm_size(MPI_COMM_WORLD, nproc, info)
call mpi_comm_rank(MPI_COMM_WORLD, me, info)

N = two**l

if (me.eq.0) then
!  delta = 0d0
!  call random_seed()
  call random_number(delta)
!  print *, 'delta = ', delta
  open(UNIT=fileid,FILE='lattice.dat')
   read(fileid,*) z
  close(fileid)
end if
call mpi_bcast(delta,d*r,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
if(info.ne.0)then;write(*,*)'mpi bcast fail for delta: ',info;stop;endif
call mpi_bcast(z,d,MPI_INTEGER,0,MPI_COMM_WORLD,info)
if(info.ne.0)then;write(*,*)'mpi bcast fail for lattice: ',info;stop;endif

Fall=0.d0
do k=1,r
 F = 0d0
 do i=one*me,N-one,one*nproc
  do j=1,d
    y(j) = dble(z(j))*dble(i)/dble(N) + delta(j,k)
    y(j) = y(j) - int(y(j))
    y(j) = 2*y(j)-1.d0
  end do
  F = F + dfunc_c(d, y, par)
 end do
 F = F/N
 Fall(k) = F
end do 
call mpi_allreduce(MPI_IN_PLACE, Fall, r, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
if(info.ne.0)then;write(*,*)'mpi allreduce fail: ',info;stop;endif
F = sum(Fall)/r
stddev = sqrt(sum((Fall-F)**2)/(r-1))/F

end subroutine

