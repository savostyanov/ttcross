module default_lib
 implicit none
 interface default
  module procedure default_i, default_d, default_z, default_a
 end interface
contains
 integer function default_i(def,opt) result(fun)
  integer,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
 double precision function default_d(def,opt) result(fun)
  double precision,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
 double complex function default_z(def,opt) result(fun)
  double complex,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
 character function default_a(def,opt) result(fun)
  character,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
 

 subroutine share(first,last,own)
  implicit none
  include 'mpif.h'
  integer,intent(in) :: first,last
  integer,intent(out) :: own(0:)
  integer :: p,me,nproc,info
  
  ! check mpi variables
  call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_size fail: ',info;stop;endif
  call mpi_comm_rank(MPI_COMM_WORLD,me,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_rank fail: ',info;stop;endif
  
  ! share load: proc p owns tasks own(p)...own(p+1)-1
  own(0)=first
  do p=1,nproc-1; own(p) = first + int(dble(last-first+1)*dble(p)/nproc); enddo
  own(nproc) = last+1
 end subroutine
end module
