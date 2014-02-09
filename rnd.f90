module rnd_lib
 implicit none
 interface random
  module procedure d1rnd,z1rnd,d2rnd,z2rnd,d3rnd,z3rnd
 end interface
contains
 
 subroutine arnd()
  implicit none
  integer s,clock,crate,cmax
  integer,allocatable :: seed(:)
  call random_seed(size=s)
  allocate(seed(s))
  call random_seed(get=seed)
  write(*,*) 'oldseed: ',seed
  call system_clock(clock,crate,cmax)
  !write(*,*)clock,crate,cmax
  seed(1)=clock 
  call random_seed(put=seed)
  write(*,*) 'newseed: ',seed
 end subroutine 

 double precision function drnd( )
  implicit none
  call random_number(drnd)
 return
 end function

 subroutine d1rnd(d)
  double precision,intent(out)  :: d(:)
  call random_number(d)
 end subroutine
 subroutine z1rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:)
  character(len=*),parameter :: subnam='z1rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine

 subroutine d2rnd(d)
  double precision,intent(out)  :: d(:,:)
  call random_number(d)
 end subroutine
 subroutine z2rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:,:)
  character(len=*),parameter :: subnam='z2rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine

 subroutine d3rnd(d)
  double precision,intent(out)  :: d(:,:,:)
  call random_number(d)
 end subroutine
 subroutine z3rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:,:,:)
  character(len=*),parameter :: subnam='z3rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine
 

 integer function irnd( maxi )
  integer,intent(in) :: maxi
  double precision :: d
  call random_number(d)
  irnd=int(d*maxi)+1
 return
 end function
 subroutine irand(maxi,ix)
  integer,intent(in)  :: maxi
  integer,intent(out) :: ix(:)
  integer :: i,n
  double precision,allocatable :: d(:)
  n=size(ix)
  allocate(d(n))
  call random_number(d)
  ix=int(d*maxi)+1
  deallocate(d)
 return
 end subroutine
 

 subroutine lottery2(npnt,m,n,wcol,wrow,points)
  implicit none
  integer,intent(in) :: npnt,m,n
  double precision,intent(in) :: wcol(m),wrow(n)
  integer,intent(out) :: points(npnt,2)
  character(len=*),parameter :: subnam='lottery2'
  double precision,allocatable :: pcol(:),prow(:),d(:,:)
  double precision :: scol,srow
  integer :: ipnt,i,j,info
  double precision,external :: dasum
  allocate(pcol(0:m),prow(0:n),d(npnt,2),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  scol=dasum(m,wcol,1); srow=dasum(n,wrow,1)
  pcol(0)=0.d0; do i=1,m;pcol(i)=pcol(i-1)+dabs(wcol(i))/scol;enddo
  prow(0)=0.d0; do j=1,n;prow(j)=prow(j-1)+dabs(wrow(j))/srow;enddo
  call random_number(d)
  do ipnt=1,npnt
   points(ipnt,1)=find_d(m+1,pcol,d(ipnt,1)); if(points(ipnt,1).gt.m)points(ipnt,1)=m
   points(ipnt,2)=find_d(n+1,prow,d(ipnt,2)); if(points(ipnt,2).gt.n)points(ipnt,2)=n
  end do
  deallocate(pcol,prow,d)
 end subroutine
 
 pure integer function find_d(n,x,y) result (pos)
  ! for sorted vector x(1) <= x(2) <= ... <= x(n) and value y find pos, s.t. x(pos) <= y < x(pos+1)
  implicit none
  integer,intent(in) :: n
  double precision,intent(in) :: x(n),y
  integer :: s,t,i
  logical :: key
  if(n.eq.0)then;pos=0;return;endif
  if(y.lt.x(1))then;pos=0;return;endif
  if(x(n).le.y)then;pos=n;return;endif
  s=1;t=n;pos=(t+s)/2
  do while(t-s.gt.1)
   if(y.lt.x(pos))then;t=pos;else;s=pos;end if
   pos=(s+t)/2
  enddo 
  return
 end function
end module
