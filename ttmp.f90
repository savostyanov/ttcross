module ttmp_lib
 use mpmodule
 implicit none
 integer,parameter :: tt_size=2048
 
 type,public :: pointi
  integer,dimension(:),pointer,contiguous :: p=>null()
 end type
 type,public :: pointi2
  integer,dimension(:,:),pointer,contiguous :: p=>null()
 end type
 type,public :: point_mp
  type(mp_real),dimension(:),pointer,contiguous :: p=>null()
 end type
 type,public :: point_mp2
  type(mp_real),dimension(:,:),pointer,contiguous :: p=>null()
 end type
 type,public :: point_mp3
  type(mp_real),dimension(:,:,:),pointer,contiguous :: p=>null()
 end type
 
 type,public:: mptt           ! multiple precision tensor train
  integer :: l=1              ! index of the leftmost core
  integer :: m=0              ! index of the rightmost core
  integer :: n(tt_size)=0     ! mode sizes (storage for matrices)
  integer :: q(tt_size)=0     ! first mode sizes (for matrices)
  integer :: s(tt_size)=0     ! second mode sizes (for matrices)
  integer :: t=0              ! data type for future sparsity etc.
  integer :: r(0:tt_size)=0   ! TT ranks
  type(point_mp3) :: u(tt_size) ! TT cores
 end type
 
 interface alloc        ! allocation and deallocation of the cores
  module procedure mptt_alloc
 end interface
 interface dealloc
  module procedure mptt_dealloc
 end interface
 interface erank        ! effective rank of tt-vector
  module procedure mptt_rank
 end interface
 interface ones
  module procedure mptt_ones
 end interface
 interface say
  module procedure mptt_say
 end interface
 interface assignment (=)
  module procedure mptt_assign
 end interface

contains

 subroutine mptt_alloc(arg)
  implicit none
  type(mptt),intent(inout) :: arg
  character(len=*),parameter :: subnam='mptt_alloc'
  integer :: i,j,k,p,info
  type(mp_real) :: zero
  zero="0.d0"
  if(arg%m.lt.arg%l)return
  if(arg%l.le.0)then;write(*,*)subnam,': %l should be > 0';stop;endif
  if(arg%m.gt.tt_size)then;write(*,*)subnam,': %m exceeds tt_size, change parameter and recompile!';stop;endif
  do p=arg%l,arg%m
   if(associated(arg%u(p)%p))deallocate(arg%u(p)%p)
   allocate(arg%u(p)%p(arg%r(p-1),arg%n(p),arg%r(p)), stat=info)
   do k=1,arg%r(p); do j=1,arg%n(p); do i=1,arg%r(p-1); arg%u(p)%p(i,j,k)=zero; enddo;enddo;enddo
   if(info.ne.0)then;write(*,*)'TT allocate fail: no memory';stop;endif
  end do
 end subroutine
 subroutine mptt_dealloc(arg)
  implicit none
  type(mptt),intent(inout) :: arg
  character(len=*),parameter :: subnam='mptt_dealloc'
  integer :: i
  do i=1,tt_size
   if(associated(arg%u(i)%p))deallocate(arg%u(i)%p)
  end do
 end subroutine
 subroutine mptt_assign(b,a)
  implicit none
  type(mptt),intent(inout) :: b
  type(mptt),intent(in) :: a
  integer :: k,l,m
  l=a%l;m=a%m
  b%l=l; b%m=m; b%n(l:m)=a%n(l:m); b%r(l-1:m)=a%r(l-1:m); call alloc(b)
  do k=l,m;call mpcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(k)%p,1); end do
 end subroutine
 subroutine mptt_ones(arg)
 ![tt] array of ones
  implicit none
  type(mptt),intent(inout) :: arg
  character(len=*),parameter :: subnam='mptt_ones'
  integer :: l,m,k
  type(mp_real) :: one
  one="1.d0"
  l=arg%l;m=arg%m
  if(m.lt.l)return
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  arg%r(l-1:m)=1
  call alloc(arg)
  do k=l,m; arg%u(k)%p=one; end do
 end subroutine 
 
 double precision function mptt_rank(arg) result (r)
  implicit none
  type(mptt),intent(in) :: arg
  integer :: l,m,i,a,b,d
  l=arg%l;m=arg%m;d=m-l+1
  if(d.le.0)then;r=-1.d0;return;endif
  if(d.eq.1)then;r= 0.d0;return;endif
  r=0.d0
  do i=arg%l,arg%m
   r=r+arg%r(i-1)*arg%n(i)*arg%r(i)
  end do
  if(r.eq.0.d0)return
  b=arg%r(l-1)*arg%n(l) + arg%n(m)*arg%r(m)
  if(d.eq.2)then;r=r/b;return;endif
  a=sum(arg%n(l+1:m-1))
  r=(dsqrt(b*b+4.d0*a*r)-b)/(2.d0*a)
  return
 end function
 subroutine mptt_say(arg)
  implicit none
  type(mptt),intent(in) :: arg
  character(len=1)d
  write(*,'(a,i2,a,i4,a,f6.2)') 'dtt[',arg%l,':', arg%m,']: rank ',erank(arg)
  if(all(arg%r(arg%l:arg%m).le.100))then;d='3';else;d='4';endif
  if(all(arg%n(arg%l+1:arg%m)==arg%n(arg%l)))then
   write(*,'(a,1x,i'//d//',a)') 'n: ',arg%n(arg%l),' for all modes'
  else 
   write(*,'(a,1x,1024i'//d//')') 'n: ',arg%n(arg%l:arg%m)
  end if
  write(*,'(a,1024i'//d//')') 'r: ',arg%r(arg%l-1:arg%m)
 end subroutine

end module
