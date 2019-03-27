module ttaux_lib
 use tt_lib
 implicit none
 double precision,parameter,private :: tpi=6.28318530717958647692528676655900577d0
 
 interface exp
  module procedure dtt_exp,ztt_exp
 end interface
 interface ones
  module procedure dtt_ones, ztt_ones
 end interface
 interface zeros
  module procedure dtt_zeros, ztt_zeros
 end interface
 interface unit
  module procedure dtt_unit_i,dtt_unit_ii,dtt_unit_d,dtt_unit_d1
  module procedure ztt_unit_i,ztt_unit_ii,ztt_unit_d,ztt_unit_d1
 end interface
 interface random
  module procedure dtt_rnd, ztt_rnd
 end interface
 interface mirror
  module procedure dtt_mirror, ztt_mirror
 end interface
 interface reverse
  module procedure dtt_rev, ztt_rev
 end interface
 interface double
  module procedure dtt_double, ztt_double
 end interface
 interface push
  module procedure dtt_push
 end interface 

contains

! ONES ZEROS
 subroutine dtt_ones(arg)
 ![tt] array of ones
  implicit none
  type(dtt),intent(inout) :: arg
  character(len=*),parameter :: subnam='dtt_ones'
  integer :: l,m,k
  l=arg%l;m=arg%m
  if(m.lt.l)return
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  arg%r(l-1:m)=1
  call alloc(arg)
  do k=l,m; arg%u(k)%p=1.d0; end do
 end subroutine 
 subroutine ztt_ones(arg)
 ![tt] array of ones
  implicit none
  type(ztt),intent(inout) :: arg
  character(len=*),parameter :: subnam='ztt_ones'
  integer :: l,m,k
  l=arg%l;m=arg%m
  if(m.lt.l)return
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  arg%r(l-1:m)=1
  call alloc(arg)
  do k=l,m; arg%u(k)%p=(1.d0,0.d0); end do
 end subroutine
 
 subroutine dtt_zeros(arg)
 ![tt] array of zeros
  implicit none
  type(dtt),intent(inout) :: arg
  character(len=*),parameter :: subnam='dtt_zeros'
  integer :: l,m,k
  l=arg%l;m=arg%m
  if(m.lt.l)return
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  arg%r(l-1:m)=1
  call alloc(arg)
  do k=l,m; arg%u(k)%p=0.d0; end do
 end subroutine 
 subroutine ztt_zeros(arg)
 ![tt] array of zeros
  implicit none
  type(ztt),intent(inout) :: arg
  character(len=*),parameter :: subnam='ztt_zeros'
  integer :: l,m,k
  l=arg%l;m=arg%m
  if(m.lt.l)return
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  arg%r(l-1:m)=1
  call alloc(arg)
  do k=l,m; arg%u(k)%p=(0.d0,0.d0); end do
 end subroutine

! UNIT
 subroutine dtt_unit_i(arg,i)
  ! unit element in position given by 1 <= i <= 2^d
  implicit none
  type(dtt),intent(inout) :: arg
  integer :: i
  character(len=*),parameter :: subnam='dtt_unit_i'
  integer :: l,m,k,ind(tt_size)
  l=arg%l;m=arg%m
  if(m.lt.l)then;write(*,*)subnam,': set arg%l and arg%m before calling'; return; endif
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  call mindex(i,arg%n(l:m),ind(1:m-l+1))
  call dtt_unit_ii(arg,ind)
 end subroutine 
 subroutine ztt_unit_i(arg,i)
  implicit none
  type(ztt),intent(inout) :: arg
  integer :: i
  character(len=*),parameter :: subnam='ztt_unit_i'
  integer :: l,m,k,ind(tt_size)
  l=arg%l;m=arg%m
  if(m.lt.l)then;write(*,*)subnam,': set arg%l and arg%m before calling'; return; endif
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  call mindex(i,arg%n(l:m),ind(1:m-l+1))
  call ztt_unit_ii(arg,ind)
 end subroutine 

 subroutine dtt_unit_ii(arg,ind)
  ! unit element in position given by ind(1:m-l+1)
  implicit none
  type(dtt),intent(inout) :: arg
  integer :: ind(:)
  character(len=*),parameter :: subnam='dtt_unit_ii'
  integer :: l,m,k
  l=arg%l;m=arg%m
  if(m.lt.l)then;write(*,*)subnam,': set arg%l and arg%m before calling'; return; endif
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  if(any(ind(1:m-l+1)<=0).or.any(ind(1:m-l+1)>arg%n(l:m)))then;write(*,*)subnam,': wrong index: ',ind;stop;endif
  arg%r(l-1:m)=1; call alloc(arg)
  do k=l,m; arg%u(k)%p=0.d0; arg%u(k)%p(1,ind(k-l+1),1)=1.d0; end do
 end subroutine 
 subroutine ztt_unit_ii(arg,ind)
  implicit none
  type(ztt),intent(inout) :: arg
  integer :: ind(:)
  character(len=*),parameter :: subnam='ztt_unit_ii'
  integer :: l,m,k
  l=arg%l;m=arg%m
  if(m.lt.l)then;write(*,*)subnam,': set arg%l and arg%m before calling'; return; endif
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  if(any(ind(1:m-l+1)<=0).or.any(ind(1:m-l+1)>arg%n(l:m)))then;write(*,*)subnam,': wrong index: ',ind;stop;endif
  arg%r(l-1:m)=1; call alloc(arg)
  do k=l,m; arg%u(k)%p=(0.d0,0.d0); arg%u(k)%p(1,ind(k-l+1),1)=(1.d0,0.d0); end do
 end subroutine
 
 subroutine dtt_unit_d(arg,x)
  implicit none
  type(dtt),intent(inout) :: arg
  double precision,intent(in) :: x
  double precision :: xx(1)
  xx=x; call dtt_unit_d1(arg,xx)
 end subroutine
 subroutine ztt_unit_d(arg,x)
  implicit none
  type(ztt),intent(inout) :: arg
  double precision,intent(in) :: x
  double precision :: xx(1)
  xx=x; call ztt_unit_d1(arg,xx)
 end subroutine

 subroutine dtt_unit_d1(arg,x)
  implicit none
  type(dtt),intent(inout) :: arg
  double precision,intent(in) :: x(:)
  character(len=*),parameter :: subnam='dtt_unit_d1'
  integer :: l,m,r(0:tt_size),n(tt_size), id,dd,i,j,k,pos,mm,ind(tt_size)=0
  double precision :: xx
  l=arg%l; m=arg%m; r=arg%r; n=arg%n; dd=size(x)
  if(m.lt.l)then;write(*,*)subnam,': set arg%l and arg%m before calling'; return; endif
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  if(l.gt.m)return
  arg%r(l-1:m)=1; call alloc(arg)
  mm=(m-l+1)/dd
  do id=1,dd
   xx=x(id)
   if(xx.lt.0.d0)return
   if(xx.gt.1.d0)xx=xx-int(xx)
   do j=1,mm
    pos=l+(id-1)*mm+mm-j
    i=int(n(pos)*xx)
    if(i.eq.n(pos))i=n(pos)-1
    ind(pos-l+1)=i+1
    xx=xx*n(pos)-i
   end do
  end do
  do k=l,m; arg%u(k)%p=0.d0; arg%u(k)%p(1,ind(k-l+1),1)=1.d0; end do
 end subroutine
 subroutine ztt_unit_d1(arg,x)
  implicit none
  type(ztt),intent(inout) :: arg
  double precision,intent(in) :: x(:)
  character(len=*),parameter :: subnam='ztt_unit_d1'
  integer :: l,m,r(0:tt_size),n(tt_size), id,dd,i,j,k,pos,mm,ind(tt_size)=0
  double precision :: xx
  l=arg%l; m=arg%m; r=arg%r; n=arg%n; dd=size(x)
  if(m.lt.l)then;write(*,*)subnam,': set arg%l and arg%m before calling'; return; endif
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  if(l.gt.m)return
  arg%r(l-1:m)=1; call alloc(arg)
  mm=(m-l+1)/dd
  do id=1,dd
   xx=x(id)
   if(xx.lt.0.d0)return
   if(xx.gt.1.d0)xx=xx-int(xx)
   do j=1,mm
    pos=l+(id-1)*mm+mm-j
    i=int(n(pos)*xx)
    if(i.eq.n(pos))i=n(pos)-1
    ind(pos-l+1)=i+1
    xx=xx*n(pos)-i
   end do
  end do
  do k=l,m; arg%u(k)%p=(0.d0,0.d0); arg%u(k)%p(1,ind(k-l+1),1)=(1.d0,0.d0); end do
 end subroutine
 

! RND
 subroutine dtt_rnd(arg)
  ![tt] random array with prescribed sizes and ranks
  use rnd_lib
  implicit none
  type(dtt),intent(inout) :: arg
  character(len=*),parameter :: subnam='dtt_rnd'
  integer :: l,m,k
  l=arg%l;m=arg%m
  if(m.lt.l)return
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  if(.not.all(arg%r(l-1:m)>0))then;write(*,*)subnam,': r: ',arg%r(l-1:m);stop;endif
  call alloc(arg)
  do k=l,m
   call random(arg%u(k)%p)
  end do
 end subroutine 
 subroutine ztt_rnd(arg)
  ![tt] random array with prescribed sizes and ranks
  use rnd_lib
  implicit none
  type(ztt),intent(inout) :: arg
  character(len=*),parameter :: subnam='ztt_rnd'
  integer :: l,m,k
  l=arg%l;m=arg%m
  if(m.lt.l)return
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  if(.not.all(arg%r(l-1:m)>0))then;write(*,*)subnam,': r: ',arg%r(l-1:m);stop;endif
  call alloc(arg)
  do k=l,m
   call random(arg%u(k)%p)
  end do
 end subroutine

!EXP SIN COS
 subroutine ztt_exp(arg,alpha,beta)
  ![bitt] exp function exp(ax+b)
  implicit none
  type(ztt),intent(inout) :: arg
  double complex,intent(in),optional :: alpha,beta
  double complex,parameter :: one=(1.d0,0.d0),zero=(0.d0,0.d0),im1=(0.d0,1.d0)
  double complex :: alp,bet
  integer :: r,k,l,m
  if(present(alpha))then;alp=alpha;else;alp=one;endif
  if(present(beta))then;bet=beta;else;bet=zero;endif
  m=arg%m; l=arg%l
  arg%n(l:m)=2;arg%r(l-1:m)=1
  call alloc(arg)
  arg%u(l)%p(1,1,1)=cdexp(bet)
  arg%u(l)%p(1,2,1)=cdexp(bet+alp)
  do k=l+1,m
   alp=alp*2
   arg%u(k)%p(1,1,1)=one
   arg%u(k)%p(1,2,1)=cdexp(alp)
  end do
 end subroutine
 subroutine dtt_exp(arg,alpha,beta)
  ![bitt] exp function exp(ax+b)
  implicit none
  type(dtt),intent(inout) :: arg
  double precision,intent(in),optional :: alpha,beta
  double precision :: alp,bet
  integer :: r,k,l,m
  if(present(alpha))then;alp=alpha;else;alp=1.d0;endif
  if(present(beta))then;bet=beta;else;bet=0.d0;endif
  m=arg%m; l=arg%l
  arg%n(l:m)=2;arg%r(l-1:m)=1
  call alloc(arg)
  arg%u(l)%p(1,1,1)=dexp(bet)
  arg%u(l)%p(1,2,1)=dexp(bet+alp)
  do k=l+1,m
   alp=alp*2
   arg%u(k)%p(1,1,1)=1.d0
   arg%u(k)%p(1,2,1)=dexp(alp)
  end do
 end subroutine

 subroutine dtt_sin(arg,alpha,beta)
  ![bitt] sine function
  implicit none
  type(dtt),intent(inout) :: arg
  double precision,intent(in),optional :: alpha,beta
  double precision :: alp,bet
  if(present(alpha))then;alp=alpha;else;alp=1.d0;endif
  if(present(beta))then;bet=beta;else;bet=0.d0;endif
  call dtt_cos(arg,alp,bet-tpi/4)
 end subroutine
 subroutine dtt_cos(arg,alpha,beta)
  ![bitt] cosine function
  implicit none
  type(dtt),intent(inout) :: arg
  double precision,intent(in),optional :: alpha,beta
  double precision :: alp,bet,ca,sa
  integer :: r,k,l,m
  if(present(alpha))then;alp=alpha;else;alp=1.d0;endif
  if(present(beta))then;bet=beta;else;bet=0.d0;endif
  m=arg%m; l=arg%l
  arg%n(l:m)=2;  arg%r(l-1)=1;arg%r(l:m-1)=2;arg%r(m)=1;
  call alloc(arg)
  arg%u(l)%p(1,:,:)=reshape( (/dcos(bet), dcos(alp+bet), -dsin(bet), -dsin(alp+bet) /), [2,2] )
  do k=l+1,m-1
   alp=alp*2; ca=dcos(alp); sa=dsin(alp)
   arg%u(k)%p=reshape( (/1.d0, 0.d0, ca, sa, 0.d0, 1.d0, -sa, ca /), [2,2,2] )
  end do
  alp=alp*2
  arg%u(m)%p(:,:,1)=reshape( (/1.d0, 0.d0, dcos(alp), dsin(alp) /), [2,2])
 end subroutine
 subroutine dtt_cossin(arg,alpha,beta)
  ![bitt] cosine and sine function, r(d)=2
  implicit none
  type(dtt),intent(inout) :: arg
  double precision,intent(in),optional :: alpha,beta
  double precision :: alp,bet,ca,sa
  integer :: r,k,l,m
  if(present(alpha))then;alp=alpha;else;alp=1.d0;endif
  if(present(beta))then;bet=beta;else;bet=0.d0;endif
  m=arg%m; l=arg%l
  arg%n(l:m)=2;  arg%r(l-1)=1;arg%r(l:m)=2
  call alloc(arg)
  arg%u(l)%p(1,:,:)=reshape( (/dcos(bet), dcos(alp+bet), -dsin(bet), -dsin(alp+bet) /), [2,2] )
  do k=l+1,m
   alp=alp*2; ca=dcos(alp); sa=dsin(alp)
   arg%u(k)%p=reshape( (/1.d0, 0.d0, ca, sa, 0.d0, 1.d0, -sa, ca /), [2,2,2] )
  end do
  arg%u(m)%p(:,:,2)=-arg%u(m)%p(:,:,2)
 end subroutine

! MIRROR
 subroutine dtt_mirror(arg,part)
  implicit none
  type(dtt),intent(inout) :: arg
  integer,intent(in),optional :: part(2)
  integer :: l,m,p,i,j,k,r(0:tt_size),n(tt_size)
  l=arg%l;m=arg%m;r=arg%r;n=arg%n
  if(present(part))then;l=part(1);m=part(2);endif
  do p=l,m
   forall(i=1:r(p-1),j=1:n(p)/2,k=1:r(p))
    ! swap arg%u(p)%p(i,j,k) <-> arg%u(p)%p(i,n(p)-j+1,k)
    arg%u(p)%p(i,j,k)=arg%u(p)%p(i,j,k)+arg%u(p)%p(i,n(p)-j+1,k)
    arg%u(p)%p(i,n(p)-j+1,k)=arg%u(p)%p(i,j,k)-arg%u(p)%p(i,n(p)-j+1,k)
    arg%u(p)%p(i,j,k)=arg%u(p)%p(i,j,k)-arg%u(p)%p(i,n(p)-j+1,k)
   end forall
  end do
 end subroutine
 subroutine ztt_mirror(arg,part)
  implicit none
  type(ztt),intent(inout) :: arg
  integer,intent(in),optional :: part(2)
  integer :: l,m,p,i,j,k,r(0:tt_size),n(tt_size)
  l=arg%l;m=arg%m;r=arg%r;n=arg%n
  if(present(part))then;l=part(1);m=part(2);endif
  do p=l,m
   forall(i=1:r(p-1),j=1:n(p)/2,k=1:r(p))
    ! swap arg%u(p)%p(i,j,k) <-> arg%u(p)%p(i,n(p)-j+1,k)
    arg%u(p)%p(i,j,k)=arg%u(p)%p(i,j,k)+arg%u(p)%p(i,n(p)-j+1,k)
    arg%u(p)%p(i,n(p)-j+1,k)=arg%u(p)%p(i,j,k)-arg%u(p)%p(i,n(p)-j+1,k)
    arg%u(p)%p(i,j,k)=arg%u(p)%p(i,j,k)-arg%u(p)%p(i,n(p)-j+1,k)
   end forall
  end do
 end subroutine

! REVERSE
 subroutine dtt_rev(arg)
  implicit none
  type(dtt),intent(inout) :: arg
  type(dtt) :: tmp
  integer :: l,m,k
  l=arg%l;m=arg%m
  if(l.gt.m)return
  tmp%l=l;tmp%m=m
  forall(k=l-1:m)tmp%r(k)=arg%r(m-k+l-1)
  forall(k=l:m)tmp%n(k)=arg%n(m-k+l)
  call alloc(tmp)
  do k=l,m; call trans(4,arg%u(k)%p,tmp%u(m-k+l)%p); end do
  arg=tmp
  call dealloc(tmp)
 end subroutine 
 subroutine ztt_rev(arg)
  implicit none
  type(ztt),intent(inout) :: arg
  type(ztt) :: tmp
  integer :: l,m,k
  l=arg%l;m=arg%m
  if(l.gt.m)return
  tmp%l=l;tmp%m=m
  forall(k=l-1:m)tmp%r(k)=arg%r(m-k+l-1)
  forall(k=l:m)tmp%n(k)=arg%n(m-k+l)
  call alloc(tmp)
  do k=l,m; call trans(4,arg%u(k)%p,tmp%u(m-k+l)%p); end do
  arg=tmp
  call dealloc(tmp)
 end subroutine 

! DOUBLE
 subroutine dtt_double(arg)
  implicit none
  ! a -> [a a]
  type(dtt),intent(inout) :: arg
  character(len=*),parameter :: subnam='dtt_double'
  integer :: l,m,r,i,j,k,info
  l=arg%l; m=arg%m; r=arg%r(m)
  if(associated(arg%u(m+1)%p))then
   !write(*,*)subnam,': next core is already allocated'
   deallocate(arg%u(m+1)%p)
  end if
  allocate(arg%u(m+1)%p(r,2,r),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': allocation fail';stop;endif
  forall(i=1:r,j=1:2,k=1:r)arg%u(m+1)%p(i,j,k)=0.d0
  forall(i=1:r,j=1:2      )arg%u(m+1)%p(i,j,i)=1.d0
  arg%m=m+1; arg%n(m+1)=2; arg%r(m+1)=r
 end subroutine
 subroutine ztt_double(arg)
  implicit none
  ! a -> [a a]
  type(ztt),intent(inout) :: arg
  character(len=*),parameter :: subnam='ztt_double'
  integer :: l,m,r,i,j,k,info
  l=arg%l; m=arg%m; r=arg%r(m)
  if(associated(arg%u(m+1)%p))then
   !write(*,*)subnam,': next core is already allocated'
   deallocate(arg%u(m+1)%p)
  end if
  allocate(arg%u(m+1)%p(r,2,r),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': allocation fail';stop;endif
  forall(i=1:r,j=1:2,k=1:r)arg%u(m+1)%p(i,j,k)=(0.d0,0.d0)
  forall(i=1:r,j=1:2      )arg%u(m+1)%p(i,j,i)=(1.d0,0.d0)
  arg%m=m+1; arg%n(m+1)=2; arg%r(m+1)=r
 end subroutine


! AP
 subroutine dtt_ap(arg,alpha,beta)
  ![bitt] arg(i) = alpha*i + beta, i=0...2^d-1
  implicit none
  type (dtt),intent(inout) :: arg
  double precision,intent(in),optional :: alpha,beta
  character(len=*),parameter :: subnam ='dtt_ap'
  double precision :: alp,bet
  integer :: l,m,k
  alp=1.d0;if(present(alpha))alp=alpha
  bet=0.d0;if(present(beta))bet=beta
  l=arg%l;m=arg%m
  if(m.lt.l)return
  arg%r(l-1)=1;arg%r(l:m-1)=2;arg%r(m)=1;arg%n(l:m)=2
  call alloc(arg)
  arg%u(l)%p=reshape( (/2.d0,2.d0,bet,alp+bet/) ,[1,2,2])
  do k=l+1,m-1
   arg%u(k)%p=reshape( (/2.d0,0.d0,2.d0,0.d0,0.d0,1.d0,alp,1.d0/) ,[2,2,2])
  end do
  arg%u(m)%p=reshape( (/0.d0,1.d0,alp,1.d0/) ,[2,2,1])
 end subroutine
 
 subroutine dtt_push(arg,dir,new)
  implicit none
  type(dtt),intent(inout),target :: arg
  character(len=1),intent(in),optional :: dir
  double precision,intent(in),optional :: new
  character(len=*),parameter :: subnam='dtt_push'
  type(dtt) :: tmp
  double precision :: val
  integer :: l,m,k,p,i,j,d,info,ind(tt_size)
  integer,pointer :: r(:),n(:)
  double precision,allocatable :: x(:),y(:)
  tmp%l=arg%l; tmp%m=arg%m; l=arg%l; m=arg%m; if(l.gt.m)return
  if(arg%r(l-1).ne.1)then;write(*,*)subnam,': left border rank should be 1';stop;endif
  if(arg%r(m).ne.1)then;write(*,*)subnam,': right border rank should be 1';stop;endif
  !if(.not.all(arg%n(l:m)==2))then;write(*,*)subnam,': not bitt on input';stop;endif
  d=+1
  if(present(dir))then
   select case(dir)
    case('r','R');d=+1
    case('l','L');d=-1
    case default;write(*,*)subnam,': illegal dir: ',dir;stop
   end select
  end if 
  if(present(new))then
   val=new
  else
   select case(d)
    case(+1);ind(l:m)=arg%n(l:m);val=tijk(arg,ind(l:m))
    case(-1);ind(l:m)=1;val=tijk(arg,ind(l:m))
    case default;write(*,*)subnam,': illegal d:',d;stop
   end select
  end if 
  r=>arg%r; n=>arg%n
  tmp%r(l:m-1)=r(l:m-1)+1;tmp%r(l-1)=r(l-1);tmp%r(m)=r(m);tmp%n(l:m)=n(l:m)
  call alloc(tmp)
  allocate(x(maxval(r(l-1:m))),y(maxval(r(l-1:m))*maxval(n(l:m))),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif

  forall(i=1:r(l-1),j=1:n(l),k=1:r(l)+1)tmp%u(l)%p(i,j,k)=0.d0
  do p=l+1,m
   forall(i=1:r(p-1),j=1:n(p),k=1:r(p))tmp%u(p)%p(i,j,k)=arg%u(p)%p(i,j,k)
  end do 
  do p=l+1,m-1
   forall(           j=1:n(p),k=1:r(p))tmp%u(p)%p(r(p-1)+1,j,k)=0.d0
   forall(i=1:r(p-1),j=1:n(p)         )tmp%u(p)%p(i,j,r(p)+1)=0.d0
   forall(           j=1:n(p)         )tmp%u(p)%p(r(p-1)+1,j,r(p)+1)=0.d0
  end do
  
  select case(d)
   case(+1)
    forall(i=1:r(l-1),j=2:n(l),k=1:r(l))tmp%u(l)%p(i,j,k)=arg%u(l)%p(i,j-1,k)
    x(1:r(l))=arg%u(l)%p(1,n(l),:)
    tmp%u(l)%p(1,1,r(l)+1)=1.d0
    do p=l+1,m-1;tmp%u(p)%p(r(p-1)+1,1,r(p)+1)=1.d0;enddo
    tmp%u(m)%p(r(m-1)+1,1,r(m))=val
    do p=l+1,m
     call dgemv('t',r(p-1),n(p)*r(p),1.d0,arg%u(p)%p,r(p-1),x,1,0.d0,y,1)
     forall(j=1:n(p)-1,k=1:r(p))tmp%u(p)%p(r(p-1)+1,j+1,k)=y(j+(k-1)*n(p))
     forall(k=1:r(p))x(k)=y(k*n(p))
    end do
   case(-1)
    forall(i=1:r(l-1),j=1:n(l)-1,k=1:r(l))tmp%u(l)%p(i,j,k)=arg%u(l)%p(i,j+1,k)
    x(1:r(l))=arg%u(l)%p(1,1,:)
    tmp%u(l)%p(1,n(l),r(l)+1)=1.d0
    do p=l+1,m-1;tmp%u(p)%p(r(p-1)+1,n(p),r(p)+1)=1.d0;enddo
    tmp%u(m)%p(r(m-1)+1,n(m),r(m))=val
    do p=l+1,m
     call dgemv('t',r(p-1),n(p)*r(p),1.d0,arg%u(p)%p,r(p-1),x,1,0.d0,y,1)
     forall(j=2:n(p),k=1:r(p))tmp%u(p)%p(r(p-1)+1,j-1,k)=y(j+(k-1)*n(p))
     forall(k=1:r(p))x(k)=y(1+(k-1)*n(p))
    end do
   case default
    write(*,*)subnam,': illegal d: ',d
    stop
  end select
  arg=tmp;call dealloc(tmp)
  deallocate(x,y)
 end subroutine

 subroutine dtt_goup(arg)
  ![bitt] add leading bit by linear interpolation
  implicit none
  type(dtt),intent(inout),target :: arg
  character(len=*),parameter :: subnam='dtt_goup'
  type(dtt) :: tmp
  integer :: k,l,m,i,j,p,info
  integer,pointer :: r(:),n(:)
  double precision,allocatable :: x(:),y(:)
  tmp%l=arg%l; tmp%m=arg%m+1; l=arg%l; m=arg%m
  if(l.gt.m)return
  if(arg%r(l-1).ne.1)then;write(*,*)subnam,': left border rank should be 1';stop;endif
  if(.not.all(arg%n(l:m)==2))then;write(*,*)subnam,': not bitt on input';stop;endif
  r=>arg%r; n=>arg%n
  tmp%r(l-1)=1;tmp%r(l:m)=arg%r(l-1:m-1)+1;tmp%r(m+1)=arg%r(m); tmp%n(l:m+1)=2
  call alloc(tmp)

  allocate(x(maxval(r(l-1:m-1))),y(maxval(r(l-1:m-1))),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate: ',info;stop;endif
  
  tmp%u(l)%p=reshape( (/1.d0, 0.5d0, 0.d0, 0.5d0/), [1,2,2])
  
  tmp%u(l+1)%p(1,1,1:r(l))=arg%u(l)%p(1,1,:); tmp%u(l+1)%p(1,1,1+r(l))=0.d0
  tmp%u(l+1)%p(2,1,1:r(l))=arg%u(l)%p(1,2,:); tmp%u(l+1)%p(2,1,1+r(l))=0.d0
  tmp%u(l+1)%p(1,2,1:r(l))=arg%u(l)%p(1,2,:); tmp%u(l+1)%p(1,2,1+r(l))=0.d0
  tmp%u(l+1)%p(2,2,1:r(l))=0.d0;              tmp%u(l+1)%p(2,2,1+r(l))=1.d0
  x(1:r(l))=arg%u(l)%p(1,1,:)

  do k=l+1,m
   forall(i=1:r(k-1),p=1:2,j=1:r(k)) tmp%u(k+1)%p(i,p,j)=arg%u(k)%p(i,p,j)
   call dgemv('t',r(k-1),r(k),1.d0,arg%u(k)%p(1,2,1),r(k-1)*2,x,1,0.d0,y,1)
   tmp%u(k+1)%p(1+r(k-1),1,1:r(k))=y(1:r(k))
   tmp%u(k+1)%p(1+r(k-1),2,1:r(k))=0.d0
   if(k.lt.m)then
    forall(i=1:r(k-1),p=1:2) tmp%u(k+1)%p(i,p,r(k)+1)=0.d0
    tmp%u(k+1)%p(1+r(k-1),1,  r(k)+1)=0.d0
    tmp%u(k+1)%p(1+r(k-1),2,  r(k)+1)=1.d0
    call dgemv('t',r(k-1),r(k),1.d0,arg%u(k)%p,r(k-1)*2,x,1,0.d0,y,1)
    x(1:r(k))=y(1:r(k))
   end if
  end do
  arg=tmp
  call dealloc(tmp)
  deallocate(x,y)
 end subroutine
 
 subroutine mindex(ind,n,i)
  implicit none
  integer,intent(in) :: ind, n(:)
  integer,intent(out) :: i(:)
  integer :: j
  if(size(i).gt.size(n)+1)then;i=-999;return;endif
  i(1)=ind-1
  do j=2,size(i)
   i(j)=i(j-1)/n(j-1)
   i(j-1)=i(j-1)-i(j)*n(j-1)
  enddo
  i=i+1
 end subroutine
end module
