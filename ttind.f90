module ttind_lib
 use tt_lib
 type,public:: ttind
  integer :: p(tt_size)=0
  integer :: n(tt_size)=0
  integer :: m=0
 end type
 interface say
  module procedure ttind_say
 end interface
 interface rnd
  module procedure ttind_rnd
 end interface
 interface reverse
  module procedure ttind_rev
 end interface


 interface assignment (=)
  module procedure ttind_assign
 end interface
 interface operator (.eq.)
  module procedure ttind_eq
 end interface 
 interface operator (.lt.)
  module procedure ttind_lt
 end interface 
 interface operator (.le.)
  module procedure ttind_le
 end interface 
 interface operator (.gt.)
  module procedure ttind_gt
 end interface 
 interface operator (.ge.)
  module procedure ttind_ge
 end interface 

 interface find
  module procedure find_ttind
 end interface
 interface push
  module procedure push_ttind
 end interface
 interface dble
  module procedure dble_ttind
 end interface
 interface int
  module procedure int_ttind
 end interface

contains
! SAY 
 subroutine ttind_say(arg)
  implicit none
  type(ttind),intent(in) :: arg
  write(*,'(a,i2)') 'ttind, m=', arg%m
  write(*,'(a,1x,128i3)') 'n: ',arg%n(1:arg%m)
  write(*,'(a,128i3)') 'p: ',arg%p(1:arg%m)
 end subroutine

! RND 
 subroutine ttind_rnd(arg)
  use rnd_lib
  implicit none
  type(ttind),intent(inout) :: arg
  character(len=*),parameter :: subnam='ttind_rnd'
  integer :: l,m,k
  l=1;m=arg%m
  if(m.lt.l)return
  if(.not.all(arg%n(l:m)>0))then;write(*,*)subnam,': n: ',arg%n(l:m);stop;endif
  do k=l,m
   arg%p(k)=irnd(arg%n(k))
  end do
 end subroutine

! REVerse
 subroutine ttind_rev(arg)
  implicit none
  type(ttind),intent(inout) :: arg
  integer :: m,k,q
  m=arg%m
  if(m.le.0)return
  do k=1,m/2
   q=m-k+1
   arg%n(k)=arg%n(k)+arg%n(q); arg%n(q)=arg%n(k)-arg%n(q); arg%n(k)=arg%n(k)-arg%n(q)
   arg%p(k)=arg%p(k)+arg%p(q); arg%p(q)=arg%p(k)-arg%p(q); arg%p(k)=arg%p(k)-arg%p(q)
  end do
 end subroutine 

! TTINDEX convertion tools
 pure type(ttind) function ttindex(i,n) result (ind)
  implicit none
  integer,intent(in) :: i,n(:)
  integer :: m
  ind%p=0; m=1; ind%p(m)=i-1;
  do while(m.lt.tt_size .and. m.le.size(n) .and. ind%p(m).ge.n(m))
   m=m+1
   ind%p(m)=ind%p(m-1)/n(m-1)
   ind%p(m-1)=ind%p(m-1)-ind%p(m)*n(m-1)
  enddo
  ind%m=m
  ind%p(1:m)=ind%p(1:m)+1
  m=min(size(n),tt_size)
  ind%n(1:m)=n(1:m)
 end function

 
 pure double precision function dble_ttind(ind) result (x)
  implicit none
  type(ttind),intent(in) :: ind
  integer :: j
  x=dble(ind%p(ind%m))-1
  do j=ind%m-1,1,-1
   x=x*ind%n(j)
   x=x+ind%p(j)-1
  end do
  x=x+1
 end function
 pure integer function int_ttind(ind) result (i)
  implicit none
  type(ttind),intent(in) :: ind
  integer :: j
  i=ind%p(ind%m)-1
  do j=ind%m-1,1,-1
   i=i*ind%n(j)
   i=i+ind%p(j)-1
  end do
  i=i+1
 end function

! FIND
 integer function find_ttind(n,x,y) result (pos)
  ! for sorted vector x(1) <= x(2) <= ... <= x(n) and value y find pos, s.t. x(pos) <= y < x(pos+1)
  implicit none
  integer,intent(in) :: n
  type(ttind),intent(in) :: x(n),y
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

! PUSH
 subroutine push_ttind(n,x)
  implicit none
  integer,intent(in) :: n
  type(ttind),intent(inout) :: x(n+1)
  integer :: i
  if(n.le.0)return
  do i=n,1,-1
   x(i+1)=x(i)
  end do
  x(1)%p=0;x(1)%m=0;x(1)%n=0
 end subroutine

! ASSIGN
 subroutine ttind_assign(b,a)
  implicit none
  type(ttind),intent(inout) :: b
  type(ttind),intent(in) :: a
  b%p=a%p; b%n=a%n; b%m=a%m
 end subroutine

! EQ LE GE LT GT
 pure logical function ttind_eq(a,b) result (key)
  type(ttind),intent(in) :: a,b
  key=(a%m.eq.b%m); if(.not.key)return
  key=all(a%p(1:a%m)==b%p(1:b%m))
 end function
 pure logical function ttind_lt(a,b) result (key)
  type(ttind),intent(in) :: a,b
  integer :: i
  if(a%m.gt.b%m)then;if(any(a%p(b%m+1:a%m)>1))then;key=.false.;return;endif;endif
  if(b%m.gt.a%m)then;if(any(b%p(a%m+1:b%m)>1))then;key=.true. ;return;endif;endif
  key=.false.; i=min(a%m,b%m)
  do while(i.gt.0 .and. a%p(i).eq.b%p(i)); i=i-1; end do
  if(i.gt.0)key=(a%p(i).lt.b%p(i))
 end function
 pure logical function ttind_le(a,b) result (key)
  type(ttind),intent(in) :: a,b
  integer :: i
  if(a%m.gt.b%m)then;if(any(a%p(b%m+1:a%m)>1))then;key=.false.;return;endif;endif
  if(b%m.gt.a%m)then;if(any(b%p(a%m+1:b%m)>1))then;key=.true. ;return;endif;endif
  key=.true.; i=min(a%m,b%m)
  do while(i.gt.0 .and. a%p(i).eq.b%p(i)); i=i-1; end do
  if(i.gt.0)key=(a%p(i).lt.b%p(i))
 end function
 pure logical function ttind_gt(a,b) result (key)
  type(ttind),intent(in) :: a,b
  integer :: i
  if(a%m.gt.b%m)then;if(any(a%p(b%m+1:a%m)>1))then;key=.true. ;return;endif;endif
  if(b%m.gt.a%m)then;if(any(b%p(a%m+1:b%m)>1))then;key=.false.;return;endif;endif
  key=.false.; i=min(a%m,b%m)
  do while(i.gt.0 .and. a%p(i).eq.b%p(i)); i=i-1; end do
  if(i.gt.0)key=(a%p(i).gt.b%p(i))
 end function
 pure logical function ttind_ge(a,b) result (key)
  type(ttind),intent(in) :: a,b
  integer :: i
  if(a%m.gt.b%m)then;if(any(a%p(b%m+1:a%m)>1))then;key=.true. ;return;endif;endif
  if(b%m.gt.a%m)then;if(any(b%p(a%m+1:b%m)>1))then;key=.false.;return;endif;endif
  key=.true.; i=min(a%m,b%m)
  do while(i.gt.0 .and. a%p(i).eq.b%p(i)); i=i-1; end do
  if(i.gt.0)key=(a%p(i).gt.b%p(i))
 end function


end module
