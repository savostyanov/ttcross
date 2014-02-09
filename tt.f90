module tt_lib
 use ptype_lib
 use default_lib
 implicit none
 integer,parameter :: tt_size=2048              ! Maximum tensor train size
 double precision,parameter :: errval = -999.d0

 type,public:: dtt            ! double precision tensor train
  integer :: l=1              ! index of the leftmost core
  integer :: m=0              ! index of the rightmost core
  integer :: n(tt_size)=0     ! mode sizes (storage for matrices)
  integer :: q(tt_size)=0     ! first mode sizes (for matrices)
  integer :: s(tt_size)=0     ! second mode sizes (for matrices)
  integer :: t=0              ! data type for future sparsity etc.
  integer :: r(0:tt_size)=0   ! TT ranks
  type(pointd3) :: u(tt_size) ! TT cores
 !contains
 ! final :: dtt_dealloc  ! finalization not yet supported by most compilers
 end type
 type,public:: ztt
  integer :: l=1
  integer :: m=0
  integer :: n(tt_size)=0
  integer :: q(tt_size)=0
  integer :: s(tt_size)=0
  integer :: t=0
  integer :: r(0:tt_size)=0
  type(pointz3) :: u(tt_size)
 !contains
 ! final :: ztt_dealloc  ! finalization not yet supported by most compilers
 end type

 interface tijk         ! particular element from tt-vector as function of index
  module procedure dtt_ijk,ztt_ijk
 end interface
 
 interface full         ! full array from tt-vector
  module procedure dtt_full, ztt_full
 end interface

 interface alloc        ! allocation and deallocation of the cores
  module procedure dtt_alloc,ztt_alloc
 end interface
 interface dealloc
  module procedure dtt_dealloc,ztt_dealloc
 end interface
 
 interface copy
  module procedure dtt_copy,ztt_copy
 end interface
 interface assignment (=)
  module procedure dtt_assign, ztt_assign, ztt_dtt_assign
 end interface

contains


 pure double precision function dtt_ijk(arg,ind) result (a)
  implicit none
  type(dtt),intent(in) :: arg
  integer,intent(in) :: ind(:)
  character(len=*),parameter :: subnam='dtt_ijk'
  integer :: info,i,l,m,n(tt_size),r(0:tt_size)
  double precision,pointer :: x(:,:),y(:,:),z(:,:)
  l=arg%l;m=arg%m;n=arg%n;r=arg%r
  if(any(ind(1:m-l+1)<=0).or.any(ind(1:m-l+1)>n(l:m)))then;a=-3.d0;return;endif
  if(r(l-1).ne.1 .or. r(m).ne.1)then;a=-4.d0;return;endif
  allocate(x(r(m-1),r(m)),stat=info)
  if(info.ne.0)then;a=-1.d0;return;endif
  x=arg%u(m)%p(:,ind(m-l+1),:)
  do i=m-1,l,-1
   allocate(y(r(i-1),r(i)),z(r(i-1),r(m)),stat=info)
   if(info.ne.0)then;a=-2.d0;return;endif
   y=arg%u(i)%p(:,ind(i-l+1),:)
   z=matmul(y,x)
   deallocate(x,y); x=>z; nullify(z)
  end do
  a=x(1,1)
  deallocate(x)
 end function
 pure double complex function ztt_ijk(arg,ind) result (a)
  implicit none
  type(ztt),intent(in) :: arg
  integer,intent(in) :: ind(:)
  character(len=*),parameter :: subnam='ztt_ijk'
  integer :: info,i,l,m,n(tt_size),r(0:tt_size)
  double complex,pointer :: x(:,:),y(:,:),z(:,:)
  double complex,parameter :: one=(1.d0,0.d0),zero=(0.d0,0.d0),im1=(0.d0,1.d0)
  l=arg%l;m=arg%m;n=arg%n;r=arg%r
  if(any(ind(1:m-l+1)<=0).or.any(ind(1:m-l+1)>n(l:m)))then;a=-3*one;return;endif
  if(r(l-1).ne.1 .or. r(m).ne.1)then;a=-4*one;return;endif
  allocate(x(r(m-1),r(m)),stat=info)
  if(info.ne.0)then;a=-one;return;endif
  x=arg%u(m)%p(:,ind(m-l+1),:)
  do i=m-1,l,-1
   allocate(y(r(i-1),r(i)),z(r(i-1),r(m)),stat=info)
   if(info.ne.0)then;a=-2*one;return;endif
   y=arg%u(i)%p(:,ind(i-l+1),:)
   z=matmul(y,x)
   deallocate(x,y); x=>z; nullify(z)
  end do
  a=x(1,1)
  deallocate(x)
 end function

 subroutine dtt_full(arg,a,alpha,beta,part,ind)
  ! A = [beta]*A + [alpha]*FULL(TT), TT = arg([l:m]), l,m=[part]
  ! A size r(l-1) *n(l)*...*n(m)* r(m)
  implicit none
  type(dtt),intent(in),target :: arg
  double precision,intent(inout) :: a(*)
  double precision,intent(in),optional :: alpha,beta
  integer,intent(in),optional :: part(2),ind(:)
  character(len=*),parameter :: subnam='dtt_full'
  double precision :: alp,bet
  type(pointd) :: p(0:1)
  integer,pointer :: r(:),n(:)
  integer :: l,m,na,nb,mem,rr,info,i,j,pp
  integer,allocatable :: ii(:),nn(:)
  double precision,allocatable :: q(:)

  alp=default(1.d0,alpha)
  bet=default(0.d0,beta)
  if(present(part))then;l=part(1);m=part(2);else;l=arg%l;m=arg%m;endif
  r=>arg%r; n=>arg%n

  if(l.gt.m)then
   write(*,*)subnam,': empty input';
   do i=1,r(l-1) * product(nn(l:m)) * r(m); a(i)=0.d0; enddo
   return
  end if

  allocate(ii(l:m),nn(l:m)); ii=0; nn(l:m)=arg%n(l:m)
  if(present(ind))then
   do i=l,m
    ii(i)=ind(i-l+1)
    if(ii(i).lt.0 .or. ii(i).gt.n(i))then;write(*,*)subnam,': invalid ind:',ind(1:m-l+1);stop;endif
    if(ii(i).ne.0)nn(i)=1
   end do
  end if

  na=r(l-1) * product(nn(l:m)) * r(m)
  mem=na; rr=1
  do i=l,m
   if(.not.(r(l-1)*product(nn(l:i))*r(i).gt.0))then;write(*,*)subnam,': oversized mem';stop;endif
   mem=max(mem, r(l-1)*product(nn(l:i))*r(i) )
   rr=max(rr,r(i-1)*r(i))
  end do
  allocate(p(0)%p(mem),p(1)%p(mem),q(rr),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory';stop;endif

  if(ii(l).eq.0)then
   call dcopy(r(l-1)*n(l)*r(l),arg%u(l)%p,1,p(0)%p,1)
  else
   do j=1,r(l);call dcopy(r(l-1),arg%u(l)%p(1,ii(l),j),1,p(0)%p(1+(j-1)*r(l-1)),1); enddo
  end if
  pp=0
  do i=l+1,m
   nb=r(l-1) * product(nn(l:i-1))
   if(ii(i).eq.0)then
    if(nb*n(i)*r(i).gt.mem)then;write(*,*)subnam,': nb-by-n-by-r > mem: ',nb,n(i),r(i),mem;stop;endif
    call dgemm('n','n',nb,n(i)*r(i),r(i-1),1.d0,p(pp)%p,nb,arg%u(i)%p,r(i-1),0.d0,p(1-pp)%p,nb)
   else
    if(nb*r(i).gt.mem)then;write(*,*)subnam,': nb-by-r > mem: ',nb,r(i),mem;stop;endif
    do j=1,r(i);call dcopy(r(i-1),arg%u(i)%p(1,ii(i),j),1,q(1+(j-1)*r(i-1)),1);enddo
    call dgemm('n','n',nb,r(i),r(i-1),1.d0,p(pp)%p,nb,q,r(i-1),0.d0,p(1-pp)%p,nb)
   end if
   pp=1-pp
  end do

  if(bet.eq.0.d0)then
   call dcopy(na,p(pp)%p,1,a,1)
   call dscal(na,alp,a,1)
  else
   call dscal(na,bet,a,1)
   call daxpy(na,alp,p(pp)%p,1,a,1)
  endif
  deallocate(p(0)%p,p(1)%p,q,ii,nn)
 end subroutine
 subroutine ztt_full(arg,a,alpha,beta,part,ind)
  ! A = [beta]*A + [alpha]*FULL(TT), TT = arg([l:m]), l,m=[part]
  ! A sizes r(l-1) *n(l)*...*n(m)* r(m)
  implicit none
  type(ztt),intent(in),target :: arg
  double complex,intent(inout) :: a(*)
  double complex,intent(in),optional :: alpha,beta
  integer,intent(in),optional :: part(2),ind(:)
  character(len=*),parameter :: subnam='ztt_full'
  double complex,parameter :: one=(1.d0,0.d0),zero=(0.d0,0.d0)
  double complex :: alp,bet
  type(pointz) :: p(0:1)
  integer,pointer :: r(:),n(:)
  integer :: l,m,na,nb,mem,rr,info,i,j,pp
  integer,allocatable :: ii(:),nn(:)
  double complex,allocatable :: q(:)

  alp=default(one,alpha)
  bet=default(zero,beta)
  if(present(part))then;l=part(1);m=part(2);else;l=arg%l;m=arg%m;endif
  r=>arg%r; n=>arg%n

  allocate(ii(l:m),nn(l:m)); ii=0; nn(l:m)=arg%n(l:m)
  if(present(ind))then
   do i=l,m
    ii(i)=ind(i-l+1)
    if(ii(i).lt.0 .or. ii(i).gt.n(i))then;write(*,*)subnam,': invalid ind:',ind(1:m-l+1);stop;endif
    if(ii(i).ne.0)nn(i)=1
   end do
  end if

  na=r(l-1) * product(nn(l:m)) * r(m)
  mem=na; rr=1
  do i=l,m
   if(.not.(r(l-1)*product(nn(l:i))*r(i).gt.0))then;write(*,*)subnam,': oversized mem';stop;endif
   mem=max(mem, r(l-1)*product(nn(l:i))*r(i) )
   rr=max(rr,r(i-1)*r(i))
  end do
  allocate(p(0)%p(mem),p(1)%p(mem),q(rr),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory';stop;endif

  if(ii(l).eq.0)then
   call zcopy(r(l-1)*n(l)*r(l),arg%u(l)%p,1,p(0)%p,1)
  else
   do j=1,r(l);call zcopy(r(l-1),arg%u(l)%p(1,ii(l),j),1,p(0)%p(1+(j-1)*r(l-1)),1); enddo
  end if
  pp=0
  do i=l+1,m
   nb=r(l-1) * product(nn(l:i-1))
   if(ii(i).eq.0)then
    if(nb*n(i)*r(i).gt.mem)then;write(*,*)subnam,': nb-by-n-by-r > mem: ',nb,n(i),r(i),mem;stop;endif
    call zgemm('n','n',nb,n(i)*r(i),r(i-1),one,p(pp)%p,nb,arg%u(i)%p,r(i-1),zero,p(1-pp)%p,nb)
   else
    if(nb*r(i).gt.mem)then;write(*,*)subnam,': nb-by-r > mem: ',nb,r(i),mem;stop;endif
    do j=1,r(i);call zcopy(r(i-1),arg%u(i)%p(1,ii(i),j),1,q(1+(j-1)*r(i-1)),1);enddo
    call dgemm('n','n',nb,r(i),r(i-1),one,p(pp)%p,nb,q,r(i-1),zero,p(1-pp)%p,nb)
   end if
   pp=1-pp
  end do

  if(bet.eq.zero)then
   call zcopy(na,p(pp)%p,1,a,1)
   call zscal(na,alp,a,1)
  else
   call zscal(na,bet,a,1)
   call zaxpy(na,alp,p(pp)%p,1,a,1)
  endif
  deallocate(p(0)%p,p(1)%p,q,ii,nn)
 end subroutine

 integer function dtt_mem(arg) result (sz)
  implicit none
  type(dtt),intent(in) :: arg
  integer :: i
  sz=0
  do i=arg%l,arg%m
   sz=sz+arg%r(i-1)*arg%n(i)*arg%r(i)
  end do
 end function
 integer function ztt_mem(arg) result (sz)
  implicit none
  type(ztt),intent(in) :: arg
  integer :: i
  sz=0
  do i=arg%l,arg%m
   sz=sz+arg%r(i-1)*arg%n(i)*arg%r(i)
  end do
 end function
 double precision function dtt_mb(arg) result (sz)
  implicit none
  type(dtt),intent(in) :: arg
  integer :: i
  sz=0.d0
  do i=arg%l,arg%m
   sz=sz+arg%r(i-1)*arg%n(i)*arg%r(i)
  end do
  sz=sz*8.d0/(2**20)
 end function
 double precision function ztt_mb(arg) result (sz)
  implicit none
  type(ztt),intent(in) :: arg
  integer :: i
  sz=0.d0
  do i=arg%l,arg%m
   sz=sz+arg%r(i-1)*arg%n(i)*arg%r(i)
  end do
  sz=sz*16.d0/(2**20)
 end function
 
 subroutine dtt_alloc(arg)
  implicit none
  type(dtt),intent(inout) :: arg
  character(len=*),parameter :: subnam='dtt_alloc'
  integer :: i,info
  if(arg%m.lt.arg%l)return
  if(arg%l.le.0)then;write(*,*)subnam,': %l should be > 0';stop;endif
  if(arg%m.gt.tt_size)then;write(*,*)subnam,': %m exceeds tt_size, change parameter and recompile!';stop;endif
  do i=arg%l,arg%m
   if(associated(arg%u(i)%p))deallocate(arg%u(i)%p)
   allocate(arg%u(i)%p(arg%r(i-1),arg%n(i),arg%r(i)), stat=info)
   if(info.ne.0)then;write(*,*)'TT allocate fail: no memory';stop;endif
  end do
 end subroutine
 subroutine ztt_alloc(arg)
  implicit none
  type(ztt),intent(inout) :: arg
  character(len=*),parameter :: subnam='ztt_alloc'
  integer :: i,info
  if(arg%m.lt.arg%l)return
  if(arg%l.le.0)then;write(*,*)subnam,': %l should be > 0';stop;endif
  if(arg%m.gt.tt_size)then;write(*,*)subnam,': %m exceeds tt_size, change parameter and recompile!';stop;endif
  do i=arg%l,arg%m
   if(associated(arg%u(i)%p))deallocate(arg%u(i)%p)
   allocate(arg%u(i)%p(arg%r(i-1),arg%n(i),arg%r(i)), stat=info)
   if(info.ne.0)then;write(*,*)'TT allocate fail: no memory';stop;endif
  end do
 end subroutine

 subroutine dtt_dealloc(arg)
  implicit none
  type(dtt),intent(inout) :: arg
  character(len=*),parameter :: subnam='dtt_dealloc'
  integer :: i
  do i=1,tt_size
   if(associated(arg%u(i)%p))deallocate(arg%u(i)%p)
  end do
 end subroutine
 subroutine ztt_dealloc(arg)
  implicit none
  type(ztt),intent(inout) :: arg
  character(len=*),parameter :: subnam='ztt_dealloc'
  integer :: i
  do i=1,tt_size
   if(associated(arg%u(i)%p))deallocate(arg%u(i)%p)
  end do
 end subroutine

 
 subroutine dtt_assign(b,a)
  implicit none
  type(dtt),intent(inout) :: b
  type(dtt),intent(in) :: a
  integer :: k,l,m
  l=a%l;m=a%m
  b%l=l; b%m=m; b%n(l:m)=a%n(l:m); b%r(l-1:m)=a%r(l-1:m); call alloc(b)
  do k=l,m; call dcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(k)%p,1); end do
 end subroutine
 subroutine ztt_assign(b,a)
  implicit none
  type(ztt),intent(inout) :: b
  type(ztt),intent(in) :: a
  integer :: k,l,m
  l=a%l;m=a%m
  b%l=l; b%m=m; b%n(l:m)=a%n(l:m); b%r(l-1:m)=a%r(l-1:m); call alloc(b)
  do k=l,m; call zcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(k)%p,1); end do
 end subroutine
 subroutine ztt_dtt_assign(z,d)
  implicit none
  type(ztt),intent(inout) :: z
  type(dtt),intent(in) :: d
  integer :: l,m,i,j,k,p,r(0:tt_size),n(tt_size)
  l=d%l; m=d%m; r=d%r; n=d%n
  if(l.gt.m)then;write(*,*)'ztt_dtt_assign: l,m: ',l,m;return;endif
  z%l=l;z%m=m;z%n(l:m)=n(l:m);z%r(l-1:m)=r(l-1:m)
  call alloc(z)
  do p=l,m
   forall(i=1:r(p-1),j=1:n(p),k=1:r(p))z%u(p)%p(i,j,k)=dcmplx(d%u(p)%p(i,j,k),0.d0)
  end do
 end subroutine

 subroutine dtt_copy(a,b,low)
  implicit none
  type(dtt),intent(in) :: a
  type(dtt),intent(inout) :: b
  integer,intent(in),optional :: low
  integer :: k,l,m,ll,mm
  l=a%l;m=a%m; ll=default(b%l,low); mm=ll-l+m
  b%l=ll; b%m=mm; b%n(ll:mm)=a%n(l:m); b%r(ll-1:mm)=a%r(l-1:m);
  if(.not.all(a%n(l:m)>0))return;if(.not.all(a%r(l-1:m)>0))return;call alloc(b)
  do k=l,m; call dcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(ll-l+k)%p,1); end do
 end subroutine
 subroutine ztt_copy(a,b,low)
  implicit none
  type(ztt),intent(in) :: a
  type(ztt),intent(inout) :: b
  integer,intent(in),optional :: low
  integer :: k,l,m,ll,mm
  l=a%l;m=a%m; ll=default(b%l,low); mm=ll-l+m
  b%l=ll; b%m=mm; b%n(ll:mm)=a%n(l:m); b%r(ll-1:mm)=a%r(l-1:m);
  if(.not.all(a%n(l:m)>0))return;if(.not.all(a%r(l-1:m)>0))return;call alloc(b)
  do k=l,m; call zcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(ll-l+k)%p,1); end do
 end subroutine

end module
