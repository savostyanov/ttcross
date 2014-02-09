module lr_lib
 use default_lib
 implicit none

 contains

 subroutine d2_lrg(m,n,r,a,u,v,ind,jnd)
 implicit none
  integer,intent(in) :: m,n,r
  double precision,intent(in) :: a(m,n)
  double precision,intent(out) :: u(m,r),v(r,n)
  integer,intent(out) :: ind(r),jnd(r)
  character(len=*),parameter :: subnam='d2_lrg'
  double precision,allocatable :: e(:,:)
  integer :: info,p,ij,i,j
  integer,external :: idamax
  allocate(e(m,n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': alloc fail';stop;endif
  call dcopy(m*n,a,1,e,1)
  do p=1,r
   ij=idamax(m*n,e,1)
   j=(ij-1)/m+1; i=ij-(j-1)*m
   ind(p)=i; jnd(p)=j
   call dcopy(m,e(1,j),1,u(1,p),1)
   call dcopy(n,e(i,1),m,v(p,1),r)
   call dscal(n,1.d0/e(i,j),v(p,1),r)
   call dger(m,n,-1.d0,u(1,p),1,v(p,1),r,e,m)
  end do
  deallocate(e)
 end subroutine
 
 subroutine d2_lug(r,a,g)
 implicit none
  integer,intent(in) :: r
  double precision,intent(in) :: a(r,r)
  double precision,intent(out) :: g(r*r)
  character(len=*),parameter :: subnam='d2_lug'
  double precision,allocatable :: e(:,:)
  integer :: info,p,pos,i
  double precision :: aa
  allocate(e(r,r),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': alloc fail';stop;endif
  call dcopy(r*r,a,1,e,1)
  pos=0
  do p=1,r
   aa=1.d0/e(p,p)
   g(pos+1)=aa; pos=pos+1
   call dscal(r-p,aa,e(p,p+1),r)
   if(p.gt.1)then
    forall(i=1:p-1)g(pos+i)=e(i,p); pos=pos+p-1
    forall(i=1:p-1)g(pos+i)=e(p,i); pos=pos+p-1
   end if
   call dger(r-p,r-p,-1.d0,e(p+1,p),1,e(p,p+1),r,e(p+1,p+1),r)
  end do
  deallocate(e)
 end subroutine
 
 subroutine d2_luac(m,r,g,col,u,from)
 implicit none
  integer,intent(in) :: m,r
  double precision,intent(in) :: g(r*r),col(m,r)
  double precision,intent(inout) :: u(m,r)
  integer,intent(in),optional :: from
  character(len=*),parameter :: subnam='d2_luac'
  integer :: info,p,pos,i,p0
  double precision :: aa
  p0=default(1,from)
  pos=(p0-1)**2
  do p=p0,r
   aa=g(pos+1); pos=pos+1
   call dcopy(m,col(1,p),1,u(1,p),1)
   if(p.gt.1)then
    call dgemv('n',m,p-1,-1.d0,u,m,g(pos+1),1,1.d0,u(1,p),1); pos=pos+p-1
    pos=pos+p-1
   end if
  end do
 end subroutine
 subroutine d2_luar(r,n,g,row,v,from)
 implicit none
  integer,intent(in) :: r,n
  double precision,intent(in) :: g(r*r),row(r,n)
  double precision,intent(inout) :: v(r,n)
  integer,intent(in),optional :: from
  character(len=*),parameter :: subnam='d2_luar'
  integer :: info,p,pos,i,p0
  double precision :: aa
  p0=default(1,from)
  pos=(p0-1)**2
  do p=p0,r
   aa=g(pos+1); pos=pos+1
   call dcopy(n,row(p,1),r,v(p,1),r)
   if(p.gt.1)then
    pos=pos+p-1
    call dgemv('t',p-1,n,-1.d0,v,r,g(pos+1),1,1.d0,v(p,1),r); pos=pos+p-1
   end if
   call dscal(n,aa,v(p,1),r)
  end do
 end subroutine
 subroutine d2_lua(m,n,r,g,col,row,u,v,from)
 implicit none
  integer,intent(in) :: m,n,r
  double precision,intent(in) :: g(r*r),col(m,r),row(r,n)
  double precision,intent(inout) :: u(m,r),v(r,n)
  integer,intent(in),optional :: from
  character(len=*),parameter :: subnam='d2_lua'
  call d2_luac(m,r,g,col,u,from)
  call d2_luar(r,n,g,row,v,from)
 end subroutine
 
 subroutine d2_cgr(m,n,r,acol,arow,aa,a,alpha,beta)
 implicit none
  integer,intent(in) :: m,n,r
  double precision,intent(in) :: acol(m,r),arow(r,n),aa(r,r)
  double precision,intent(inout) :: a(m,n)
  double precision,intent(in),optional :: alpha,beta
  character(len=*),parameter :: subnam='d2_cgr'
  integer :: info
  double precision :: alp,bet
  double precision,allocatable :: u(:,:),v(:,:),g(:)
  alp=default(1.d0,alpha) 
  bet=default(0.d0,beta)
  allocate(u(m,r),v(r,n),g(r*r),stat=info)
  call d2_lug(r,aa,g)
  call d2_lua(m,n,r,g,acol,arow,u,v)
  call dgemm('n','n',m,n,r,alp,u,m,v,r,bet,a,m)
  deallocate(u,v,g)
 end subroutine
end module
