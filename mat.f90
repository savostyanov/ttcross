module mat_lib
 use nan_lib
 implicit none

 interface matinv
  module procedure matinv_d,matinv_z
 end interface
 private :: matinv_d,matinv_z,matinv_svd_d,matinv_svd_z

 interface svd
  ! A = USV
  module procedure d_svd, z_svd
 end interface
 interface eye
  module procedure eye_d2,eye_z2
 end interface
 interface laplace
  module procedure laplace_d2,laplace_z2
 end interface

contains
!MATINV 
 subroutine matinv_d(a,ainv,alg,tol)
  implicit none
  double precision,intent(inout)  :: a(:,:)
  double precision,intent(out),optional :: ainv(size(a,2),size(a,1))
  character(len=1),intent(in),optional :: alg
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_d'
  if(.not.present(alg) .or. size(a,1).ne.size(a,2))then
   call matinv_svd_d(a,ainv,tol)
  else
   select case(alg)
    case('s','S')
     call matinv_svd_d(a,ainv,tol)
    case('t','T')
     call matinv_lu_d(a,ainv)
    case default
     write(*,'(3a)')subnam,': unknown alg: ',alg
     stop
   end select
  end if
  return
 end subroutine 
 subroutine matinv_z(a,ainv,alg,tol)
  implicit none
  double complex,intent(inout)  :: a(:,:)
  double complex,intent(out),optional :: ainv(size(a,2),size(a,1))
  character(len=1),intent(in),optional :: alg
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_z'
  if(.not.present(alg) .or. size(a,1).ne.size(a,2))then
   call matinv_svd_z(a,ainv,tol)
  else
   select case(alg)
    case('s','S')
     call matinv_svd_z(a,ainv,tol)
    case('t','T')
     call matinv_lu_z(a,ainv)
    case default
     write(*,'(3a)')subnam,': unknown alg: ',alg
     stop
   end select
  end if
  return
 end subroutine 

 subroutine matinv_svd_d(a,ainv,tol)
  implicit none
  double precision,intent(inout)  :: a(:,:)
  double precision,intent(out),optional :: ainv(size(a,2),size(a,1))
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_svd_d'
  double precision,allocatable :: b(:,:),u(:,:),v(:,:),s(:),work(:)
  integer :: r,lwork,info,i,ntrunc,m,n
  double precision :: si,s1,small
  character(len=180) :: str
  
  m=size(a,1); n=size(a,2)
  if(present(tol))then;small=tol;else;small=1.d-14;endif
  write(str,'(a,a,i8,1x,i8)')subnam,': m,n:',m,n 
  !call plog(2,str)
  r=min(m,n); lwork=64*max(m,n)
  allocate(b(m,n),u(m,r),v(r,n),s(r),work(lwork))
  b=a
  call dgesvd('s','s',m,n,b,m,s,u,m,v,r, work,lwork,info)
  if(info.ne.0)then
   write(*,'(a,a,i10)')subnam,': dgesvd info: ',info
   stop
  end if
  
  s1=s(1)
  !write(str,'(a,80(e9.3,1x))')'sv:',(fmem(s+i)/s1,i=1,r)
  !call plog(1,str)

  ntrunc=0
  do i=1,r
   si=s(i)
    if(si.lt.small*s1)then
     s(i)=0.d0; ntrunc=ntrunc+1
    else
     s(i)=1.d0/s(i)
    end if
  end do
  
  forall(i=1:r)u(:,i)=u(:,i)*s(i)
  if(present(ainv))then
   call dgemm('t','t',n,m,r,1.d0,v,r,u,m,0.d0,ainv,n)
  else
   call dgemm('t','t',n,m,r,1.d0,v,r,u,m,0.d0,a,n)
  end if 
  if(ntrunc.gt.0)then
   write(*,'(a,a,i5,a,i5)')subnam,': truncate: ',ntrunc,' of: ',r
   !call plog(1,str)
  end if

  deallocate(b,u,v,s,work)
  return
 end subroutine
 subroutine matinv_svd_z(a,ainv,tol)
  implicit none
  double complex,intent(in)  :: a(:,:)
  double complex,intent(out) :: ainv(size(a,2),size(a,1))
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_svd_z'
  double complex,allocatable   :: b(:,:),u(:,:),v(:,:),work(:)
  double precision,allocatable :: s(:),rwork(:)
  integer :: r,lwork,lrwork,info,i,ntrunc,m,n
  double precision :: small=1.d-14
  double precision :: si,s1
  character(len=180) :: str
  
  m=size(a,1); n=size(a,2)
  if(present(tol))small=tol
  write(str,'(a,a,i8,1x,i8)')subnam,': m,n:',m,n 
  !call plog(2,str)
  r=min(m,n); lwork=64*max(m,n); lrwork=5*min(m,n)
  allocate(b(m,n),u(m,r),v(n,r),s(r),work(lwork),rwork(lrwork))
  b=a
  call zgesvd('s','s',m,n,b,m,s,u,m,v,n,work,lwork,rwork,info)
  if(info.ne.0)then
   write(*,'(a,a,i10)')subnam,': zgesvd info: ',info
   stop
  end if
  
  s1=s(1)
  !write(str,'(a,80(e9.3,1x))')'sv:',(fmem(s+i)/s1,i=1,r)
  !call plog(1,str)

  ntrunc=0
  do i=1,r
   si=s(i)
    if(si.lt.small*s1)then
     s(i)=0.d0; ntrunc=ntrunc+1
    else
     s(i)=1.d0/s(i)
    end if
  end do
  
  forall(i=1:r)u(:,i)=u(:,i)*s(i)
  call zgemm('c','c',n,m,r,(1.d0,0.d0),v,r,u,m,(0.d0,0.d0),ainv,n)
  if(ntrunc.gt.0)then
   write(str,'(a,a,i5,a,i5)')subnam,': truncate: ',ntrunc,' of: ',r
   !call plog(1,str)
  end if

  deallocate(b,u,v,s,work,rwork)
  return
 end subroutine

 subroutine matinv_lu_d(a,ainv)
  implicit none
  double precision,intent(inout),target  :: a(:,:)
  double precision,intent(out),optional,target :: ainv(size(a,2),size(a,1))
  character(len=*),parameter   :: subnam='matinv_lu_d'
  double precision,allocatable :: work(:)
  double precision,dimension(:,:),pointer :: aa
  integer :: lwork,m,n,info
  integer,allocatable :: piv(:)
  m=size(a,1); n=size(a,2)
  if(m.ne.n)then
   write(*,*)subnam,': matrix not square: ',m,n
   stop
  end if 
  if(present(ainv))then
   call dcopy(m*n,a,1,ainv,1)
   aa=>ainv
  else
   aa=>a
  end if
 
  lwork=64*n
  allocate(work(lwork),piv(n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call dgetrf(n,n,aa,n,piv,info)
  if(info.ne.0)then; write(*,*)subnam,': dgertf: info: ',info; stop; end if
  call dgetri(n,aa,n,piv,work,lwork,info)
  if(info.ne.0)then; write(*,*)subnam,': dgerti: info: ',info; stop; end if

  nullify(aa)
  deallocate(work,piv)
  return
 end subroutine
 subroutine matinv_lu_z(a,ainv)
  implicit none
  double complex,intent(inout),target  :: a(:,:)
  double complex,intent(out),optional,target :: ainv(size(a,2),size(a,1))
  character(len=*),parameter   :: subnam='matinv_lu_z'
  double complex,allocatable :: work(:)
  double complex,dimension(:,:),pointer :: aa
  integer :: lwork,m,n,info
  integer,allocatable :: piv(:)
  m=size(a,1); n=size(a,2)
  if(m.ne.n)then
   write(*,*)subnam,': matrix not square: ',m,n
   stop
  end if 
  if(present(ainv))then
   call zcopy(m*n,a,1,ainv,1)
   aa=>ainv
  else
   aa=>a
  end if

  lwork=64*n
  allocate(work(lwork),piv(n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif

  call zgetrf(n,n,aa,n,piv,info)
  if(info.ne.0)then; write(*,*)subnam,': dgertf: info: ',info; stop; end if
  call zgetri(n,aa,n,piv,work,lwork,info)
  if(info.ne.0)then; write(*,*)subnam,': dgerti: info: ',info; stop; end if

  deallocate(work,piv)
  return
 end subroutine

!EYE
 subroutine eye_d2(a)
  implicit none
  double precision,intent(inout) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  a=0.d0
  forall(i=1:min(m,n))a(i,i)=1.d0
  return
 end subroutine 
 subroutine eye_z2(a)
  implicit none
  double complex,intent(inout) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  a=(0.d0,0.d0)
  forall(i=1:min(m,n))a(i,i)=(1.d0,0.d0)
  return
 end subroutine 
 subroutine d2eye(a,n)
  implicit none
  double precision,intent(out) :: a(n,n)
  integer,intent(in) :: n
  integer :: i
  call dscal(n*n,0.d0,a,1)
  forall(i=1:n)a(i,i)=1.d0
 end subroutine 
 subroutine z2eye(n,a)
  implicit none
  integer,intent(in) :: n
  double complex,intent(out) :: a(n,n)
  integer :: i
  call zdscal(n*n,0.d0,a,1)
  forall(i=1:n)a(i,i)=(1.d0,0.d0)
 end subroutine 

!LAPLACE
 subroutine laplace_d2(a)
  implicit none
  double precision,intent(out) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  call dscal(m*n,0.d0,a,1)
  forall(i=1:min(m,n))a(i,i)=2.d0
  forall(i=1:min(m-1,n))a(i+1,i)=-1.d0
  forall(i=1:min(m,n-1))a(i,i+1)=-1.d0
  return
 end subroutine 
 subroutine laplace_z2(a)
  implicit none
  double complex,intent(out) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  call zdscal(m*n,0.d0,a,1)
  forall(i=1:min(m,n))a(i,i)=(2.d0,0.d0)
  forall(i=1:min(m-1,n))a(i+1,i)=-(1.d0,0.d0)
  forall(i=1:min(m,n-1))a(i,i+1)=-(1.d0,0.d0)
  return
 end subroutine 
 

!SUBMAT
 subroutine d2submat(m,n,a,lda,b)
  implicit none
  integer,intent(in) :: m,n,lda
  double precision,intent(in) :: a(lda,n)
  double precision,intent(out) :: b(m,n)
  integer :: i,j
  if(m.gt.lda)then;write(*,*)'d2submat: lda,m: ',lda,m;stop;endif
  if(m.eq.lda)then;call dcopy(m*n,a,1,b,1);return;endif
  forall(i=1:m,j=1:n)b(i,j)=a(i,j)
 end subroutine
 subroutine z2submat(m,n,a,lda,b)
  implicit none
  integer,intent(in) :: m,n,lda
  double complex,intent(in) :: a(lda,n)
  double complex,intent(out) :: b(m,n)
  integer :: i,j
  if(m.gt.lda)then;write(*,*)'d2submat: lda,m: ',lda,m;stop;endif
  if(m.eq.lda)then;call zcopy(m*n,a,1,b,1);return;endif
  forall(i=1:m,j=1:n)b(i,j)=a(i,j)
 end subroutine
 subroutine d2subset(m,n,r,ind,jnd,a,b)
  implicit none
  integer,intent(in) :: m,n,r
  integer,intent(in) :: ind(r),jnd(r)
  double precision,intent(in) :: a(m,n)
  double precision,intent(out) :: b(r,r)
  integer :: i,j
  forall(i=1:r,j=1:r)b(i,j)=a(ind(i),jnd(j))
 end subroutine
 subroutine z2subset(m,n,r,ind,jnd,a,b)
  implicit none
  integer,intent(in) :: m,n,r
  integer,intent(in) :: ind(r),jnd(r)
  double complex,intent(in) :: a(m,n)
  double complex,intent(out) :: b(r,r)
  integer :: i,j
  forall(i=1:r,j=1:r)b(i,j)=a(ind(i),jnd(j))
 end subroutine

! SVD
 subroutine d_svd(a,u,s,v,tol,rmax,err,info)
  implicit none
  double precision,intent(in) :: a(:,:)
  double precision,pointer :: u(:,:),s(:),v(:,:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='d_svd'
  double precision,allocatable :: b(:,:),work(:),ss(:),uu(:,:),vv(:,:)
  integer,allocatable :: iwork(:)
  integer :: m,n,mn,mx,lwork,ierr,r,i,j
  double precision,external :: dnrm2
  m=size(a,1); n=size(a,2); mn=min(m,n); mx=max(m,n); lwork=-1
  allocate(work(1), b(m,n),uu(m,mn),vv(mn,n),ss(mn))
  ! call dgesdd('s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,iwork,ierr)
  call dgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,ierr)
  lwork=int(work(1))+1
  deallocate(work)
  allocate(work(lwork),iwork(8*mn), stat=ierr)
  if(ierr.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call dcopy(m*n,a,1,b,1)
  call dgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,ierr)
  !call dgesdd('s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,iwork,ierr)
  if(present(info))info=ierr
  deallocate(b,work,iwork)
  if(ierr.ne.0)then
   write(*,*)subnam,': dgesvd info: ',ierr
   if(ierr.lt.0)stop
   if(nan(a))then
    write(*,*) subnam,': NaNs detected in the input array'
    stop
   else
    write(*,*) subnam,': min/max element of input: ',minval(a),maxval(a)
   end if 
   u=>null(); v=>null(); s=>null()
  else 
   r=chop(ss,tol,rmax,err)
   if(present(err))err=err/dnrm2(mn,ss,1)
   allocate(u(m,r),s(r),v(r,n))
   call dcopy(m*r,uu,1,u,1)
   call dcopy(r,ss,1,s,1)
   call d2submat(r,n,vv,mn,v)
  end if 
  deallocate(uu,vv,ss)
 end subroutine 
 subroutine z_svd(a,u,s,v,tol,rmax,err,info)
  implicit none
  double complex,intent(in) :: a(:,:)
  double complex,pointer :: u(:,:),v(:,:)
  double precision,pointer :: s(:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='z_svd'
  double complex,allocatable :: b(:,:),work(:),uu(:,:),vv(:,:)
  double precision,allocatable :: ss(:),rwork(:)
  integer,allocatable :: iwork(:)
  integer :: m,n,mn,mx,lwork,lrwork,ierr,r,i,j
  double precision,external :: dnrm2
  m=size(a,1); n=size(a,2); mn=min(m,n); mx=max(m,n); lwork=-1; lrwork=mn*max(5*(mn+1),2*(mn+mx+1))
  allocate(work(1),rwork(lrwork),b(m,n), uu(m,mn),vv(mn,n),ss(mn))
  ! call zgesdd('s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,rwork,iwork,ierr)
  call zgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,rwork,ierr)
  lwork=int(real(work(1)))+1
  deallocate(work)
  allocate(work(lwork),iwork(8*mn), stat=ierr)
  if(ierr.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call zcopy(m*n,a,1,b,1)
  call zgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,rwork,ierr)
  !call zgesdd('s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,rwork,iwork,ierr)
  if(present(info))info=ierr
  deallocate(b,work,rwork,iwork)
  if(ierr.ne.0)then
   write(*,*)subnam,': zgesvd info: ',ierr
   if(ierr.lt.0)stop
   if (nan(a)) then
    write(*,*) subnam,': NaNs detected in the input array'
    stop
   end if 
   u=>null(); v=>null(); s=>null()
  else 
   r=chop(ss,tol,rmax,err)
   if(present(err))err=err/dnrm2(mn,ss,1)
   allocate(u(m,r),s(r),v(r,n))
   call zcopy(m*r,uu,1,u,1)
   call dcopy(r,ss,1,s,1)
   call z2submat(r,n,vv,mn,v)
  end if
  deallocate(uu,vv,ss)
 end subroutine 

 integer function chop(s,tol,rmax,err) result (r)
  implicit none
  double precision,intent(in) :: s(:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  double precision :: nrm,er,er2,bound
  double precision,external :: dnrm2
  r=size(s); er2=0.d0
  if(present(rmax))then
   if(rmax.lt.r)then
    er2=dot_product(s(rmax+1:r),s(rmax+1:r))
    r=rmax
   end if 
  end if 
  if(present(tol))then
   nrm=dnrm2(size(s),s,1)
   bound=tol*tol*nrm*nrm 
   er=er2+s(r)*s(r)
   do while(er.lt.bound)
    er2=er; r=r-1; er=er+s(r)*s(r)
   end do
  end if
  if(present(err))err=dsqrt(er2)
  return
 end function

 
! NORM2 
 double precision function norm2_d(a) result(nrm)
  implicit none
  double precision,intent(in)  :: a(:,:)
  character(len=*),parameter :: subnam='norm2_d'
  double precision,pointer :: u(:,:),v(:,:),s(:)
  integer :: info
  call d_svd(a,u,s,v,tol=1.d-3,info=info)
  if(info.ne.0)then;write(*,'(a,a,i10)')subnam,': svd info: ',info;stop;end if
  nrm=s(1)
  deallocate(u,v,s)
 end function 
 
 double precision function norm2p_d(a,x) result(nrm)
  implicit none
  double precision,intent(in)  :: a(:,:)
  double precision,intent(inout) :: x(size(a,2))
  character(len=*),parameter :: subnam='norm2p_d'
  double precision,allocatable :: y(:)
  double precision :: xnrm,ynrm
  integer :: m,n,i,info
  double precision,external :: dnrm2

  m=size(a,1); n=size(a,2)
  allocate(y(m),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  
  xnrm=dnrm2(n,x,1)
  if(xnrm.le.1.d-14)then
   write(*,*)subnam,': please provide non-zero x!'
   nrm=-1
   return
  end if 
  if(xnrm.ne.1.d0)call dscal(n,1.d0/xnrm,x,1)

  do i=1,32
   call dgemv('n',m,n,1.d0,a,m,x,1,0.d0,y,1)
   ynrm=dnrm2(m,y,1)
   if(ynrm.ne.1.d0)call dscal(m,1.d0/ynrm,y,1)
   call dgemv('t',m,n,1.d0,a,m,y,1,0.d0,x,1)
   xnrm=dnrm2(n,x,1)
   if(xnrm.ne.1.d0)call dscal(n,1.d0/xnrm,x,1)
   nrm=xnrm
  end do

  deallocate(y)
 end function 
end module      
