!
! Tensor train types and basic functionality
!----------------------------------------------
!
! Written by Dmitry Savostyanov, dmitry.savostyanov@gmail.com
! Contributed by Sergey Dolgov, sergey.v.dolgov@gmail.com
!
module tt_lib
 use ptype_lib
 use default_lib
 use mat_lib
 use trans_lib
 use time_lib
 use say_lib
 implicit none
 integer,parameter :: tt_size=2048

 type,public:: dtt            ! double precision tensor train
  integer :: l=1              ! index of the leftmost core
  integer :: m=0              ! index of the rightmost core
  integer :: n(tt_size)=0     ! mode sizes (storage for matrices)
  integer :: q(tt_size)=0     ! first mode sizes (for matrices)
  integer :: s(tt_size)=0     ! second mode sizes (for matrices)
  integer :: t=0              ! data type for future sparsity etc.
  integer :: r(0:tt_size)=0   ! TT ranks
  type(pointd3) :: u(tt_size) ! TT cores
 contains
  !procedure :: assign1 => dtt_assign,dtt_mul_dt,dtt_dealloc
  !generic :: assignment(=) => assign1
  !generic :: operator(*)   => dtt_mul_dt

   ! finalisation is not implemented in GCC<4.9
   ! and in GCC<5 it is buggy
   ! So, turn it off and NEVER use norm(a+b) or stuff like that
   
  !final :: dtt_dealloc 
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
 contains
  !procedure :: ztt_assign,ztt_mul_zt,ztt_dealloc
  !generic :: assignment(=) => ztt_assign
  !generic :: operator(*)   => ztt_mul_zt
  !final :: ztt_dealloc 
 end type

 interface tijk         ! particular element from tt-vector as function of index
  module procedure dtt_ijk,ztt_ijk
 end interface
 interface value        ! value of tt-vector as function of its coordinates on [0:1]
  module procedure dtt_value,dtt_value0,ztt_value,ztt_value0
 end interface
 interface elem         ! particular element from tt-vector as subroutine
  module procedure dtt_elem
 end interface
 interface sumall       ! sum of all elements of tt-vector
  module procedure dtt_sumall,ztt_sumall
 end interface
 interface numel        ! number of elements of tt-vector
  module procedure dtt_numel,ztt_numel
 end interface
 interface norm         ! frobenius norm of tt-vector
  module procedure dtt_norm,ztt_norm
 end interface
 interface lognrm       ! log10 frobenius norm of tt-vector
  module procedure dtt_lognrm,ztt_lognrm
 end interface
 interface erank        ! effective rank of tt-vector
  module procedure dtt_rank,ztt_rank
 end interface
 interface memory          ! memory to keep all cores
  module procedure dtt_mem,ztt_mem
 end interface

 interface alloc        ! allocation and deallocation of the cores
  module procedure dtt_alloc,ztt_alloc
 end interface
 interface dealloc
  module procedure dtt_dealloc,ztt_dealloc
 end interface
 interface ready
  module procedure dtt_ready,ztt_ready
 end interface

 interface ort          ! tt-orthogonalization
  module procedure dtt_ort, ztt_ort
 end interface
 interface svd          ! tt-svd for rank truncation (compression, rounding)
  module procedure dtt_svd,dtt_svd0, ztt_svd,ztt_svd0
 end interface
 interface dot_product
  module procedure dtt_dot,ztt_dot
 end interface

 interface say    ! output some information about tt-vector
  module procedure dtt_say, ztt_say
 end interface
 
 interface ones
  module procedure dtt_ones, ztt_ones
 end interface
 interface zeros
  module procedure dtt_zeros, ztt_zeros
 end interface

 interface copy
  module procedure dtt_copy,ztt_copy
 end interface
 interface assignment (=)
  module procedure dtt_assign, ztt_assign, ztt_dtt_assign
 end interface
 interface operator (+)
  module procedure dtt_plus_dtt, ztt_plus_ztt, dtt_plus_d
 end interface
 interface operator (*)
  module procedure dtt_mul_dt, ztt_mul_zt
 end interface

contains


! ORT
 subroutine dtt_ort(arg)
  ![tt] ortogonalize from left
  implicit none
  type(dtt),intent(inout),target :: arg
  character(len=*),parameter :: subnam='dtt_ort'
  integer :: l,m,k,i,j,lwork,info,nn,rr,mn,mm,kk
  integer,pointer :: r(:),n(:)
  double precision,allocatable :: work(:),tau(:),mat(:),u(:)
  double precision :: err,nrm,lognrm
  double precision,external :: dnrm2

  l=arg%l; m=arg%m
  if(m.lt.l)return
  r=>arg%r; n=>arg%n

  nn=maxval(n(l:m)); rr=maxval(r(l-1:m))
  lwork=128*nn*rr
  allocate(work(lwork),tau(nn*rr), mat(rr*nn*rr),u(rr*nn*rr), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': no memory';stop;endif

  lognrm=0.d0
  do k=l,m-1
   mm=r(k-1)*n(k); nn=r(k); mn=min(mm,nn); kk=n(k+1)*r(k+1)
   call dcopy(mm*nn, arg%u(k)%p,1,u,1)
   call dgeqrf(mm,nn, u,mm,tau,work,lwork,info)
   if(info.ne.0)then; write(*,*) subnam,': dgeqrf info: ',info; stop; end if
   do j=1,nn
    forall(i=1:min(j,mm))    mat(i+(j-1)*mn)=u(i+(j-1)*mm)
    forall(i=min(j,mm)+1:mn) mat(i+(j-1)*mn)=0.d0
   end do
   nrm=dnrm2(mn*nn,mat,1)
   if(nrm.ne.0.d0)then
    call dscal(mn*nn,1.d0/nrm,mat,1)
    lognrm=lognrm+dlog(nrm)
   endif

   call dorgqr(mm,mn,mn,u,mm,tau,work,lwork,info)
   if(info.ne.0)then; write(*,*) subnam,': dorgqr info: ',info; stop; end if
!
!   nrm=dnrm2(mm*nn,arg%u(k)%p,1)
!   call dgemm('n','n',mm,nn,mn,1.d0,u,mm,mat,mn,-1.d0,arg%u(k)%p,mm)
!   err=dnrm2(mm*nn,arg%u(k)%p,1)
!   if(err.gt.1.d-10*nrm)then; write(*,*)subnam,': qr error: m,n: ',mm,nn;stop;endif
!
   call dcopy(mm*mn, u,1,arg%u(k)%p,1)
   call dgemm('n','n',mn,kk,nn,1.d0,mat,mn,arg%u(k+1)%p,nn,0.d0,u,mn)
   if(r(k).ne.mn)then
    call dcopy(mm*mn, arg%u(k)%p,1,mat,1)
    deallocate(arg%u(k)%p,arg%u(k+1)%p)
    r(k)=mn
    allocate(arg%u(k)%p(r(k-1),n(k),r(k)),arg%u(k+1)%p(r(k),n(k+1),r(k+1)))
    call dcopy(mm*mn, mat,1,arg%u(k)%p,1)
   end if
   call dcopy(mn*kk, u,1,arg%u(k+1)%p,1)
  end do
  deallocate(work,tau,mat,u)

  nrm=dnrm2(r(m-1)*n(m)*r(m),arg%u(m)%p,1)
  if(nrm.ne.0.d0)then
   call dscal(r(m-1)*n(m)*r(m),1.d0/nrm,arg%u(m)%p,1)
   lognrm=lognrm+dlog(nrm)
  endif

  lognrm=lognrm/(m-l+1)
  nrm=dexp(lognrm)
  do k=l,m
   call dscal(r(k-1)*n(k)*r(k),nrm,arg%u(k)%p,1)
  end do
 end subroutine
 subroutine ztt_ort(arg)
  ![tt] ortogonalize from left
  implicit none
  type(ztt),intent(inout),target :: arg
  character(len=*),parameter :: subnam='ztt_ort'
  integer :: l,m,k,i,j,lwork,info,nn,rr,mn,mm,kk
  integer,pointer :: r(:),n(:)
  double complex,allocatable :: work(:),tau(:),mat(:),u(:)
  double precision :: err,nrm,lognrm
  double precision,external :: dznrm2

  l=arg%l; m=arg%m
  if(m.lt.l)return
  r=>arg%r; n=>arg%n

  nn=maxval(n(l:m)); rr=maxval(r(l-1:m))
  lwork=128*nn*rr
  allocate(work(lwork),tau(nn*rr), mat(rr*nn*rr),u(rr*nn*rr), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': no memory';stop;endif

  do k=l,m-1
   mm=r(k-1)*n(k); nn=r(k); mn=min(mm,nn); kk=n(k+1)*r(k+1)
   call zcopy(mm*nn, arg%u(k)%p,1,u,1)
   call zgeqrf(mm,nn, u,mm,tau,work,lwork,info)
   if(info.ne.0)then; write(*,*) subnam,': zgeqrf info: ',info; stop; end if
   do j=1,nn
    forall(i=1:min(j,mm))    mat(i+(j-1)*mn)=u(i+(j-1)*mm)
    forall(i=min(j,mm)+1:mn) mat(i+(j-1)*mn)=(0.d0,0.d0)
   end do
   nrm=dznrm2(mn*nn,mat,1)
   if(nrm.ne.0.d0)then
    call zdscal(mn*nn,1.d0/nrm,mat,1)
    lognrm=lognrm+dlog(nrm)
   endif

   call zungqr(mm,mn,mn,u,mm,tau,work,lwork,info)
   if(info.ne.0)then; write(*,*) subnam,': zungqr info: ',info; stop; end if
!
!   nrm=dznrm2(mm*nn,arg%u(k)%p,1)
!   call zgemm('n','n',mm,nn,mn,(1.d0,0.d0),u,mm,mat,mn,(-1.d0,0.d0),arg%u(k)%p,mm)
!   err=dznrm2(mm*nn,arg%u(k)%p,1)
!   if(err.gt.1.d-10*nrm)then
!    write(*,*)subnam,': qr error: m,n: ',mm,nn
!    write(*,*)subnam,': err: ',err,' nrm ',nrm
!    stop
!   endif
!
   call zcopy(mm*mn, u,1,arg%u(k)%p,1)
   call zgemm('n','n',mn,kk,nn,(1.d0,0.d0),mat,mn,arg%u(k+1)%p,nn,(0.d0,0.d0),u,mn)
   if(r(k).ne.mn)then
    call zcopy(mm*mn, arg%u(k)%p,1,mat,1)
    deallocate(arg%u(k)%p,arg%u(k+1)%p,stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot deallocate data';stop;endif
    arg%r(k)=mn
    allocate(arg%u(k)%p(r(k-1),n(k),r(k)),arg%u(k+1)%p(r(k),n(k+1),r(k+1)))
    call zcopy(mm*mn, mat,1,arg%u(k)%p,1)
   end if
   call zcopy(mn*kk, u,1,arg%u(k+1)%p,1)
  end do
  deallocate(work,tau,mat,u)
  
  nrm=dznrm2(r(m-1)*n(m)*r(m),arg%u(m)%p,1)
  if(nrm.ne.0.d0)then
   call zdscal(r(m-1)*n(m)*r(m),1.d0/nrm,arg%u(m)%p,1)
   lognrm=lognrm+dlog(nrm)
  endif

  lognrm=lognrm/(m-l+1)
  nrm=dexp(lognrm)
  do k=l,m
   call zdscal(r(k-1)*n(k)*r(k),nrm,arg%u(k)%p,1)
  end do
 end subroutine

 subroutine dtt_normalize(arg)
  ![tt] ortogonalize from left and set norm=1
  implicit none
  type(dtt),intent(inout) :: arg
  character(len=*),parameter :: subnam='dtt_normalize'
  integer :: l,m,k
  integer :: r(0:tt_size),n(tt_size)
  double precision :: nrm
  double precision,external :: dnrm2
  call dtt_ort(arg)
  l=arg%l;m=arg%m; r=arg%r;n=arg%n
  nrm=dnrm2(r(m-1)*n(m)*r(m),arg%u(m)%p,1)
  do k=l,m
   call dscal(r(k-1)*n(k)*r(k),1.d0/nrm,arg%u(k)%p,1)
  end do
 end subroutine 
 subroutine ztt_normalize(arg)
  ![tt] ortogonalize from left and set norm=1
  implicit none
  type(ztt),intent(inout) :: arg
  character(len=*),parameter :: subnam='ztt_normalize'
  integer :: l,m,k
  integer :: r(0:tt_size),n(tt_size)
  double precision :: nrm
  double precision,external :: dznrm2
  call ztt_ort(arg)
  l=arg%l;m=arg%m; r=arg%r;n=arg%n
  nrm=dznrm2(r(m-1)*n(m)*r(m),arg%u(m)%p,1)
  do k=l,m
   call zdscal(r(k-1)*n(k)*r(k),1.d0/nrm,arg%u(k)%p,1)
  end do
 end subroutine

! SVD
 subroutine dtt_svd(arg,tol,rmax)
  implicit none
  type(dtt),intent(inout),target :: arg
  double precision,intent(in) :: tol
  integer,intent(in),optional :: rmax
  character(len=*),parameter :: subnam='dtt_svd'
  integer :: l,m,k,i,j,info,nn,rr,mn,mm,kk
  integer,pointer :: r(:),n(:)
  double precision,pointer,contiguous :: mat(:,:),u(:,:),s(:),v(:,:),tmp(:)
  double precision :: err,nrm,lognrm
  double precision,external :: dnrm2

  l=arg%l; m=arg%m
  if(m.le.l)return
  r=>arg%r; n=>arg%n

  nn=maxval(n(l:m)); rr=maxval(r(l-1:m))
  allocate(tmp(rr*nn*rr), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': no memory';stop;endif

  call dtt_ort(arg)
  lognrm=0.d0

  do k=m,l+1,-1
   mm=r(k-1); nn=n(k)*r(k); mn=min(mm,nn); kk=r(k-2)*n(k-1)
   mat(1:mm,1:nn) => arg%u(k)%p
   call svd(mat,u,s,v,tol,rmax,err,info)
   if(info.ne.0)then; write(*,*)subnam,': svd info: ',info; stop; endif
   rr=size(s); nrm=dnrm2(rr,s,1)
   if(nrm.ne.0.d0)then
    call dscal(rr,1.d0/nrm,s,1)
    lognrm=lognrm+dlog(nrm)
   end if
   do j=1,rr
    forall(i=1:mm) u(i,j)=u(i,j)*s(j)
   enddo

   call dgemm('n','n',kk,rr,mm,1.d0,arg%u(k-1)%p,kk,u,mm,0.d0,tmp,kk)
   if(r(k-1).ne.rr)then
    deallocate(arg%u(k-1)%p,arg%u(k)%p)
    r(k-1)=rr
    allocate(arg%u(k-1)%p(r(k-2),n(k-1),r(k-1)),arg%u(k)%p(r(k-1),n(k),r(k)))
   end if
   call dcopy(kk*rr,tmp,1,arg%u(k-1)%p,1)
   call dcopy(rr*nn,v,1,arg%u(k)%p,1)
   deallocate(u,s,v)
   nullify(mat)
  end do
  
  nrm=dnrm2(r(l-1)*n(l)*r(l),arg%u(l)%p,1)
  if(nrm.ne.0.d0)then
   call dscal(r(l-1)*n(l)*r(l),1.d0/nrm,arg%u(l)%p,1)
   lognrm=lognrm+dlog(nrm)
  endif

  lognrm=lognrm/(m-l+1)
  nrm=dexp(lognrm)
  do k=l,m
   call dscal(r(k-1)*n(k)*r(k),nrm,arg%u(k)%p,1)
  end do
  deallocate(tmp)
 end subroutine
 subroutine ztt_svd(arg,tol,rmax)
  implicit none
  type(ztt),intent(inout),target :: arg
  double precision,intent(in) :: tol
  integer,intent(in),optional :: rmax
  character(len=*),parameter :: subnam='ztt_svd'
  double complex,parameter :: zero=(0.d0,0.d0), one=(1.d0,0.d0)
  integer :: l,m,k,i,j,lwork,info,nn,rr,mn,mm,kk
  integer,pointer :: r(:),n(:)
  double precision,pointer :: s(:)
  double complex,pointer,contiguous :: mat(:,:),u(:,:),v(:,:),tmp(:)
  double precision :: err,nrm,lognrm
  double precision,external :: dnrm2,dznrm2

  l=arg%l; m=arg%m
  if(m.le.l)return
  r=>arg%r; n=>arg%n

  nn=maxval(n(l:m)); rr=maxval(r(l-1:m))
  allocate(tmp(rr*nn*rr), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': no memory';stop;endif

  call ztt_ort(arg)
  lognrm=0.d0

  do k=m,l+1,-1
   mm=r(k-1); nn=n(k)*r(k); mn=min(mm,nn); kk=r(k-2)*n(k-1)
   mat(1:mm,1:nn) => arg%u(k)%p
   call svd(mat,u,s,v,tol,rmax,err,info)
   if(info.ne.0)then; write(*,*)subnam,': svd info: ',info; stop; end if
   rr=size(s); nrm=dnrm2(rr,s,1)
   if(nrm.ne.0.d0)then
    call dscal(rr,1.d0/nrm,s,1)
    lognrm=lognrm+dlog(nrm)
   end if 
   do j=1,rr
    forall(i=1:mm) u(i,j)=u(i,j)*s(j)
   enddo

   call zgemm('n','n',kk,rr,mm,one,arg%u(k-1)%p,kk,u,mm,zero,tmp,kk)
   if(r(k-1).ne.rr)then
    deallocate(arg%u(k-1)%p,arg%u(k)%p)
    r(k-1)=rr
    allocate(arg%u(k-1)%p(r(k-2),n(k-1),r(k-1)),arg%u(k)%p(r(k-1),n(k),r(k)))
   end if
   call zcopy(kk*rr,tmp,1,arg%u(k-1)%p,1)
   call zcopy(rr*nn,v,1,arg%u(k)%p,1)
   deallocate(u,s,v)
   nullify(mat)
  end do
  
  nrm=dznrm2(r(l-1)*n(l)*r(l),arg%u(l)%p,1)
  if(nrm.ne.0.d0)then
   call zdscal(r(l-1)*n(l)*r(l),1.d0/nrm,arg%u(l)%p,1)
   lognrm=lognrm+dlog(nrm)
  endif

  lognrm=lognrm/(m-l+1)
  nrm=dexp(lognrm)
  do k=l,m
   call zdscal(r(k-1)*n(k)*r(k),nrm,arg%u(k)%p,1)
  end do
  deallocate(tmp)
 end subroutine

 subroutine dtt_svd0(n,a,tt,tol,rmax)
  implicit none
  integer,intent(in) :: n(:)
  double precision,intent(in) :: a(*)
  type(dtt),intent(inout),target :: tt
  double precision,intent(in) :: tol
  integer,intent(in),optional :: rmax
  character(len=*),parameter :: subnam='dtt_svd0'
  integer :: i,j,l,m,k,nn,mm,mn,info
  double precision,pointer,contiguous :: b(:,:),u(:,:),v(:,:),tmp(:)
  double precision,pointer :: s(:)
  integer,pointer :: r(:)
  double precision,external :: dnrm2
  
  l=tt%l; m=l+size(n)-1;tt%m=m;tt%n(l:m)=n; r=>tt%r; r(l-1)=1;r(m)=1
  nn=product(n)
  if(.not.8*nn.gt.0)then
   write(*,*)subnam,': nn too big: ',nn
   write(*,*)'kind(int):', kind(nn)
   stop
  endif
  allocate(tmp(nn),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory';stop;endif
  if(dnrm2(nn,a,1).eq.0.d0)then;write(*,*)subnam,': zero norm input';tt%r=0;call dealloc(tt);return;endif
  call dcopy(nn,a,1,tmp,1)
  
  do k=m,l+1,-1
   mm=r(l-1)*product(tt%n(l:k-1)); nn=tt%n(k)*r(k); mn=min(mm,nn)
   b(1:mm,1:nn)=>tmp
   call svd(b,u,s,v,tol,rmax,info=info)
   if(info.ne.0)then;write(*,*)subnam,': svd info: ',info;stop;endif
   r(k-1)=size(s)
   do j=1,r(k-1)
    forall(i=1:mm) u(i,j)=u(i,j)*s(j)
   enddo
   if(associated(tt%u(k)%p))deallocate(tt%u(k)%p)
   allocate(tt%u(k)%p(r(k-1),tt%n(k),r(k)))
   call dcopy(r(k-1)*tt%n(k)*r(k),v,1,tt%u(k)%p,1)
   call dcopy(mm*r(k-1),u,1,tmp,1)
   deallocate(u,s,v)
  end do
  if(associated(tt%u(l)%p))deallocate(tt%u(l)%p)
  allocate(tt%u(l)%p(r(l-1),tt%n(l),r(l)))
  call dcopy(r(l-1)*tt%n(l)*r(l),tmp,1,tt%u(l)%p,1)
  deallocate(tmp)
 end subroutine
 subroutine ztt_svd0(n,a,tt,tol,rmax)
  implicit none
  integer,intent(in) :: n(:)
  double complex,intent(in) :: a(*)
  type(ztt),intent(inout),target :: tt
  double precision,intent(in) :: tol
  integer,intent(in),optional :: rmax
  character(len=*),parameter :: subnam='ztt_svd0'
  integer :: i,j,l,m,k,nn,mm,mn,info
  double complex,pointer,contiguous :: b(:,:),u(:,:),v(:,:),tmp(:)
  double precision,pointer :: s(:)
  integer,pointer :: r(:)
  double precision,external :: dznrm2
  
  l=tt%l; m=l+size(n)-1;tt%m=m;tt%n(l:m)=n; r=>tt%r; r(l-1)=1;r(m)=1
  nn=product(n)
  if(.not.8*nn.gt.0)then;write(*,*)subnam,': nn too big: ',nn;stop;endif
  allocate(tmp(nn),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory';stop;endif
  if(dznrm2(nn,a,1).eq.0.d0)then;write(*,*)subnam,': zero norm input';tt%r=0;call dealloc(tt);return;endif
  call zcopy(nn,a,1,tmp,1)
  
  do k=m,l+1,-1
   mm=r(l-1)*product(tt%n(l:k-1)); nn=tt%n(k)*r(k); mn=min(mm,nn)
   b(1:mm,1:nn)=>tmp
   call svd(b,u,s,v,tol,rmax,info=info)
   if(info.ne.0)then;write(*,*)subnam,': svd info: ',info;stop;endif
   r(k-1)=size(s)
   do j=1,r(k-1)
    forall(i=1:mm) u(i,j)=u(i,j)*s(j)
   enddo
   if(associated(tt%u(k)%p))deallocate(tt%u(k)%p)
   allocate(tt%u(k)%p(r(k-1),tt%n(k),r(k)))
   call zcopy(r(k-1)*tt%n(k)*r(k),v,1,tt%u(k)%p,1)
   call zcopy(mm*r(k-1),u,1,tmp,1)
   deallocate(u,s,v)
   nullify(b)
  end do
  if(associated(tt%u(l)%p))deallocate(tt%u(l)%p)
  allocate(tt%u(l)%p(r(l-1),tt%n(l),r(l)))
  call zcopy(r(l-1)*tt%n(l)*r(l),tmp,1,tt%u(l)%p,1)
  deallocate(tmp)
 end subroutine



! GROUP
 subroutine dtt_group(arg,grp,side)
  !grp=[grp arg]
  implicit none
  type(dtt),intent(in),target :: arg
  type(dtt),intent(inout),target :: grp
  integer,intent(in),optional :: side
  character(len=*),parameter :: subnam='dtt_group'
  type(dtt) :: z
  integer,pointer :: r(:),q(:),n(:)
  integer :: l,m,k,sid,mm,ll,i,j

  if(arg%l.ne.grp%l .or. arg%m.ne.grp%m)then;write(*,*)subnam,': length mismatch';stop;endif
  l=arg%l; m=arg%m
  if(.not.all(arg%n(l:m)==grp%n(l:m)))then;write(*,*)subnam,': size mismatch';stop;endif
  n=>arg%n;r=>grp%r;q=>arg%r
  if(present(side))then;sid=side;else;if(r(l-1).ge.r(m))then;sid=0;else;sid=1;endif;endif

  z%l=l; z%m=m; z%n=n; z%r=0
  select case(sid)
   case(0)
    if(r(m).ne.q(m))then;write(*,*)subnam,': right border ranks mismatch:',r(m),q(m);stop;endif
    z%r(l-1:m-1)=r(l-1:m-1)+q(l-1:m-1); z%r(m)=r(m)
   case(1)
    if(r(l-1).ne.q(l-1))then;write(*,*)subnam,': left border ranks mismatch:',r(l-1),q(l-1);stop;endif
    z%r(l-1)=r(l-1); z%r(l:m)=r(l:m)+q(l:m)
   case default
    write(*,*)subnam,': illegal side:',sid; stop
  end select
  call alloc(z)

  if(sid.eq.1)then
   forall(j=1:r(l)) z%u(l)%p(:,:,     j)=grp%u(l)%p(:,:,j)
   forall(j=1:q(l)) z%u(l)%p(:,:,r(l)+j)=arg%u(l)%p(:,:,j)
   ll=l+1;mm=m
  endif
  if(sid.eq.0)then
   forall(i=1:r(m-1)) z%u(m)%p(       i,:,:)=grp%u(m)%p(i,:,:)
   forall(i=1:q(m-1)) z%u(m)%p(r(m-1)+i,:,:)=arg%u(m)%p(i,:,:)
   ll=l;mm=m-1
  endif

  do k=ll,mm
   z%u(k)%p=0.d0
   forall(i=1:r(k-1),j=1:r(k)) z%u(k)%p(       i,:,     j)=grp%u(k)%p(i,:,j)
   forall(i=1:q(k-1),j=1:q(k)) z%u(k)%p(r(k-1)+i,:,r(k)+j)=arg%u(k)%p(i,:,j)
  end do
  grp=z
  call dealloc(z)
 end subroutine

 subroutine ztt_group(arg,grp,side)
  !grp=[grp arg]
  implicit none
  type(ztt),intent(in),target :: arg
  type(ztt),intent(inout),target :: grp
  integer,intent(in),optional :: side
  character(len=*),parameter :: subnam='ztt_group'
  type(ztt) :: z
  integer,pointer :: r(:),q(:),n(:)
  integer :: l,m,k,sid,mm,ll,i,j

  if(arg%l.ne.grp%l .or. arg%m.ne.grp%m)then;write(*,*)subnam,': length mismatch';stop;endif
  l=arg%l; m=arg%m
  if(.not.all(arg%n(l:m)==grp%n(l:m)))then;write(*,*)subnam,': size mismatch';stop;endif
  n=>arg%n;r=>grp%r;q=>arg%r
  if(present(side))then;sid=side;else;if(r(l-1).ge.r(m))then;sid=0;else;sid=1;endif;endif

  z%l=l; z%m=m; z%n=n; z%r=0
  select case(sid)
   case(0)
    if(r(m).ne.q(m))then;write(*,*)subnam,': right border ranks mismatch:',r(m),q(m);stop;endif
    z%r(l-1:m-1)=r(l-1:m-1)+q(l-1:m-1); z%r(m)=r(m)
   case(1)
    if(r(l-1).ne.q(l-1))then;write(*,*)subnam,': left border ranks mismatch:',r(l-1),q(l-1);stop;endif
    z%r(l-1)=r(l-1); z%r(l:m)=r(l:m)+q(l:m)
   case default
    write(*,*)subnam,': illegal side:',sid; stop
  end select
  call alloc(z)

  if(sid.eq.1)then
   forall(j=1:r(l)) z%u(l)%p(:,:,     j)=grp%u(l)%p(:,:,j)
   forall(j=1:q(l)) z%u(l)%p(:,:,r(l)+j)=arg%u(l)%p(:,:,j)
   ll=l+1;mm=m
  endif
  if(sid.eq.0)then
   forall(i=1:r(m-1)) z%u(m)%p(       i,:,:)=grp%u(m)%p(i,:,:)
   forall(i=1:q(m-1)) z%u(m)%p(r(m-1)+i,:,:)=arg%u(m)%p(i,:,:)
   ll=l;mm=m-1
  endif

  do k=ll,mm
   z%u(k)%p=0.d0
   forall(i=1:r(k-1),j=1:r(k)) z%u(k)%p(       i,:,     j)=grp%u(k)%p(i,:,j)
   forall(i=1:q(k-1),j=1:q(k)) z%u(k)%p(r(k-1)+i,:,r(k)+j)=arg%u(k)%p(i,:,j)
  end do
  grp=z
  call dealloc(z)
 end subroutine



! TIJK ELEM VALUE
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

 subroutine dtt_elem(arg,ind,a)
  implicit none
  type(dtt),intent(in):: arg
  integer,intent(in) :: ind(:)
  double precision,intent(out) :: a(*)
  character(len=*),parameter :: subnam='dtt_elem'
  integer :: info,i,l,m,n(tt_size),r(0:tt_size)
  double precision,pointer :: x(:,:),y(:,:),z(:,:)
  l=arg%l;m=arg%m;n=arg%n;r=arg%r
  if(any(ind(1:m-l+1)<=0).or.any(ind(1:m-l+1)>n(l:m)))then;write(*,*)subnam,': wrong index: ',ind;stop;endif
  allocate(x(r(l-1),r(l)),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': allocation error: ',info;stop;endif
  x=arg%u(l)%p(:,ind(1),:)
  do i=l+1,m
   allocate(y(r(i-1),r(i)),z(r(l-1),r(i)),stat=info)
   if(info.ne.0)then;write(*,*)subnam,': allocation error: ',info;stop;endif
   y=arg%u(i)%p(:,ind(i-l+1),:)
   z=matmul(x,y)
   deallocate(x,y); x=>z; nullify(z)
  end do
  call dcopy(r(l-1)*r(m),x,1,a,1)
  deallocate(x)
 end subroutine

 pure double precision function dtt_value(arg,x) result (val)
  implicit none
  type(dtt),intent(in) :: arg
  double precision,intent(in) :: x(:)
  integer :: l,m,r(0:tt_size),n(tt_size), id,dd,i,j,pos,mm,ind(tt_size)
  double precision :: xx
  ind=0; val=0.d0
  l=arg%l; m=arg%m; r=arg%r; n=arg%n; dd=size(x)
  if(l.gt.m)return
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
  val=dtt_ijk(arg,ind)
  !write(*,'(3f20.12,1x,e20.12)') x,val
  !write(*,'(127i1)')(mod(i,10),i=1,127)
  !write(*,'(127i1)')ind
 end function
 pure double complex function ztt_value(arg,x) result (val)
  implicit none
  type(ztt),intent(in) :: arg
  double precision,intent(in) :: x(:)
  integer :: l,m,r(0:tt_size),n(tt_size), id,dd,i,j,pos,mm,ind(tt_size)
  double precision :: xx
  ind=0;val=0.d0
  l=arg%l; m=arg%m; r=arg%r; n=arg%n; dd=size(x)
  if(l.gt.m)return
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
  val=ztt_ijk(arg,ind)
 end function

 pure double precision function dtt_value0(arg,x) result (val)
  implicit none
  type(dtt),intent(in) :: arg
  double precision,intent(in) :: x
  double precision :: xx(1)
  xx=x; val=dtt_value(arg,xx)
 end function
 pure double complex function ztt_value0(arg,x) result (val)
  implicit none
  type(ztt),intent(in) :: arg
  double precision,intent(in) :: x
  double precision :: xx(1)
  xx=x; val=ztt_value(arg,xx)
 end function

! SUMALL
 pure double precision function dtt_sumall(arg) result(val)
  implicit none
  type(dtt),intent(in):: arg
  character(len=*),parameter :: subnam='dtt_sumall'
  integer :: info,i,j,p,q,l,m,n(tt_size),r(0:tt_size)
  double precision,pointer :: x(:,:),y(:,:),z(:,:)
  l=arg%l;m=arg%m;n=arg%n;r=arg%r; val=0.d0
  if(r(l-1).gt.1 .or. r(m).gt.1) then; val=-1.d0;return;endif !write(*,*)subnam,': matrix-valued output does not fit!'
  allocate(x(r(l-1),r(l)),stat=info)
  if(info.ne.0)then;val=-2.d0;return;endif  !write(*,*)subnam,': allocation error: ',info
  forall(i=1:r(l-1),j=1:r(l))x(i,j)=0.d0
  do p=1,n(l); forall(i=1:r(l-1),j=1:r(l))x(i,j)=x(i,j)+arg%u(l)%p(i,p,j); enddo
  do q=l+1,m
   allocate(y(r(q-1),r(q)),z(r(l-1),r(q)),stat=info)
   if(info.ne.0)then;val=-3.d0;return;endif  !write(*,*)subnam,': allocation error: ',info
   forall(i=1:r(q-1),j=1:r(q))y(i,j)=0.d0
   do p=1,n(q); forall(i=1:r(q-1),j=1:r(q))y(i,j)=y(i,j)+arg%u(q)%p(i,p,j); enddo
   z=matmul(x,y)
   deallocate(x,y); x=>z; nullify(z)
  end do
  val=x(1,1)
  deallocate(x)
 end function
 pure double complex function ztt_sumall(arg) result(val)
  implicit none
  type(ztt),intent(in):: arg
  character(len=*),parameter :: subnam='ztt_sumall'
  integer :: info,i,j,p,q,l,m,n(tt_size),r(0:tt_size)
  double complex,pointer :: x(:,:),y(:,:),z(:,:)
  l=arg%l;m=arg%m;n=arg%n;r=arg%r; val=(0.d0,0.d0)
  if(r(l-1).gt.1 .or. r(m).gt.1) then; val=-1.d0;return;endif !write(*,*)subnam,': matrix-valued output does not fit!'
  allocate(x(r(l-1),r(l)),stat=info)
  if(info.ne.0)then;val=-2.d0;return;endif  !write(*,*)subnam,': allocation error: ',info
  forall(i=1:r(l-1),j=1:r(l))x(i,j)=(0.d0,0.d0)
  do p=1,n(l); forall(i=1:r(l-1),j=1:r(l))x(i,j)=x(i,j)+arg%u(l)%p(i,p,j); enddo
  do q=l+1,m
   allocate(y(r(q-1),r(q)),z(r(l-1),r(q)),stat=info)
   if(info.ne.0)then;val=-3.d0;return;endif  !write(*,*)subnam,': allocation error: ',info
   forall(i=1:r(q-1),j=1:r(q))y(i,j)=(0.d0,0.d0)
   do p=1,n(q); forall(i=1:r(q-1),j=1:r(q))y(i,j)=y(i,j)+arg%u(q)%p(i,p,j); enddo
   z=matmul(x,y)
   deallocate(x,y); x=>z; nullify(z)
  end do
  val=x(1,1)
 end function

! NUMEL
 double precision function dtt_numel(arg) result (s)
  implicit none
  type(dtt),intent(in) :: arg
  integer :: l,m,i
  s=0.d0; l=arg%l; m=arg%m; if(l.gt.m)return
  s=1.d0; do i=l,m;s=s*arg%n(i); enddo
  return
 end function
 double precision function ztt_numel(arg) result (s)
  implicit none
  type(ztt),intent(in) :: arg
  integer :: l,m,i
  s=0.d0; l=arg%l; m=arg%m; if(l.gt.m)return
  s=1.d0; do i=l,m;s=s*arg%n(i); enddo
  return
 end function


!ALLOC
 subroutine dtt_memchk(arg)
  implicit none
  type(dtt),intent(inout) :: arg
  integer :: i
  character(len=tt_size) :: a
  if(arg%m<arg%l)return
  do i=arg%l,arg%m
   if(associated(arg%u(i)%p))then
    if(size(arg%u(i)%p).gt.arg%r(i-1)*arg%n(i)*arg%r(i))then
     a(i:i)='>'
    else if (size(arg%u(i)%p).lt.arg%r(i-1)*arg%n(i)*arg%r(i))then
     a(i:i)='<'
    else
     a(i:i)='='
    endif
   else
    a(i:i)='-'
   end if
  end do
  write(*,'(a,i2,a,i2,2a)') 'dtt[',arg%l,':', arg%m,']: ',a(arg%l:arg%m)
 end subroutine
 subroutine ztt_memchk(arg)
  implicit none
  type(ztt),intent(inout) :: arg
  integer :: i
  character(len=tt_size) :: a
  if(arg%m<arg%l)return
  do i=arg%l,arg%m
   if(associated(arg%u(i)%p))then
    if(size(arg%u(i)%p).gt.arg%r(i-1)*arg%n(i)*arg%r(i))then
     a(i:i)='>'
    else if (size(arg%u(i)%p).lt.arg%r(i-1)*arg%n(i)*arg%r(i))then
     a(i:i)='<'
    else
     a(i:i)='='
    endif
   else
    a(i:i)='-'
   end if
  end do
  write(*,'(a,i2,a,i2,2a)') 'ztt[',arg%l,':', arg%m,']: ',a(arg%l:arg%m)
 end subroutine

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

! PLUS
 type(dtt) function dtt_plus_dtt(a,b) result(c)
  type(dtt),intent(in),target :: a,b
  integer :: i,j,k,l,m,p
  integer,pointer :: n(:),r(:),q(:)
  if(.not.(ready(a).and.ready(b)))return
  l=a%l; m=a%m; n=>a%n; r=>a%r; q=>b%r
  c%l=a%l; c%m=a%m; c%n=a%n; c%r=0; c%r(l-1)=a%r(l-1); c%r(m)=a%r(m); c%r(l:m-1)=a%r(l:m-1)+b%r(l:m-1)
  call alloc(c)
  forall(i=1:r(l-1),j=1:n(l),k=1:r(l)) c%u(l)%p(i,j,     k)=a%u(l)%p(i,j,k)
  forall(i=1:r(l-1),j=1:n(l),k=1:q(l)) c%u(l)%p(i,j,r(l)+k)=b%u(l)%p(i,j,k)
  do p=l+1,m-1
   forall(i=1:r(p-1),j=1:n(p),k=1:r(p)) c%u(p)%p(       i,j,     k)=a%u(p)%p(i,j,k)
   forall(i=1:r(p-1),j=1:n(p),k=1:q(p)) c%u(p)%p(       i,j,r(p)+k)=0.d0
   forall(i=1:q(p-1),j=1:n(p),k=1:r(p)) c%u(p)%p(r(p-1)+i,j,     k)=0.d0
   forall(i=1:q(p-1),j=1:n(p),k=1:q(p)) c%u(p)%p(r(p-1)+i,j,r(p)+k)=b%u(p)%p(i,j,k)
  end do
  forall(i=1:r(m-1),j=1:n(m),k=1:r(m)) c%u(m)%p(       i,j,k)=a%u(m)%p(i,j,k)
  forall(i=1:q(m-1),j=1:n(m),k=1:r(m)) c%u(m)%p(r(m-1)+i,j,k)=b%u(m)%p(i,j,k)
 end function 
 type(ztt) function ztt_plus_ztt(a,b) result(c)
  type(ztt),intent(in),target :: a,b
  integer :: i,j,k,l,m,p
  integer,pointer :: n(:),r(:),q(:)
  if(.not.(ready(a).and.ready(b)))return
  l=a%l; m=a%m; n=>a%n; r=>a%r; q=>b%r
  c%l=a%l; c%m=a%m; c%n=a%n; c%r=0; c%r(l-1)=a%r(l-1); c%r(m)=a%r(m); c%r(l:m-1)=a%r(l:m-1)+b%r(l:m-1)
  call alloc(c)
  forall(i=1:r(l-1),j=1:n(l),k=1:r(l)) c%u(l)%p(i,j,     k)=a%u(l)%p(i,j,k)
  forall(i=1:r(l-1),j=1:n(l),k=1:q(l)) c%u(l)%p(i,j,r(l)+k)=b%u(l)%p(i,j,k)
  do p=l+1,m-1
   forall(i=1:r(p-1),j=1:n(p),k=1:r(p)) c%u(p)%p(       i,j,     k)=a%u(p)%p(i,j,k)
   forall(i=1:r(p-1),j=1:n(p),k=1:q(p)) c%u(p)%p(       i,j,r(p)+k)=(0.d0,0.d0)
   forall(i=1:q(p-1),j=1:n(p),k=1:r(p)) c%u(p)%p(r(p-1)+i,j,     k)=(0.d0,0.d0)
   forall(i=1:q(p-1),j=1:n(p),k=1:q(p)) c%u(p)%p(r(p-1)+i,j,r(p)+k)=b%u(p)%p(i,j,k)
  end do
  forall(i=1:r(m-1),j=1:n(m),k=1:r(m)) c%u(m)%p(       i,j,k)=a%u(m)%p(i,j,k)
  forall(i=1:q(m-1),j=1:n(m),k=1:r(m)) c%u(m)%p(r(m-1)+i,j,k)=b%u(m)%p(i,j,k)
 end function 
 type(dtt) function dtt_plus_d(a,b) result(c)
  type(dtt),intent(in),target :: a
  double precision,intent(in) :: b
  integer :: i,j,k,l,m,p
  integer,pointer :: n(:),r(:)
  if(.not.(ready(a)))return
  l=a%l; m=a%m; n=>a%n; r=>a%r
  c%l=a%l; c%m=a%m; c%n=a%n; c%r=0; c%r(l-1)=a%r(l-1); c%r(m)=a%r(m); c%r(l:m-1)=a%r(l:m-1)+1
  call alloc(c)
  forall(i=1:r(l-1),j=1:n(l),k=1:r(l)) c%u(l)%p(i,j,     k)=a%u(l)%p(i,j,k)
  forall(i=1:r(l-1),j=1:n(l))          c%u(l)%p(i,j,r(l)+1)=b
  do p=l+1,m-1
   forall(i=1:r(p-1),j=1:n(p),k=1:r(p)) c%u(p)%p(       i,j,     k)=a%u(p)%p(i,j,k)
   forall(i=1:r(p-1),j=1:n(p))          c%u(p)%p(       i,j,r(p)+1)=0.d0
   forall(           j=1:n(p),k=1:r(p)) c%u(p)%p(r(p-1)+1,j,     k)=0.d0
   forall(           j=1:n(p))          c%u(p)%p(r(p-1)+1,j,r(p)+1)=1.d0
  end do
  forall(i=1:r(m-1),j=1:n(m),k=1:r(m)) c%u(m)%p(       i,j,k)=a%u(m)%p(i,j,k)
  forall(           j=1:n(m),k=1:r(m)) c%u(m)%p(r(m-1)+1,j,k)=1.d0
 end function 


! MUL
 type(dtt) function dtt_mul_dt(a,b) result(c)
  double precision,intent(in) :: a
  type(dtt),intent(in) :: b
  integer :: k,l,m
  l=b%l; m=b%m; c%l=l; c%m=m; c%n=b%n; c%r=b%r; call alloc(c)
  do k=l,m
   call dcopy(b%r(k-1)*b%n(k)*b%r(k),b%u(k)%p,1,c%u(k)%p,1)
  end do
  call dscal(b%r(l-1)*b%n(l)*b%r(l),a,c%u(l)%p,1)
 end function
 type(ztt) function ztt_mul_zt(a,b) result(c)
  double complex,intent(in) :: a
  type(ztt),intent(in) :: b
  integer :: k,l,m, p,i,q
  l=b%l; m=b%m; c%l=l; c%m=m; c%n=b%n; c%r=b%r; call alloc(c)
  do k=l,m
   !call zcopy(b%r(k-1)*b%n(k)*b%r(k),b%u(k)%p,1,c%u(k)%p,1)
   forall(p=1:b%r(k-1),i=1:b%n(k),q=1:b%r(k)) c%u(k)%p(p,i,q)=a*b%u(k)%p(p,i,q)
  end do
  !call zscal(b%r(l-1)*b%n(l)*b%r(l),a,c%u(l)%p,1)
 end function

! ASSIGN COPY
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
  integer :: k,l,m, p,i,q
  l=a%l;m=a%m
  b%l=l; b%m=m; b%n(l:m)=a%n(l:m); b%r(l-1:m)=a%r(l-1:m); call alloc(b)
  do k=l,m
   !call zcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(k)%p,1);
   forall(p=1:a%r(k-1),i=1:a%n(k),q=1:a%r(k)) b%u(k)%p(p,i,q)=a%u(k)%p(p,i,q)
  end do
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
  integer :: k,l,m,ll,mm, p,i,q
  l=a%l;m=a%m; ll=default(b%l,low); mm=ll-l+m
  b%l=ll; b%m=mm; b%n(ll:mm)=a%n(l:m); b%r(ll-1:mm)=a%r(l-1:m);
  if(.not.all(a%n(l:m)>0))return;if(.not.all(a%r(l-1:m)>0))return;call alloc(b)
  do k=l,m
   !call zcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(ll-l+k)%p,1)
   forall(p=1:a%r(k-1),i=1:a%n(k),q=1:a%r(k)) b%u(ll-l+k)%p(p,i,q)=a%u(k)%p(p,i,q)
  end do
 end subroutine

! NORM
 double precision function dtt_norm(arg,tol) result (nrm)
  implicit none
  type(dtt),intent(in) :: arg
  double precision,intent(in),optional :: tol
  integer :: l,m
  type(dtt) :: tmp
  double precision,external :: dnrm2
  l=arg%l;m=arg%m; nrm=0.d0
  tmp=arg
  if(present(tol))then
   call svd(tmp,tol)
   nrm=dnrm2(size(tmp%u(l)%p),tmp%u(l)%p,1)
  else
   call ort(tmp)
   nrm=dnrm2(size(tmp%u(m)%p),tmp%u(m)%p,1)
  end if
  nrm=nrm**(m-l+1)
  call dealloc(tmp)
 end function
 double precision function ztt_norm(arg,tol) result (nrm)
  implicit none
  type(ztt),intent(in) :: arg
  double precision,intent(in),optional :: tol
  integer :: l,m
  type(ztt) :: tmp
  double precision,external :: dznrm2
  l=arg%l;m=arg%m; nrm=0.d0
  tmp=arg
  if(present(tol))then
   call svd(tmp,tol)
   nrm=dznrm2(size(tmp%u(l)%p),tmp%u(l)%p,1)
  else
   call ort(tmp)
   nrm=dznrm2(size(tmp%u(m)%p),tmp%u(m)%p,1)
  end if
  nrm=nrm**(m-l+1)
  call dealloc(tmp)
 end function


 double precision function dtt_lognrm(arg,tol) result (nrm)
  implicit none
  type(dtt),intent(in) :: arg
  double precision,intent(in),optional :: tol
  integer :: l,m
  type(dtt) :: tmp
  double precision,external :: dnrm2
  l=arg%l;m=arg%m; nrm=0.d0
  tmp=arg
  if(present(tol))then
   call svd(tmp,tol)
   nrm=dnrm2(size(tmp%u(l)%p),tmp%u(l)%p,1)
  else
   call dtt_ort(tmp)
   nrm=dnrm2(size(tmp%u(m)%p),tmp%u(m)%p,1)
  end if
  call dealloc(tmp)
  nrm=(dlog(nrm)/dlog(10.d0))*(m-l+1)
 end function
 double precision function ztt_lognrm(arg,tol) result (nrm)
  implicit none
  type(ztt),intent(in) :: arg
  double precision,intent(in),optional :: tol
  integer :: l,m
  type(ztt) :: tmp
  double precision,external :: dznrm2
  l=arg%l;m=arg%m; nrm=0.d0
  tmp=arg
  if(present(tol))then
   call svd(tmp,tol)
   nrm=dznrm2(size(tmp%u(l)%p),tmp%u(l)%p,1)
  else
   call ztt_ort(tmp)
   nrm=dznrm2(size(tmp%u(m)%p),tmp%u(m)%p,1)
  end if
  call dealloc(tmp)
  nrm=(dlog(nrm)/dlog(10.d0))*(m-l+1)
 end function


! DOT
 double precision function dtt_dot(x,y) result(dot)
 implicit none
 type(dtt),intent(in) :: x,y
 character(len=*),parameter :: subnam='dtt_dot'
 integer :: i,l,m,info,rx(0:tt_size),ry(0:tt_size),n(tt_size)
 double precision, allocatable :: phi(:), res(:) 
 if(x%l.ne.y%l .or. x%m.ne.y%m)then;write(*,*)subnam,': dimensions not match';stop;endif
 if(.not.all(x%n(x%l:x%m)==y%n(y%l:y%m)))then;write(*,*)subnam,': sizes not match';stop;endif
 l=x%l;m=x%m; rx=x%r;ry=y%r; n=x%n

 allocate(phi(maxval(rx(l-1:m)*ry(l-1:m))), res(maxval(rx(l-1:m-1)*n(l:m)*ry(l:m))), stat=info)
 if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
 
 phi(1)=1.d0
 do i = l,m
    call dgemm('n','n', rx(i-1),n(i)*ry(i),ry(i-1), 1d0,phi,rx(i-1), y%u(i)%p,ry(i-1), 0d0,res,rx(i-1))
    call dgemm('t','n', rx(i),ry(i),rx(i-1)*n(i), 1d0,x%u(i)%p,rx(i-1)*n(i), res,rx(i-1)*n(i), 0d0,phi,rx(i))
 end do
 dot=phi(1)
 deallocate(phi,res)
 end function
 double complex function ztt_dot(x,y) result(dot)
 implicit none
 type(ztt),intent(in) :: x,y
 character(len=*),parameter :: subnam='ztt_dot'
 double complex,parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
 integer :: i,l,m,info,rx(0:tt_size),ry(0:tt_size),n(tt_size)
 double complex, allocatable :: phi(:), res(:) 
 if(x%l.ne.y%l .or. x%m.ne.y%m)then;write(*,*)subnam,': dimensions not match';stop;endif
 if(.not.all(x%n(x%l:x%m)==y%n(y%l:y%m)))then;write(*,*)subnam,': sizes not match';stop;endif
 l=x%l;m=x%m; rx=x%r;ry=y%r; n=x%n

 allocate(phi(maxval(rx(l-1:m)*ry(l-1:m))), res(maxval(rx(l-1:m-1)*n(l:m)*ry(l:m))), stat=info)
 if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif

 phi(1)=one
 do i = l,m
    call zgemm('n','n',rx(i-1),n(i)*ry(i),ry(i-1), one,phi,rx(i-1), y%u(i)%p,ry(i-1), zero,res,rx(i-1))
    call zgemm('c','n',rx(i),ry(i),rx(i-1)*n(i), one,x%u(i)%p,rx(i-1)*n(i), res,rx(i-1)*n(i), zero,phi,rx(i))
 end do
 dot=phi(1)
 deallocate(phi,res)
 end function

! SAY
 subroutine dtt_say(arg)
  implicit none
  type(dtt),intent(in) :: arg
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
 subroutine ztt_say(arg)
  implicit none
  type(ztt),intent(in) :: arg
  character(len=1)d
  write(*,'(a,i2,a,i4,a,f6.2)') 'ztt[',arg%l,':', arg%m,']: rank ',erank(arg)
  if(all(arg%r(arg%l:arg%m).le.100))then;d='3';else;d='4';endif
  if(all(arg%n(arg%l+1:arg%m)==arg%n(arg%l)))then
   write(*,'(a,1x,i'//d//',a)') 'n: ',arg%n(arg%l),' for all modes'
  else 
   write(*,'(a,1x,1024i'//d//')') 'n: ',arg%n(arg%l:arg%m)
  end if
  write(*,'(a,1024i'//d//')') 'r: ',arg%r(arg%l-1:arg%m)
 end subroutine

! RANK
 double precision function dtt_rank(arg) result (r)
  implicit none
  type(dtt),intent(in) :: arg
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
 double precision function ztt_rank(arg) result (r)
  implicit none
  type(ztt),intent(in) :: arg
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

! MEM
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

! READY
 logical function dtt_ready(arg) result(l)
  implicit none
  type(dtt),intent(in) :: arg
  integer :: i
  l=(arg%l .le. arg%m)
  if(.not.l)return
  l=all(arg%n(arg%l:arg%m)>0)
  if(.not.l)return
  l=all(arg%r(arg%l-1:arg%m)>0)
  if(.not.l)return
  do i=arg%l,arg%m; l=l.and.associated(arg%u(i)%p); enddo
  if(.not.l)return
  do i=arg%l,arg%m 
   l=l.and.(size(arg%u(i)%p,1).eq.arg%r(i-1))
   l=l.and.(size(arg%u(i)%p,2).eq.arg%n(i))
   l=l.and.(size(arg%u(i)%p,3).eq.arg%r(i))
  enddo 
  if(.not.l)return
  return
 end function
 logical function ztt_ready(arg) result(l)
  implicit none
  type(ztt),intent(in) :: arg
  integer :: i
  l=(arg%l .le. arg%m)
  if(.not.l)return
  l=all(arg%n(arg%l:arg%m)>0)
  if(.not.l)return
  l=all(arg%r(arg%l-1:arg%m)>0)
  if(.not.l)return
  do i=arg%l,arg%m; l=l.and.associated(arg%u(i)%p); enddo
  if(.not.l)return
  do i=arg%l,arg%m 
   l=l.and.(size(arg%u(i)%p,1).eq.arg%r(i-1))
   l=l.and.(size(arg%u(i)%p,2).eq.arg%n(i))
   l=l.and.(size(arg%u(i)%p,3).eq.arg%r(i))
  enddo 
  if(.not.l)return
  return
 end function

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

end module
