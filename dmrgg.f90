module dmrgg_lib
 use tt_lib
 use lr_lib 
 use rnd_lib
 use default_lib
 use time_lib
 use ptype_lib
 implicit none
contains
 subroutine dtt_dmrgg(arg,fun,par,accuracy,maxrank,mybonds,pivoting,neval,quad,tru)
  ! approximate [fun] in TT format using dmrg greedy cross interpolation
  !use mpi
  implicit none
  include 'mpif.h'
  type(dtt),intent(inout),target :: arg
  double precision,external :: fun
  double precision,intent(in),optional :: par(*)
  double precision,intent(in),optional :: accuracy
  integer,intent(in),optional :: maxrank
  integer,intent(in),optional :: mybonds(0:)
  integer,intent(in),optional :: pivoting
  integer(kind=8),intent(out),optional :: neval
  type(dtt),intent(in),optional :: quad
  double precision,intent(in),optional :: tru

  character(len=*),parameter :: subnam='dtt_dmrgg'
  integer,parameter :: lft=1,rgt=2,smin=8
  double precision :: small_element, small_pivot
  integer,pointer :: r(:),n(:)
  integer :: info,l,m,i,j,k,p,q,ii,jj,kk,pp,qq,nn,s,t,it,crs,dir,ijkq,ij,jk,kq,nlot,ilot,piv,snum
  integer :: me,nproc,they,stat(MPI_STATUS_SIZE), typed,ihave,ineed,strike, first,last
  integer(kind=8) :: nevalloc,nevalall
  double precision :: t1,t2,t3,amax,pivot,pivotmax,pivotmin,val,err, pivotmax_prev,val_prev
  double precision,allocatable :: a(:,:,:,:),b(:),c(:),d(:),bcol(:,:,:),brow(:,:,:),bcol1(:,:),brow1(:,:),acol1(:,:),arow1(:,:)
  real*8 :: bb(4),cc(4)
  integer,allocatable :: lot(:,:),shifts(:),own(:)
  type(dtt) :: col,row,ttqq
  type(pointd) :: inv(0:tt_size)
  type(pointi2) :: vip(0:tt_size)
  integer :: ind(tt_size),rr(0:tt_size),tape(4,0:tt_size),tmpp(4,0:tt_size)
  logical :: ready,done,haverow,havecol,skipcol,needsend,needrecv, upd(0:tt_size)
  character(len=2) :: sdir
  character(len=1024) :: str,stmp,sfmt
  integer,external :: idamax
  double precision,external :: dnrm2,ddot
  
  t1=timef()
  piv=default(3,pivoting)
  ! check how our double precision type compares to hard-coded MPI types
  select case(storage_size(1.d0))
   case(32);  typed=MPI_REAL4;   small_element= 5*epsilon(1.d0);  small_pivot=1.d-3; sfmt='e14.8'
   case(64);  typed=MPI_REAL8;   small_element=10*epsilon(1.d0);  small_pivot=1.d-5; sfmt='e20.14'
   case(128); typed=MPI_REAL16;  small_element=50*epsilon(1.d0);  small_pivot=1.d-7; sfmt='e39.33'
   case default
    if(me.eq.0)write(*,'(a,a,i5,a)')subnam,': unknown storage_size(1.d0): ',storage_size(1.d0),' guessing 64'
    typed=MPI_REAL8; small_element=5*epsilon(1.d0);  small_pivot=1.d-5; sfmt='e20.14'
  end select  
  
  ! initialise dimensions, mode sizes and bond ranks
  if(arg%l.gt.arg%m)then;write(*,*)subnam,': l,m: ',arg%l,arg%m;stop;endif
  l=arg%l; m=arg%m; r=>arg%r; n=>arg%n
  if(any(r(l-1:m).gt.1))then;if(me.eq.0)write(*,*)'arg came in non-empty, wiping it clear';endif
  r(l-1:m)=1; call alloc(arg)
 
  ! check mpi variables
  call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_size fail: ',info;stop;endif
  call mpi_comm_rank(MPI_COMM_WORLD,me,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_rank fail: ',info;stop;endif
  if(nproc.ge.m)then;if(me.eq.0)write(*,*)'nproc exceeds or equal dimension, cannot proceed';stop;endif
 
  ! distribute bonds (ranks) between procs
  allocate(own(0:nproc), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate own bonds';stop;endif
  if(present(mybonds))then
   own(0:nproc)=mybonds(0:nproc)
  else 
   call share(l,m-1,own)
  end if 
  !write(*,'(a,i3,a,32i4)')'[',me,']: own: ',own(0:nproc)
  
  ! allocating vip: i=vip(p)%p(1,r); j=vip(p)%p(2,r); k=vip(p)%p(3,r); q=vip(p)%p(4,r)
  allocate(shifts(0:nproc), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate init shifts';stop;endif
  do p=l-1,m 
   allocate(inv(p)%p(1),vip(p)%p(4,1), stat=info)
   if(info.ne.0)then;write(*,*)subnam,': cannot allocate init inv and vip';stop;endif
   inv(p)%p(1)=1.d0
  end do

  !if(me.eq.0)write(*,*) 'locating initial cross'
  snum=max(smin,nproc) ! number of shifts
  do p=0,nproc-1
   shifts(p) = int(dble(snum)*dble(p)/nproc)
  end do
  shifts(nproc) = snum
  !write(*,'(a,i2,a,32i4)')'[',me,']: shifts: ',shifts(0:nproc)
  nn=minval(n(l:m))
  nlot=nn*(shifts(me+1)-shifts(me))
  allocate(b(nlot), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate b';stop;endif
  ihave=shifts(me+1)-shifts(me)
  do s=0,ihave-1
!$OMP PARALLEL DO PRIVATE(p,ind)
   do k=1,nn
    forall(p=l:m)ind(p)=mod(k-1+(s+shifts(me))*(p-1),n(p))+1
    b(k+s*nn)=fun(m,ind,n,par)
   end do 
!$OMP END PARALLEL DO
  end do
  ilot=idamax(nlot,b,1); amax=dabs(b(ilot))
  nevalloc=nlot
  ilot=ilot+nn*shifts(me)
  deallocate(b)

  s=(ilot-1)/nn; k=mod(ilot-1,nn)+1; forall(p=l:m)ind(p)=mod(k-1+s*(p-1),n(p))+1
  !write(*,'(a,i2,a,e10.3,a,512i4)')'[',me,'] local max  ',amax,' at ',ind(l:m)
  if(nproc.gt.1)then
   bb(1)=amax;bb(2)=ilot
   call mpi_allreduce(bb,cc,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_WORLD,info)
   if(info.ne.0)then;write(*,*)'mpi maxloc error:',info;stop;endif
   amax=cc(1);ilot=int(cc(2))
  endif 
  s=(ilot-1)/nn; k=mod(ilot-1,nn)+1; forall(p=l:m)ind(p)=mod(k-1+s*(p-1),n(p))+1
  !write(*,'(a,i2,a,e10.3,a,512i4)')'[',me,'] global max  ',amax,' at ',ind(l:m)
  !write(*,'(a,2(1x,i4),a,e20.13)')'PIVOT ',ind(l:m),' : ',amax
  vip(l-1)%p(:,1)=(/1,1,1,1/)
  forall(p=l:m-1)vip(p)%p(:,1)=(/ 1, ind(p),ind(p+1), 1  /)
  vip(m)%p(:,1)=(/1,1,1,1/)
  
  !if(me.eq.0)write(*,*)' computing initial cross'
  do p=own(me),own(me+1)
!$OMP PARALLEL DO
   do j=1,n(p)
     arg%u(p)%p(1,j,1)=dmrgg_fun(1,j,ind(p+1),1, p, l,m,vip,r,n, fun, par)
   enddo
!$OMP END PARALLEL DO
   nevalloc=nevalloc+n(p)
   do j=1,n(p); amax=dmax1(amax,dabs(arg%u(p)%p(1,j,1))); enddo
   !write(*,'(a,i2,a,i2,a,50e10.3)')'[',me,']{',p,'}: fiber: ', arg%u(p)%p(1,:,1)
  end do 
  pivotmax_prev = amax
  do p=own(me),own(me+1)-1
   pivot=arg%u(p)%p(1,ind(p),1)
   inv(p)%p(1)=pivot
   !write(*,'(a,i2,a,i2,a,e20.15)')'[',me,']{',p,'}: cross entry: ', arg%u(p)%p(1,ind(p),1)
  end do
  
  !if(me.eq.0)write(*,*)'initialise cols and rows of approximation'
  col=arg; row=arg
  do p=own(me),own(me+1)-1
   call d2_lual(n(p),r(p),inv(p)%p,col%u(p)%p)
   call d2_luar(n(p),r(p),inv(p)%p,row%u(p+1)%p)
  end do

  if(present(quad))then
   !if(me.eq.0)write(*,*)'computing quadrature'
   val=1.d0
   do p=own(me),own(me+1)-1
    val=val*ddot(n(p),arg%u(p)%p,1,quad%u(p)%p,1) / inv(p)%p(1)
   end do
   if(me.eq.nproc-1)val=val*ddot(n(m),arg%u(m)%p,1,quad%u(m)%p,1)
   if(nproc.gt.1)then
    bb(1)=val
    call mpi_allreduce(bb,cc,1,MPI_REAL8,MPI_PROD,MPI_COMM_WORLD,info)
    if(info.ne.0)then;write(*,*)subnam,': mpi allreduce fail for val:',info;stop;endif
    val=cc(1)
   endif 
   val_prev=val
   !if(me.eq.0)write(*,'(a,e20.15)')'initial value:',val
  end if
   
  call mpi_reduce(nevalloc,nevalall,1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,info)
  if(info.ne.0)then;write(*,*)subnam,': mpi reduce neval failed: ',info;stop;endif

  !if(me.eq.0)write(*,*)'clearing the tape'
  forall(p=l-1:m)tape(:,p)=(/ -1, -1, -1, -1 /)
  forall(p=l-1:m)tmpp(:,p)=(/ -2, -2, -2, -2 /)
  forall(p=l-1:m)upd(p)=.false.
  
  !report initial results
  if(me.eq.0)then
   t2 = timef()
   write(str,'(i3,a2,a,f5.1,a,2f6.3,a)') &
    0,'::',' rank', erank(arg), &
    ' log10t,n ',dlog(t2-t1)/dlog(10.d0),dlog(dble(nevalall))/dlog(10.d0),' ...  .......'
   if(present(quad))then
    write(stmp,'(a,'//sfmt//')')' val ',val
    str=trim(str)//trim(stmp)
   end if
   write(*,'(a)')trim(str)
  end if 
  
  it=0; strike=0; ready=.false.; if(present(maxrank))ready=(it+1.ge.maxrank)
  do while(.not.ready) ! main loop increasing the ranks and hopefully accuracy
   it=it+1
   dir=2-mod(it,2); if(dir.eq.1)then;sdir='>>';else if(dir.eq.2)then;sdir='<<';else;sdir='??';endif
   rr(l-1:m)=r(l-1:m);  pivotmax=-1.d0;pivotmin=-1.d0

   do pp=1,own(me+1)-own(me)  ! sweep over processor own's bonds
    p=own(me)+pp-1; if(dir.eq.2)p=own(me+1)-pp
    !write(*,'(a,i2,a,i3,a2,a,i2)') '[,me,']: ' it,sdir,' bond ',p
    
    allocate(acol1(r(p-1),n(p)),arow1(n(p+1),r(p+1)),stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot allocate acol/arow';stop;endif

    if(piv.eq.-1)then ! full pivoting
     allocate(a(r(p-1),n(p),n(p+1),r(p+1)),b(r(p-1)*n(p)*n(p+1)*r(p+1)), bcol(r(p-1),n(p),r(p)),brow(r(p),n(p+1),r(p+1)), &
             stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate superblock';stop;endif
      ! compute long indices and form rhs
!$OMP PARALLEL DO private(i,j,k,q)
     do ijkq=1,r(p-1)*n(p)*n(p+1)*r(p+1)
      i=ijkq-1
      q=i/(r(p-1)*n(p)*n(p+1)); i=mod(i,r(p-1)*n(p)*n(p+1))
      k=i/(r(p-1)*n(p));        i=mod(i,r(p-1)*n(p))
      j=i/r(p-1);               i=mod(i,r(p-1))
      i=i+1;j=j+1;k=k+1;q=q+1
      a(i,j,k,q)=dmrgg_fun(i,j,k,q, p, l,m,vip,r,n, fun, par)
     end do
!$OMP END PARALLEL DO
     nevalloc=nevalloc+r(p-1)*n(p)*n(p+1)*r(p+1)
     ijkq=idamax(r(p-1)*n(p)*n(p+1)*r(p+1),a,1)-1
     q=ijkq/(r(p-1)*n(p)*n(p+1))+1; ijkq=mod(ijkq,r(p-1)*n(p)*n(p+1))
     k=ijkq/(r(p-1)*n(p))+1;        ijkq=mod(ijkq,r(p-1)*n(p))
     j=ijkq/r(p-1)+1;               i=mod(ijkq,r(p-1))+1
     amax=dmax1(amax,dabs(a(i,j,k,q)))

     ! compute residual and do full pivoting
     call dcopy(r(p-1)*n(p)*n(p+1)*r(p+1),a,1,b,1)
     call dgemm('n','n',r(p-1)*n(p),n(p+1)*r(p+1),r(p),-1.d0,col%u(p)%p,r(p-1)*n(p),row%u(p+1)%p,r(p),1.d0,b,r(p-1)*n(p))

     ijkq=idamax(r(p-1)*n(p)*n(p+1)*r(p+1),b,1)-1
     qq=ijkq/(r(p-1)*n(p)*n(p+1))+1; ijkq=mod(ijkq,r(p-1)*n(p)*n(p+1))
     kk=ijkq/(r(p-1)*n(p))+1;        ijkq=mod(ijkq,r(p-1)*n(p))
     jj=ijkq/r(p-1)+1;               ii=mod(ijkq,r(p-1))+1
     pivot=b(ii + r(p-1)*( jj-1 + n(p)*( kk-1 + n(p+1)*( qq-1 ) ) ) )
     !write(*,'(a,i2,a,i2,a,4i5,a,e10.3)')'[',me,']{',p,'} pivot at ',ii,jj,kk,qq,' : ',pivot

     ! copy column and row to containers
     forall(i=1:r(p-1),j=1:n(p))  acol1(i,j)=a(i,j,kk,qq)
     forall(k=1:n(p+1),q=1:r(p+1))arow1(k,q)=a(ii,jj,k,q)
     !do j=1,n(p);do i=1,r(p-1);   acol1(i,j)=a(i,j,kk,qq); enddo;enddo
     !do q=1,r(p+1);do k=1,n(p+1); arow1(k,q)=a(ii,jj,k,q); enddo;enddo

     ! deallocate superblocks
     deallocate(a,b,bcol,brow)

    else if(piv.ge.0)then ! partial pivoting on (r*n+n*r) random elements
     !nlot=r(p-1)*n(p)+n(p+1)*r(p+1)
     nlot=r(p-1)+n(p)+n(p+1)+r(p+1)
     !nlot=1

     !write(*,'(a,i2,a,i2,a,4i4)')'[',me,']{',p,'}: superblock ',r(p-1),n(p),n(p+1),r(p+1)
     !write(*,'(a,i2,a,i2,a,i4,a,i7)')'[',me,']{',p,'}: pivoting ',piv,' nlot ',nlot
     allocate(bcol1(r(p-1),n(p)),brow1(n(p+1),r(p+1)),lot(nlot,4),b(nlot), stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate bcol1/brow1 for pivoting';stop;endif

     ! distribute lottery weights for indices columns and rows
     forall(i=1:r(p-1),j=1:n(p))bcol1(i,j)=1.d0
     forall(k=1:n(p+1),q=1:r(p+1))brow1(k,q)=1.d0
     do s=1,r(p)
      i=vip(p)%p(1,s);j=vip(p)%p(2,s);k=vip(p)%p(3,s);q=vip(p)%p(4,s)
      bcol1(i,j)=0.d0; brow1(k,q)=0.d0
     end do
     ! get indices from lottery
     call lottery2(nlot,r(p-1)*n(p),n(p+1)*r(p+1),bcol1,brow1,lot)
     forall(i=1:nlot)lot(i,3)=lot(i,2)
     do ilot=1,nlot
      lot(ilot,2)=(lot(ilot,1)-1)/r(p-1)+1; lot(ilot,1)=mod(lot(ilot,1)-1,r(p-1))+1  ! ij
      lot(ilot,4)=(lot(ilot,3)-1)/n(p+1)+1; lot(ilot,3)=mod(lot(ilot,3)-1,n(p+1))+1  ! kq
     end do

     ! get matrix elements and residuals for these indices
!$OMP PARALLEL DO private(i,j,k,q)
     do ilot=1,nlot
      i=lot(ilot,1); j=lot(ilot,2); k=lot(ilot,3); q=lot(ilot,4)
      b(ilot)=dmrgg_fun(i,j,k,q, p, l,m,vip,r,n, fun, par)
      !write(*,'(a,i2,a,i2,a,4i4,a,4i4,a,e10.3)')'[',me,']{',p,'} elt: ',i,j,k,q,' : ',i,j,k,q,' : ',b(ilot)
     end do
!$OMP END PARALLEL DO
     nevalloc=nevalloc+nlot
     ilot=idamax(nlot,b,1); amax=dmax1(amax,dabs(b(ilot)))
     do ilot=1,nlot
      i=lot(ilot,1); j=lot(ilot,2); k=lot(ilot,3); q=lot(ilot,4)
      b(ilot)=b(ilot)-ddot(r(p),col%u(p)%p(i,j,1),r(p-1)*n(p),row%u(p+1)%p(1,k,q),1)
      !write(*,'(a,i2,a,i2,a,4i4,a,4i4,a,e10.3)')'[',me,']{',p,'} err: ',i,j,k,q,' : ',i,j,k,q,' : ',b(ilot)
     end do

     ! find maximum pivot
     ilot=idamax(nlot,b,1)
     ii=lot(ilot,1); jj=lot(ilot,2); kk=lot(ilot,3); qq=lot(ilot,4)
     pivot=b(ilot)
     !write(*,'(a,i2,a,i4,a,4i4,a,e20.13)') '[',me,']{',p,'} rnd pivot ',ii,jj,kk,qq,' : ',pivot

     done=.false.; havecol=.false.; haverow=.false.;
     if(piv.eq.0)then ! compute row and column
!$OMP PARALLEL private(i,j,k,q)
!$OMP DO
       do ij=0,r(p-1)*n(p)-1
        j=ij/r(p-1)+1; i=mod(ij,r(p-1))+1
        acol1(i,j)=dmrgg_fun(i,j,kk,qq, p, l,m,vip,r,n, fun, par)
       end do
!$OMP END DO
!$OMP DO
       do kq=0,n(p+1)*r(p+1)-1
        q=kq/n(p+1)+1; k=mod(kq,n(p+1))+1
        arow1(k,q)=dmrgg_fun(ii,jj,k,q, p, l,m,vip,r,n, fun, par)
       end do
!$OMP END DO
!$OMP END PARALLEL
      nevalloc=nevalloc+r(p-1)*n(p)+n(p+1)*r(p+1)
      done=.true.; havecol=.true.; haverow=.true.;
     end if

     ! rook pivoting to increase abs(pivot)
     crs=0; skipcol=(dir.eq.2)
     do while(.not.done)
      if(.not.skipcol)then
!$OMP PARALLEL DO private(i,j)
       do ij=0,r(p-1)*n(p)-1
        j=ij/r(p-1)+1; i=mod(ij,r(p-1))+1
        acol1(i,j)=dmrgg_fun(i,j,kk,qq, p, l,m,vip,r,n, fun, par)
       end do
!$OMP END PARALLEL DO
       nevalloc=nevalloc+r(p-1)*n(p)
       ij=idamax(r(p-1)*n(p),acol1,1)-1; j=ij/r(p-1)+1; i=mod(ij,r(p-1))+1; amax=dmax1(amax,dabs(acol1(i,j)))
       havecol=.true.; crs=crs+1; done=havecol.and.haverow.and.(crs.ge.2*piv)
       if(.not.done)then
        call dcopy(r(p-1)*n(p),acol1,1,bcol1,1)
        call dgemv('n',r(p-1)*n(p),r(p),-1.d0,col%u(p)%p,r(p-1)*n(p),row%u(p+1)%p(1,kk,qq),1,1.d0,bcol1,1)
        ij=idamax(r(p-1)*n(p),bcol1,1)-1; j=ij/r(p-1)+1; i=mod(ij,r(p-1))+1
        done=havecol.and.haverow.and. (i.eq.ii .and. j.eq.jj)
        ii=i;jj=j;pivot=bcol1(ii,jj)
       end if
       !write(*,'(a,i2,a,i4,a,4i4,a,e20.13)')'[',me,']{',p,'} col pivot ',ii,jj,kk,qq,' : ',pivot
      end if
      skipcol=.false.
      if(.not.done)then
!$OMP PARALLEL DO private(k,q)
       do kq=0,n(p+1)*r(p+1)-1
        q=kq/n(p+1)+1; k=mod(kq,n(p+1))+1
        arow1(k,q)=dmrgg_fun(ii,jj,k,q, p, l,m,vip,r,n, fun, par)
       end do
!$OMP END PARALLEL DO
       nevalloc=nevalloc+n(p+1)*r(p+1)
       kq=idamax(n(p+1)*r(p+1),arow1,1)-1; q=kq/n(p+1)+1; k=mod(kq,n(p+1))+1; amax=dmax1(amax,dabs(arow1(k,q)))
       haverow=.true.; crs=crs+1; done=havecol.and.haverow.and.(crs.ge.2*piv)
       if(.not.done)then
        call dcopy(n(p+1)*r(p+1),arow1,1,brow1,1)
        call dgemv('t',r(p),n(p+1)*r(p+1),-1.d0,row%u(p+1)%p,r(p),col%u(p)%p(ii,jj,1),r(p-1)*n(p),1.d0,brow1,1)
        kq=idamax(n(p+1)*r(p+1),brow1,1)-1; q=kq/n(p+1)+1; k=mod(kq,n(p+1))+1
        done=havecol.and.haverow.and. (k.eq.kk .and. q.eq.qq)
        qq=q;kk=k;pivot=brow1(kk,qq)
       end if
       !write(*,'(a,i2,a,i4,a,4i4,a,e20.13)')'[',me,']{',p,'} row pivot ',ii,jj,kk,qq,' : ',pivot
      end if
     end do
     deallocate(bcol1,brow1,lot,b, stat=info)
     if(info.ne.0)then;write(*,*)'fail to deallocate after pivoting';stop;endif
    else
     write(*,*) subnam,': unknown pivoting: ',piv; stop
    end if !(pivoting)
    !write(*,'(a,i2,a,i4,a,4i4,a,e20.13)')'[',me,']{',p,'} fin pivot ',ii,jj,kk,qq,' : ',pivot
    !write(*,'(a,2(1x,i4),a,e20.13)')'PIVOT ',jj,kk,' : ',pivot

    tape(:,p)=(/-1,-1,-1,-1/)
    upd(p)=(dabs(pivot).gt.small_element*amax) .and. (dabs(pivot).gt.small_pivot*pivotmax_prev) 
    if( upd(p) )then
     ! put new pivot in sets
     tape(:,p)=(/ ii,jj,kk,qq /)
     allocate(lot(4,r(p)+1),stat=info); if(info.ne.0)then;write(*,*)'allocate lot fail:',info;stop;endif
     forall(i=1:4,s=1:r(p))lot(i,s)=vip(p)%p(i,s)
     deallocate(vip(p)%p); allocate(vip(p)%p(4,r(p)+1),stat=info); if(info.ne.0)then;write(*,*)'reallocate vip fail:',info;stop;endif
     forall(i=1:4,s=1:r(p))vip(p)%p(i,s)=lot(i,s)
     vip(p)%p(:,r(p)+1)=tape(:,p)
     deallocate(lot)
     if(pivotmax.lt.0.d0)then;pivotmax=dabs(pivot);else;pivotmax=dmax1(pivotmax,dabs(pivot));endif
     if(pivotmin.lt.0.d0)then;pivotmin=dabs(pivot);else;pivotmin=dmin1(pivotmin,dabs(pivot));endif

      !write(*,'(a,i2,a,i2,a,4i5)')'[',me,']{',p,'}: reallocating superblock with sizes: ',r(p-1),n(p),n(p+1),r(p+1)
     allocate(bcol(r(p-1),n(p),r(p)+1),brow(r(p)+1,n(p+1),r(p+1)), bcol1(r(p-1),n(p)),brow1(n(p+1),r(p+1)),b(r(p)**2), stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot reallocate blocks for update';stop;endif

     !if(me.eq.0)write(*,*)'updating inverse..'
     call dcopy(r(p)**2,inv(p)%p,1,b,1)
     deallocate(inv(p)%p); allocate(inv(p)%p((r(p)+1)**2),stat=info)
     call dcopy(r(p)**2,b,1,inv(p)%p,1)
     if(info.ne.0)then;write(*,'(a,i2,a,i5)')'[',me,']: cannot reallocate inverse: ',p;stop;endif
     call dcopy(r(p),col%u(p)%p(ii,jj,1),r(p-1)*n(p),inv(p)%p(r(p)**2+1),1)
     call dcopy(r(p),row%u(p+1)%p(1,kk,qq),1,inv(p)%p(r(p)**2+r(p)+1),1)
     inv(p)%p((r(p)+1)**2)=pivot

     !write(*,'(a,i2,a,i2,a)') '[',me,']{',p,'}(arg): copying blocks'
     forall(i=1:r(p-1),j=1:n(p),s=1:r(p))bcol(i,j,s)=arg%u(p)%p(i,j,s)
     forall(i=1:r(p-1),j=1:n(p)         )bcol(i,j,r(p)+1)=acol1(i,j)
     forall(s=1:r(p),k=1:n(p+1),q=1:r(p+1))brow(s,k,q)=arg%u(p+1)%p(s,k,q)
     forall(         k=1:n(p+1),q=1:r(p+1))brow(r(p)+1,k,q)=arow1(k,q)
     !do s=1,r(p); do j=1,n(p); do i=1,r(p-1); bcol(i,j,s)=arg%u(p)%p(i,j,s); enddo; enddo; enddo
     !             do j=1,n(p); do i=1,r(p-1); bcol(i,j,r(p)+1)=acol1(i,j); enddo; enddo
     !do q=1,r(p+1);do k=1,n(p+1);do s=1,r(p); brow(s,k,q)=arg%u(p+1)%p(s,k,q);enddo;enddo;enddo
     !do q=1,r(p+1);do k=1,n(p+1);             brow(r(p)+1,k,q)=arow1(k,q);enddo;enddo
     deallocate(arg%u(p)%p, arg%u(p+1)%p); allocate(arg%u(p)%p(r(p-1),n(p),r(p)+1),arg%u(p+1)%p(r(p)+1,n(p+1),r(p+1)),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': not enough memory for new blocks';stop;endif
     call dcopy(r(p-1)*n(p)*(r(p)+1),bcol,1,arg%u(p)%p,1)
     call dcopy((r(p)+1)*n(p+1)*r(p+1),brow,1,arg%u(p+1)%p,1)

     !write(*,'(a,i2,a,i2,a)') '[',me,']{',p,'}(uv): updating factors'
     forall(i=1:r(p-1),j=1:n(p),s=1:r(p))bcol(i,j,s)=col%u(p)%p(i,j,s)
     forall(i=1:r(p-1),j=1:n(p)         )bcol(i,j,r(p)+1)=acol1(i,j)
     forall(s=1:r(p),k=1:n(p+1),q=1:r(p+1))brow(s,k,q)=row%u(p+1)%p(s,k,q)
     forall(         k=1:n(p+1),q=1:r(p+1))brow(r(p)+1,k,q)=arow1(k,q)
     call d2_lual(r(p-1)*n(p),  r(p)+1,inv(p)%p,bcol,from=r(p)+1)
     call d2_luar(n(p+1)*r(p+1),r(p)+1,inv(p)%p,brow,from=r(p)+1)
     deallocate(col%u(p)%p,row%u(p+1)%p); allocate(col%u(p)%p(r(p-1),n(p),r(p)+1),row%u(p+1)%p(r(p)+1,n(p+1),r(p+1)),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': not enough memory for new factors';stop;endif
     call dcopy(r(p-1)*n(p)*(r(p)+1),bcol,1,col%u(p)%p,1)
     call dcopy((r(p)+1)*n(p+1)*r(p+1),brow,1,row%u(p+1)%p,1)

     if(p.gt.own(me))then
      !write(*,'(a,i2,a,i2,a)') '[',me,']{',p,'}(v): updating left rows'
      call dcopy(r(p-1)*n(p)*r(p),row%u(p)%p,1,bcol,1)
      call dcopy(r(p-1)*n(p),arg%u(p)%p(1,1,r(p)+1),1,bcol1,1)
      call d2_luar(n(p),r(p-1),inv(p-1)%p,bcol1)
      call dcopy(r(p-1)*n(p),bcol1,1,bcol(1,1,r(p)+1),1)
      deallocate(row%u(p)%p); allocate(row%u(p)%p(r(p-1),n(p),r(p)+1), stat=info)
      if(info.ne.0)then;write(*,*)subnam,': not enough memory for left rows';stop;endif
      call dcopy(r(p-1)*n(p)*(r(p)+1),bcol,1,row%u(p)%p,1)
     end if
     if(p.lt.own(me+1)-1)then
      !write(*,'(a,i2,a,i2,a)') '[',me,']{',p,'}(u): updating right cols'
      forall(s=1:r(p),k=1:n(p+1),q=1:r(p+1))brow(s,k,q)=col%u(p+1)%p(s,k,q)
      forall(         k=1:n(p+1),q=1:r(p+1))brow1(k,q)=arg%u(p+1)%p(r(p)+1,k,q)
      call d2_lual(n(p+1),r(p+1),inv(p+1)%p,brow1)
      forall(         k=1:n(p+1),q=1:r(p+1))brow(r(p)+1,k,q)=brow1(k,q)
      deallocate(col%u(p+1)%p); allocate(col%u(p+1)%p(r(p)+1,n(p+1),r(p+1)),stat=info)
      if(info.ne.0)then;write(*,*)subnam,': not enough memory for right cols';stop;endif
      call dcopy((r(p)+1)*n(p+1)*r(p+1),brow,1,col%u(p+1)%p,1) 
     end if

     ! increase the rank and move on
     r(p)=r(p)+1
     deallocate(bcol,brow,bcol1,brow1,b,stat=info)
     if(info.ne.0)then;write(*,*)'fail to deallocate after update';stop;endif
    end if ! (update)
    deallocate(acol1,arow1)
   end do !(loop over my cores)

   if(nproc.gt.1)then
    forall(p=l-1:m)tmpp(:,p)=(/ -2, -2, -2, -2 /)
    !if(me.eq.0)write(*,*)'propagating tape to the right'
    ihave=own(me+1)-l
    ineed=own(me  )-l
    if(me.eq.0)then
     !write(*,'(a,i2,a,i6)')'[',me,'] sending tape right:',ihave
     call mpi_send(tape(1,l),4*ihave,MPI_INTEGER,me+1,rgt,MPI_COMM_WORLD,info)
     if(info.ne.0)then;write(*,'(a,i2,a,i6)')'[',me,'] mpi send fails for tape going right:',ihave;stop;endif
    else if(me.eq.nproc-1)then
     !write(*,'(a,i2,a,i6)')'[',me,'] receiving going right:',ineed
     call mpi_recv(tmpp(1,l),4*ineed,MPI_INTEGER,me-1,rgt,MPI_COMM_WORLD,stat,info)
     if(info.ne.0)then;write(*,'(a,i2,a,i6)')'[',me,'] mpi recv fails for tape going right:',ineed;stop;endif
    else
     !write(*,'(a,i2,a,2i6)')'[',me,'] sending and receiving tape going right:',ihave,ineed
     call mpi_sendrecv(tape(1,l),4*ihave,MPI_INTEGER,me+1,rgt, &
                       tmpp(1,l),4*ineed,MPI_INTEGER,me-1,rgt, &
                       MPI_COMM_WORLD,stat,info)
     if(info.ne.0)then;write(*,'(a,i2,a,2i6)')'[',me,'] mpi sendrecv fails for tape going right:',ihave,ineed;stop;endif
    end if
    !if(me.eq.0)write(*,*)'propagating tape to the left'
    ihave=m-own(me);   q=own(me)
    ineed=m-own(me+1); p=own(me+1)
    if(me.eq.0)then
     !write(*,'(a,i2,a,i6)')'[',me,'] receiving tape going left:',ineed
     call mpi_recv(tmpp(1,p),4*ineed,MPI_INTEGER,me+1,lft,MPI_COMM_WORLD,stat,info)
     if(info.ne.0)then;write(*,'(a,i2,a,i6)')'[',me,'] mpi recv fails for tape going left:',ineed;stop;endif
    else if(me.eq.nproc-1)then
     !write(*,'(a,i2,a,i6)')'[',me,'] sending tape to the left:',ihave
     call mpi_send(tape(1,q),4*ihave,MPI_INTEGER,me-1,lft,MPI_COMM_WORLD,info)
     if(info.ne.0)then;write(*,'(a,i2,a,i6)')'[',me,'] mpi send fails for tape going left:',ihave;stop;endif
    else
     !write(*,'(a,i2,a,2i6)')'[',me,'] sending and receiving tape going left:',ihave,ineed
     call mpi_sendrecv(tape(1,q),4*ihave,MPI_INTEGER,me-1,lft, &
                       tmpp(1,p),4*ineed,MPI_INTEGER,me+1,lft, &
                       MPI_COMM_WORLD,stat,info)
     if(info.ne.0)then;write(*,'(a,i2,a,2i6)')'[',me,'] mpi sendrecv fails for tape going left:',ihave,ineed;stop;endif
    end if

    !write(*,'(a,i2,a)')'[',me,'] reading from the tape and updating world info'
    tape(:,l-1:m)=tmpp(:,l-1:m)
    do p=l,m-1
     if(.not.(own(me).le.p .and. p.le.own(me+1)-1))then
      upd(p)=(tape(1,p).gt.0)
      if(upd(p))then
       allocate(lot(4,r(p)+1),stat=info); if(info.ne.0)then;write(*,*)'allocate lot fail:',info;stop;endif
       forall(i=1:4,s=1:r(p))lot(i,s)=vip(p)%p(i,s)
       deallocate(vip(p)%p); allocate(vip(p)%p(4,r(p)+1),stat=info); if(info.ne.0)then;write(*,*)'reallocate vip fail:',info;stop;endif
       forall(i=1:4,s=1:r(p))vip(p)%p(i,s)=lot(i,s)
       vip(p)%p(:,r(p)+1)=tape(:,p)
       deallocate(lot)
       r(p)=r(p)+1
      end if
     end if
    end do
   
   ! GLOBAL allreduce to communicate amax and pivots
    bb(1)=amax
    bb(2)=pivotmax
    if(pivotmin.gt.0.d0)then;bb(3)=-pivotmin;else;bb(3)=-999d9;endif
    call mpi_allreduce(bb,cc,3,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,info)
    if(info.ne.0)then;write(*,*)subnam,': mpi allreduce amax failed: ',info;stop;endif
    amax=cc(1)
    pivotmax=cc(2)
    pivotmin=-cc(3); if(pivotmin.eq.999d9)pivotmin=-1.d0

    ! share blocks to the LEFT
    q=own(me); p=own(me+1)-1
    needsend=upd(q).and.(me.gt.0); needrecv=upd(p+1).and.(me.lt.nproc-1)
    allocate(acol1(rr(q-1),n(q)),arow1(rr(p),n(p+1)),brow1(r(p),n(p+1)),brow(r(p),n(p+1),rr(p+1)), stat=info)
    if(info.ne.0)then;write(*,*)'fail to allocate for moving left';stop;endif
    if(needsend.and.needrecv)then
     !write(*,'(a,i2,2(a,i2,a,2i4))')'[',me,'] left: sending block ',q,':',rr(q-1),n(q),' receiving block ',p,':',rr(p),n(p+1)
     call dcopy(rr(q-1)*n(q),arg%u(q)%p(1,1,r(q)),1,acol1,1)
     call mpi_sendrecv(acol1,rr(q-1)*n(q),typed,me-1,lft, &
                       arow1,rr(p)*n(p+1),typed,me+1,lft, &
                       MPI_COMM_WORLD,stat,info)
     if(info.ne.0)then;write(*,*)'mpi: sendrecv going left fail: ',info,me;stop;endif
    else if(needsend)then
     !write(*,'(a,i2,a,i2,a,2i4)')'[',me,']: sending block ',q,' to the left: ',rr(q-1),n(q)
     call dcopy(rr(q-1)*n(q),arg%u(q)%p(1,1,r(q)),1,acol1,1)
     call mpi_send(acol1,rr(q-1)*n(q),typed,me-1,lft,MPI_COMM_WORLD,info)
     if(info.ne.0)then;write(*,*)'mpi: send going left fail: ',info,me;stop;endif
    else if(needrecv)then
     !write(*,'(a,i2,a,i2,a,2i4)')'[',me,']: receiving block ',p+1,' going left: ',rr(p),n(p+1)
     call mpi_recv(arow1,rr(p)*n(p+1),typed,me+1,lft,MPI_COMM_WORLD,stat,info)
     if(info.ne.0)then;write(*,*)'mpi: recv going left fail: ',info,me;stop;endif
    else
     continue
    end if
    if(needrecv)then
     call dcopy(r(p)*n(p+1)*rr(p+1),arg%u(p+1)%p,1,brow,1)
     deallocate(arg%u(p+1)%p); allocate(arg%u(p+1)%p(r(p),n(p+1),r(p+1)),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': fail to allocate memory to expand block';stop;endif
     call dcopy(r(p)*n(p+1)*rr(p+1),brow,1,arg%u(p+1)%p,1)
     forall(j=1:rr(p),k=1:n(p+1))arg%u(p+1)%p(j,k,r(p+1))=arow1(j,k)
     if(upd(p))then
      !write(*,'(a,i2,a,i2,a)') '[',me,']: expanding block ',p,' to match the right neighbour'
      ii=vip(p)%p(1,r(p)); jj=vip(p)%p(2,r(p))
!$OMP PARALLEL DO
      do k=1,n(p+1)
       arg%u(p+1)%p(r(p),k,r(p+1))=dmrgg_fun(ii,jj,k,r(p+1), p, l,m,vip,r,n, fun, par)
      end do
!$OMP END PARALLEL DO
      do k=1,n(p+1); amax=dmax1(amax,dabs(arg%u(p+1)%p(r(p),k,r(p+1)))); enddo
      nevalloc=nevalloc+n(p+1)
     end if ! upd(p)
     !write(*,'(a,i2,a,i2,a,3i4,a,3i4)')'[',me,']{',p,'}: expand rows ',r(p),n(p+1),rr(p+1),' -> ',r(p),n(p+1),r(p+1)
     call dcopy(r(p)*n(p+1)*rr(p+1),row%u(p+1)%p,1,brow,1)
     call dcopy(r(p)*n(p+1),arg%u(p+1)%p(1,1,r(p+1)),1,brow1,1)
     call d2_luar(n(p+1),r(p),inv(p)%p,brow1)
     deallocate(row%u(p+1)%p); allocate(row%u(p+1)%p(r(p),n(p+1),r(p+1)),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': not enough memory for new rows';stop;endif
     call dcopy(r(p)*n(p+1)*rr(p+1),brow,1,row%u(p+1)%p,1)
     call dcopy(r(p)*n(p+1),brow1,1,row%u(p+1)%p(1,1,r(p+1)),1)
    end if 
    deallocate(acol1,arow1,brow1,brow, stat=info)
    if(info.ne.0)then;write(*,*)'fail to deallocate after moving left';stop;endif
   
    ! share blocks to the RIGHT
    q=own(me+1)-1; p=own(me)
    needsend=upd(q).and.(me.lt.nproc-1); needrecv=upd(p-1).and.(me.gt.0)
    allocate(arow1(n(q+1),rr(q+1)),acol1(n(p),rr(p)),bcol1(n(p),r(p)),bcol(rr(p-1),n(p),r(p)), stat=info)
    if(info.ne.0)then;write(*,*)'fail to allocate for moving right';stop;endif
    if(needsend.and.needrecv)then
     !write(*,'(a,i2,2(a,i2,a,2i4))')'[',me,'] right: sending block ',q+1,':',n(q+1),rr(q+1),' receiving block ',p,':',n(p),rr(p)
     call dcopy(n(q+1)*rr(q+1),arg%u(q+1)%p(r(q),1,1),r(q),arow1,1)
     call mpi_sendrecv(arow1,n(q+1)*rr(q+1),typed,me+1,rgt, &
                       acol1,n(p)*rr(p),typed,me-1,rgt, &
                       MPI_COMM_WORLD,stat,info)
     if(info.ne.0)then;write(*,*)'mpi: sendrecv going right fail: ',info,me;stop;endif
    else if(needsend)then
     !write(*,'(a,i2,a,i2,a,2i4)')'[',me,']: sending block ',q+1,' to the right: ',n(q+1),rr(q+1)
     call dcopy(n(q+1)*rr(q+1),arg%u(q+1)%p(r(q),1,1),r(q),arow1,1)
     call mpi_send(arow1,n(q+1)*rr(q+1),typed,me+1,rgt,MPI_COMM_WORLD,info)
     if(info.ne.0)then;write(*,*)'mpi: send going right fail: ',info,me;stop;endif
    else if(needrecv)then
     !write(*,'(a,i2,a,i2,a,2i4)')'[',me,']: receiving block ',p,' going right: ',n(p),rr(p)
     call mpi_recv(acol1,n(p)*rr(p),typed,me-1,rgt,MPI_COMM_WORLD,stat,info)
     if(info.ne.0)then;write(*,*)'mpi: recv going right fail: ',info,me;stop;endif
    else
     continue
    end if 
    if(needrecv)then
     call dcopy(rr(p-1)*n(p)*r(p),arg%u(p)%p,1,bcol,1)
     deallocate(arg%u(p)%p); allocate(arg%u(p)%p(r(p-1),n(p),r(p)),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': fail to allocate memory to expand block';stop;endif
     forall(i=1:rr(p-1),j=1:n(p),k=1:r(p))arg%u(p)%p(i,j,k)=bcol(i,j,k)
     forall(            j=1:n(p),k=1:rr(p))arg%u(p)%p(r(p-1),j,k)=acol1(j,k)
     !do k=1,r(p); do j=1,n(p); do i=1,rr(p-1); arg%u(p)%p(i,j,k)=bcol(i,j,k); enddo;enddo;enddo
     !do k=1,rr(p); do j=1,n(p); arg%u(p)%p(rr(p-1)+1,j,k)=acol1(j,k); enddo;enddo
     if(upd(p))then ! r(p)=rr(p)+1
      !write(*,'(a,i2,a,i2,a)') '[',me,']: expanding block ',p,' to match the left neighbour'
      kk=vip(p)%p(3,r(p));qq=vip(p)%p(4,r(p))
!$OMP PARALLEL DO
      do j=1,n(p)
       arg%u(p)%p(r(p-1),j,r(p))=dmrgg_fun(r(p-1),j,kk,qq, p, l,m,vip,r,n, fun, par)
      end do
!$OMP END PARALLEL DO
      do j=1,n(p); amax=dmax1(amax,dabs(arg%u(p)%p(r(p-1),j,r(p)))); enddo
      nevalloc=nevalloc+n(p)
     end if ! upd(p)
     !write(*,'(a,i2,a,i2,a,3i4,a,3i4)') '[',me,']{',p,'}: expand cols ',rr(p-1),n(p),r(p),' -> ',r(p-1),n(p),r(p)
     forall(i=1:rr(p-1),j=1:n(p),k=1:r(p))bcol(i,j,k)=col%u(p)%p(i,j,k)
     forall(            j=1:n(p),k=1:r(p))bcol1(j,k)=arg%u(p)%p(r(p-1),j,k)
     call d2_lual(n(p),r(p),inv(p)%p,bcol1)
     deallocate(col%u(p)%p); allocate(col%u(p)%p(r(p-1),n(p),r(p)),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': not enough memory for new cols';stop;endif
     forall(i=1:rr(p-1),j=1:n(p),k=1:r(p))col%u(p)%p(i,j,k)=bcol(i,j,k)
     forall(            j=1:n(p),k=1:r(p))col%u(p)%p(r(p-1),j,k)=bcol1(j,k)
    end if
    deallocate(arow1,acol1,bcol1,bcol, stat=info)
    if(info.ne.0)then;write(*,*)'fail to deallocate after moving right';stop;endif
   end if ! (nproc>1)
   pivotmax_prev = pivotmax
   
   call mpi_reduce(nevalloc,nevalall,1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,info)
   if(info.ne.0)then;write(*,*)subnam,': mpi reduce neval failed: ',info;stop;endif

   ! REPORT current progress
   t2 = timef()
   if(me.eq.0)write(str,'(i3,a2,a,f5.1,a,2f6.3)') &
     it,sdir,' rank', erank(arg),&
     ' log10t,n ',dlog(t2-t1)/dlog(10.d0),dlog(dble(nevalall))/dlog(10.d0)
   
   if(present(quad))then
    ttqq%l=l;ttqq%m=m;ttqq%r(l-1:m)=1;ttqq%n(l:m)=1
    ttqq%r(own(me)-1:own(me+1))=r(own(me)-1:own(me+1))
    call alloc(ttqq)
    first=own(me); last=own(me+1)-1; if(me.eq.nproc-1)last=m
    do p=first,last
     do k=1,r(p)
      call dgemv('n',r(p-1),n(p),1.d0,arg%u(p)%p(1,1,k),r(p-1),quad%u(p)%p,1,0.d0,ttqq%u(p)%p(1,1,k),1)
     end do 
    end do
    call dtt_lua(ttqq,inv,own)
    val=dtt_quad(ttqq,mybonds=own)
    call dealloc(ttqq)
   
    ! print error or internal conv
    if(me.eq.0)then
     if(present(tru))then
      write(stmp,'(a,e8.3,a,'//sfmt//')')' err ',dabs(1.d0-val/tru),' val ',val
     else 
      write(stmp,'(a,e8.3,a,'//sfmt//')')' cnv ',dabs(1.d0-val/val_prev),' val ',val
     end if
     str=trim(str)//trim(stmp)
    end if
    val_prev=val
   end if !(quad)
   if(me.eq.0)write(*,'(a)')trim(str)

   ! check exit conditions
   if(present(maxrank))ready=ready.or.(it+1.ge.maxrank)
   if(present(accuracy))then
    if(pivotmax.le.accuracy*amax)then;strike=strike+1;else;strike=0;endif
    ready=ready.or.(strike.ge.3)
   end if 
  end do  !(loop over ranks)
  
  !write(*,'(a,i1,a,2f15.7)') '[',me,']: ',dnrm2(r(p-1)*n(p)*r(p),arg%u(p)%p,1),dnrm2(r(p)*n(p+1)*r(p+1),arg%u(p+1)%p,1)
  !write(*,'(a,i1,a,2f15.7)') '[',me,']: ',minval(arg%u(p)%p),minval(arg%u(p+1)%p)
  call mpi_barrier(MPI_COMM_WORLD,info)
  if(info.ne.0)then;write(*,*)subnam,': mpi barrier fail: ',info;stop;endif
  !write(*,'(a,i2,a,10i3)')'[',me,']: ',r(l-1:m)

  call dtt_lua(arg,inv,own)
  
  do p=l-1,m
   deallocate(inv(p)%p,vip(p)%p)
  end do
  call dealloc(col); call dealloc(row)

  if(present(neval))then 
   if(nproc.gt.1)then
    call mpi_allreduce(nevalloc,nevalall,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
    if(info.ne.0)then;write(*,*)subnam,': mpi allreduce neval failed: ',info;stop;endif
    neval=nevalall
   else
    neval=nevalloc
   endif 
  end if  
 end subroutine 
 
 
 double precision function dmrgg_fun(i,j,k,q, p, l,m,vip,r,n, fun, par) result(f)
  implicit none
  integer,intent(in) :: i,j,k,q, p, l,m
  integer,intent(in) :: r(l-1:m),n(l:m)
  type(pointi2),intent(in) :: vip(0:tt_size)
  double precision,external :: fun
  double precision,intent(in),optional :: par(*)
  integer :: t,s, ind(tt_size)
  t=i; do s=p-1,l,-1; ind(s-l+1)=vip(s)%p(2,t); t=vip(s)%p(1,t); enddo
  ind(p  )=j; ind(p+1)=k
  t=q; do s=p+1,m-1,+1; ind(s+1-l+1)=vip(s)%p(3,t); t=vip(s)%p(4,t); enddo
  f=fun(m,ind,n,par)
 end function 

 subroutine dtt_accchk(nlot,arg,einf,efro,ainf,afro,fun,par,pivot)
  implicit none
  include "mpif.h"
  integer,intent(in) :: nlot
  type(dtt),intent(in) :: arg
  double precision,intent(out) :: einf,efro,ainf,afro
  double precision,external :: fun
  double precision,intent(inout),optional :: par(*)
  integer,intent(out),optional :: pivot(tt_size)
  character(len=*),parameter :: subnam='dtt_accchk'
  integer :: ilot,i,l,m,ind(tt_size),piv(tt_size),me,nproc,info
  double precision :: aval,bval, einfall,efroall,ainfall,afroall
  ! check mpi vars
  call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_size fail: ',info;stop;endif
  call mpi_comm_rank(MPI_COMM_WORLD,me,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_rank fail: ',info;stop;endif
  
  l=arg%l;m=arg%m
  einf=0.d0;efro=0.d0;ainf=0.d0;afro=0.d0
  do ilot=1,nlot/nproc
   do i=l,m; ind(i-l+1)=irnd(arg%n(i)); enddo
   aval=fun(m-l+1,ind,arg%n(l),par)
   bval=dtt_ijk(arg,ind)
   if(einf.lt.dabs(aval-bval))then;einf=dabs(aval-bval);piv(1:m-l+1)=ind(1:m-l+1);endif
   efro=efro+(aval-bval)**2
   ainf=dmax1(ainf,aval)
   afro=afro+aval**2
  end do
  
  call mpi_allreduce(einf,einfall,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  if(info.ne.0)then;write(*,*)subnam,': mpi allreduce einf failed: ',info;stop;endif
  call mpi_allreduce(ainf,ainfall,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  if(info.ne.0)then;write(*,*)subnam,': mpi allreduce ainf failed: ',info;stop;endif
  call mpi_allreduce(efro,efroall,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  if(info.ne.0)then;write(*,*)subnam,': mpi allreduce efro failed: ',info;stop;endif
  call mpi_allreduce(afro,afroall,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  if(info.ne.0)then;write(*,*)subnam,': mpi allreduce afro failed: ',info;stop;endif
  
  einf=einfall
  ainf=ainfall
  efro=dsqrt(efroall)
  afro=dsqrt(afroall)
  if(present(pivot))pivot(1:m-l+1)=piv(1:m-l+1)
 end subroutine
 
 subroutine dtt_lua(arg,inv,own)
  implicit none
  include 'mpif.h'
  type(dtt),intent(in) :: arg
  type(pointd) :: inv(0:tt_size)
  integer,intent(in) :: own(0:)

  character(len=*),parameter :: subnam='dtt_lua'
  integer,parameter :: tag=4
  integer :: me,nproc,p,q,m, typed,stat(MPI_STATUS_SIZE),info,r(0:tt_size),n(1:tt_size)
  
  ! check how our double precision type compares to hard-coded MPI types
  select case(storage_size(1.d0))
   case(32);  typed=MPI_REAL4;
   case(64);  typed=MPI_REAL8;
   case(128); typed=MPI_REAL16;
   case default
    if(me.eq.0)write(*,'(a,a,i5,a)')subnam,': unknown storage_size(1.d0): ',storage_size(1.d0),' guessing 64'
    typed=MPI_REAL8;
  end select  
 
  ! check mpi vars
  call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_size fail: ',info;stop;endif
  call mpi_comm_rank(MPI_COMM_WORLD,me,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_rank fail: ',info;stop;endif

  r=arg%r; n=arg%n

  if(nproc.gt.1)then
   ! local SEND RECV to share the rightmost inv with the neighbour
   if(me.gt.0)then
    p = own(me)-1
    deallocate(inv(p)%p);allocate(inv(p)%p(r(p)**2), stat=info)
    if(info.ne.0)then;write(*,*)'fail to allocate inv to receive going right';stop;endif
   end if 
  
   p = own(me)-1; q = own(me+1)-1
   if(me.eq.0)then
    !write(*,'(a,i2,a,i2,a,i6)') '[',me,']: sending inv ',q,' to the right: ',r(q)**2
    call mpi_send(inv(q)%p,r(q)**2,typed,me+1,tag,MPI_COMM_WORLD,info)
    if(info.ne.0)then;write(*,*)'mpi: send inv to the right fail: ',info,me;stop;endif
   else if(me.eq.nproc-1)then
    !write(*,'(a,i2,a,i2,a,i6)') '[',me,']: receiving inv ',p,' going right: ',r(p)**2
    call mpi_recv(inv(p)%p,r(p)**2,typed,me-1,tag,MPI_COMM_WORLD,stat,info)
    if(info.ne.0)then;write(*,*)'mpi: recv inv going right fail: ',info,me;stop;endif
   else
    !write(*,'(a,i2,a,i2,a,i6,a,i2,a,i6)') '[',me,']: sending inv ',q,' of size ',r(q)**2,' recving ',p,' of size ',r(p)**2
    call mpi_sendrecv(inv(q)%p,r(q)**2,typed,me+1,tag, inv(p)%p,r(p)**2,typed,me-1,tag,  MPI_COMM_WORLD,stat,info)
    if(info.ne.0)then;write(*,*)'mpi: sendrecv inv going right fail: ',info,me;stop;endif
   end if
  end if
  
 !if(me.eq.0)write(*,*)'applying inverses'
  do p=own(me),own(me+1)-1
   !write(*,'(a,i2,a,i2,a)')'[',me,']{',p,'} applying inverses from both sides'
   call d2_luar(n(p)*r(p),r(p-1),inv(p-1)%p,arg%u(p)%p)
   call d2_lual(r(p-1)*n(p),r(p),inv(p)%p,arg%u(p)%p)
  end do
  if(me.eq.nproc-1)then
   m=own(me+1)
   !write(*,'(a,i2,a,i2,a)')'[',me,']{',m,'} applying inverses from the left side'
   call d2_luar(n(m)*r(m),r(m-1),inv(m-1)%p,arg%u(m)%p)
  end if
 end subroutine

 double precision function dtt_quad(arg,quad,mybonds) result(val) ! val is returned only on proc=0
  implicit none
  include 'mpif.h'
  type(dtt),intent(in),target :: arg
  type(dtt),intent(in),optional :: quad
  integer,intent(in),optional,target :: mybonds(0:)
  
  character(len=*),parameter :: subnam='dtt_quad'
  integer,parameter :: tagsize=11,tagdata=12
  integer,pointer :: r(:),n(:),own(:)
  integer :: me,her,him,nproc,info,stat(MPI_STATUS_SIZE),typed, mym,myn,herm,hern
  integer :: first,last, l,m,p,q,i,j,k, mn(2)
  double precision,pointer :: prev(:,:),curr(:,:),next(:,:)
  
  val=0.d0
  ! check how our double precision type compares to hard-coded MPI types
  select case(storage_size(1.d0))
   case(32);  typed=MPI_REAL4
   case(64);  typed=MPI_REAL8
   case(128); typed=MPI_REAL16
   case default
    if(me.eq.0)write(*,'(a,a,i5,a)')subnam,': unknown storage_size(1.d0): ',storage_size(1.d0),' guessing 64'
    typed=MPI_REAL8
  end select

  ! check mpi vars
  call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_size fail: ',info;stop;endif
  call mpi_comm_rank(MPI_COMM_WORLD,me,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_rank fail: ',info;stop;endif
  
  l=arg%l; m=arg%m; r=>arg%r; n=>arg%n
  ! distribute bonds (ranks) between procs
  if(present(mybonds))then
   own=>mybonds
  else 
   allocate(own(0:nproc), stat=info)
   if(info.ne.0)then;write(*,*)subnam,': cannot allocate own bonds';stop;endif
   call share(l,m-1,own)
  end if 
  !write(*,'(a,i3,a,32i4)')'[',me,']: own: ',own(0:nproc)
  
  !if(me.eq.0)write(*,*)'compute local products'
  first=own(me)
  last=own(me+1)-1; if(me.eq.nproc-1)last=m
  do p=first,last
   allocate( curr(r(p-1), r(p)) )
   if(present(quad))then
    do k=1,r(p)
     call dgemv('n',r(p-1),n(p),1.d0,arg%u(p)%p(1,1,k),r(p-1),quad%u(p)%p,1,0.d0,curr(1,k),1)
    end do
   else
    forall(i=1:r(p-1),k=1:r(p))curr(i,k)=sum(arg%u(p)%p(i,:,k))
   end if
   if(p.eq.first)then
    prev=>curr; nullify(curr)
   else 
    allocate(next(r(first-1),r(p)))
    call dgemm('n','n',r(first-1),r(p),r(p-1),1.d0,prev,r(first-1),curr,r(p-1),0.d0,next,r(first-1))
    deallocate(prev,curr)
    prev=>next; nullify(next)
   end if
   !write(*,'(a,i2,a,i2,a,2i5,a,2i5)')'[',me,']{',p,'} block shape ',shape(prev),' r(first-1) r(p): ',r(first-1),r(p)
  end do

  !write(*,'(a,i2,a)')'[',me,'] assembling along the binary tree'
  mym=r(first-1); myn=r(last)
  if(size(prev,1).ne.mym .or. size(prev,2).ne.myn)then
   write(*,'(a,i2,a,2i5,a,2i5)')'[',me,'] size mismatch: mym myn: ',mym,myn,' shape: ',shape(prev)
   stop
  endif
  q=1
  do while(q.lt.nproc)
   if(mod(me,2*q).eq.0)then
    her=me+q
    if(her.lt.nproc)then
     ! have a right neighbour - need to receive
     call mpi_recv(mn,2,MPI_INTEGER,her,tagsize,MPI_COMM_WORLD,stat,info)
     if(info.ne.0)then;write(*,*)subnam,': mpi recv for sizes fail:',info;stop;endif
     !write(*,'(a,i2,a,i2,a,8i5)')'[',me,'] received sizes from ',her,': ',mn
     herm=mn(1); hern=mn(2)
     if(myn.ne.herm)then
      write(*,'(a,i2,a,2i5,a,2i5)')'[',me,'] size mismatch: me ',mym,myn,' her ',herm,hern
      stop
     end if  
     allocate(curr(herm,hern), next(mym,hern),stat=info)
     if(info.ne.0)then;write(*,'(a,i2,a,4i10)')'[',me,'] cannot allocate block:',mym,myn,herm,hern;stop;endif
     !write(*,'(a,i2,a,i7,1x,i7,a,i2)')'[',me,']: receiving ',herm,hern,' from ',her
     call mpi_recv(curr,herm*hern,typed,her,tagdata,MPI_COMM_WORLD,stat,info)
     if(info.ne.0)then;write(*,*)subnam,': mpi recv for data fail:',info;stop;endif
     call dgemm('n','n',mym,hern,myn,1.d0,prev,mym,curr,herm,0.d0,next,mym)
     deallocate(prev,curr)
     prev=>next; nullify(next)
     myn=hern
    end if
   else
    if(mod(me,q).eq.0)then
     ! need to send
     him=me-q
     mn(1)=mym; mn(2)=myn
     !write(*,'(a,i2,a,i2,a,8i5)')'[',me,'] sending sizes to ',him,': ',mn
     call mpi_send(mn,2,MPI_INTEGER,him,tagsize,MPI_COMM_WORLD,info)
     if(info.ne.0)then;write(*,*)subnam,': mpi send for sizes fail:',info;stop;endif
     !write(*,'(a,i2,a,i7,1x,i7,a,i2)')'[',me,']: sending ',mym,myn,' to ',him
     call mpi_send(prev,mym*myn,typed,him,tagdata,MPI_COMM_WORLD,info)
     if(info.ne.0)then;write(*,*)subnam,': mpi send for data fail:',info;stop;endif
    end if
   end if
   ! walk up the tree
   q=q*2
  end do
  call mpi_barrier(mpi_comm_world,info)
  if(info.ne.0)then;write(*,*)subnam,': error at barrier: ',info;stop;endif
  if(me.eq.0)val=prev(1,1)
  deallocate(prev)
 end function
end module
