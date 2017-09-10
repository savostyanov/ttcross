module dmrgg_lib
 use rnd_lib
 use ptype_lib
 use tt_lib
 use ttind_lib
 use lr_lib
 implicit none
 integer,parameter :: lft=1,rgt=2
contains
 subroutine dtt_dmrgg(arg, rnk, pivpar, par)
  ![tt] approximate [dmrgg_fun] using dmrg with greedy cross interpolation algorithm 
  implicit none
  include "mpif.h"
  type(dtt),intent(inout),target :: arg
  integer,intent(in) :: rnk
  integer,intent(in),optional :: pivpar  ! -1=full, 1=partial
  double precision,intent(in),optional :: par(*)
  character(len=*),parameter :: subnam='dtt_dmrgg'
  double precision,parameter :: small=1.d-12
  integer,pointer :: r(:),n(:)
  integer :: info,l,m,i,j,k,p,q,ii,jj,kk,qq,s,t,pp,it,dir,pos,ijkq,ij,jk,kq,nlot,ilot,pivoting
  integer :: me,nproc,stat(MPI_STATUS_SIZE)
  double precision :: aval,pivot, t1,t2,t3
  double precision,allocatable :: a(:,:,:,:),b(:)
  double precision,allocatable :: brow(:,:,:),bcol(:,:,:),bcol1(:,:),brow1(:,:),acol1(:,:),arow1(:,:),inv(:,:)
  integer,allocatable :: vip(:,:,:,:),newvip(:,:), lot(:,:)
  integer :: rr(0:tt_size)
  type(pointd3) :: u(1:tt_size),v(1:tt_size)
  character(len=2) :: sdir
  logical :: havecol,haverow
  logical,allocatable :: upd(:)
  integer,external :: idamax
  double precision,external :: dnrm2,ddot
  pivoting=default(-1,pivpar)
  
  if(arg%l.gt.arg%m)then;write(*,*)subnam,': l,m: ',arg%l,arg%m;return;endif
  if(arg%l.ne.1)then;write(*,*)subnam,': shift l to 1';arg%l=1;arg%m=arg%m-arg%l+1;endif
  l=arg%l; m=arg%m
  if(arg%r(l-1).ne.1)then;write(*,*)subnam,': left border rank should be 1';stop;endif
  if(arg%r(m).ne.1)then;write(*,*)subnam,': right border rank should be 1';stop;endif
  r=>arg%r; n=>arg%n
 
  ! check mpi vars
  call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_size fail: ',info;stop;endif
  call mpi_comm_rank(MPI_COMM_WORLD,me,info)
  if(info.ne.0)then;write(*,*)'mpi: comm_rank fail: ',info;stop;endif
  
  ! init
  allocate(vip(l-1:m,2,rnk,2),newvip(4,l-1:m),upd(l-1:m),inv(rnk**2,l-1:m),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate init';stop;endif
  vip=0; vip(:,:,1,:)=1
  arg%r(0:m-l+1)=1; call alloc(arg)
  
  do p=l,m 
   do j=1,n(p); arg%u(p)%p(1,j,1)=dmrgg_fun(1,j,1,1, p, l,m,rnk,vip,r,n, par); enddo
  end do
  do p=l,m-1; inv(1,p)=1.d0/arg%u(p)%p(1,1,1); enddo
  
  do p=l,m-1
   allocate(u(p)%p(1,n(p),1),stat=info)
   if(info.ne.0)then;write(*,*)subnam,': allocation fail for u ',p;stop;endif
   call dcopy(n(p),arg%u(p)%p,1,u(p)%p,1)
  enddo
  do p=l+1,m
   allocate(v(p)%p(1,n(p),1),stat=info)
   if(info.ne.0)then;write(*,*)subnam,': allocation fail for v ',p;stop;endif
   call dcopy(n(p),arg%u(p)%p,1,v(p)%p,1)
   call dscal(n(p),inv(1,p-1),v(p)%p,1)
  enddo

  ! iterate
  do it=1,rnk-1
    !dir=it/2; dir=it-2*dir; dir=2-dir
    !if(dir.eq.1)then;sdir='>>';else if(dir.eq.2)then;sdir='<<';else;sdir='??';endif
    !write(*,'(i3,a2)') it,sdir
    if(me.eq.0)write(*,'(a,i3)')'iter: ',it
    !do pp=l,m-1
    p=me+l ! one proc one (super)block
     !p=pp; if(dir.eq.2)p=m-pp+l-1
     !write(*,'(i3,a2,i2)') it,sdir,p
     
     allocate(acol1(r(p-1),n(p)),arow1(n(p+1),r(p+1)),bcol1(r(p-1),n(p)),brow1(n(p+1),r(p+1)), stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate acol/arow';stop;endif
     havecol=.false.; haverow=.false.
     
     select case(pivoting)
     case(-1)  ! full pivoting
      
      allocate(a(r(p-1),n(p),n(p+1),r(p+1)), b(r(p-1)*n(p)*n(p+1)*r(p+1)),stat=info)
      if(info.ne.0)then;write(*,*)subnam,': cannot allocate superblock';stop;endif
     
      ! compute long indices and form rhs
      do ijkq=1,r(p-1)*n(p)*n(p+1)*r(p+1)
       i=ijkq-1
       q=i/(r(p-1)*n(p)*n(p+1)); i=i-q*r(p-1)*n(p)*n(p+1)
       k=i/(r(p-1)*n(p)); i=i-k*r(p-1)*n(p)
       j=i/r(p-1); i=i-j*r(p-1)
       i=i+1;j=j+1;k=k+1;q=q+1
       a(i,j,k,q)=dmrgg_fun(i,j,k,q, p, l,m,rnk,vip,r,n, par)
      end do
      call dcopy(r(p-1)*n(p)*n(p+1)*r(p+1),a,1,b,1)

      !call say(b2)
      !write(*,'(a,4i4)') 'b size: ',r(p-1),n(p),n(p+1),r(p+1)
      !write(*,*)'u'; call say(u(p)%p)
      !write(*,*)'v'; call say(v(p+1)%p)
      
      ! compute residual and do full pivoting
      call dgemm('n','n',r(p-1)*n(p),n(p+1)*r(p+1),r(p),-1.d0,u(p)%p,r(p-1)*n(p),v(p+1)%p,r(p),1.d0,b,r(p-1)*n(p))
     
      ijkq=idamax(r(p-1)*n(p)*n(p+1)*r(p+1),b,1)-1
      q=ijkq/(r(p-1)*n(p)*n(p+1)); ijkq=ijkq-q*r(p-1)*n(p)*n(p+1)
      k=ijkq/(r(p-1)*n(p)); ijkq=ijkq-k*r(p-1)*n(p)
      j=ijkq/r(p-1); ijkq=ijkq-j*r(p-1)
      i=ijkq; ii=i+1;jj=j+1;kk=k+1;qq=q+1
      ijkq=ii+r(p-1)*(jj-1+n(p)*(kk-1+n(p+1)*(qq-1)))
      pivot=b(ijkq)
      
      ! copy column and row to containers
      forall(i=1:r(p-1),j=1:n(p))  acol1(i,j)=a(i,j,kk,qq)
      forall(k=1:n(p+1),q=1:r(p+1))arow1(k,q)=a(ii,jj,k,q)
      havecol=.true.; haverow=.true.

      ! deallocate superblocks
      deallocate(a,b)

     case(1,2,3) ! partial pivoting on O(r^a n) random elements, a=pivoting
      
      nlot=(r(p-1)*n(p)+n(p+1)*r(p+1))*(r(p)**(pivoting-1))
      allocate(lot(nlot,4),b(nlot), stat=info)
      if(info.ne.0)then;write(*,*)subnam,': cannot allocate lot';stop;endif
      
      ! distribute lottery weights for indices columns and rows 
      bcol1=1.d0; brow1=1.d0
      
      ! remove cross indices from lottery
      do s=1,r(p)
       i=vip(p,lft,s,1);j=vip(p,lft,s,2)
       k=vip(p,rgt,s,2);q=vip(p,rgt,s,1)
       bcol1(i,j)=0.d0; brow1(k,q)=0.d0
      end do
      
      ! get indices from lottery
      call lottery2(nlot,r(p-1)*n(p),n(p+1)*r(p+1),bcol1,brow1,lot)
      call dcopy(nlot,lot(1,2),1,lot(1,3),1)
      do ilot=1,nlot
       lot(ilot,2)=(lot(ilot,1)-1)/r(p-1)+1; lot(ilot,1)=lot(ilot,1)-(lot(ilot,2)-1)*r(p-1)
       lot(ilot,4)=(lot(ilot,3)-1)/n(p+1)+1; lot(ilot,3)=lot(ilot,3)-(lot(ilot,4)-1)*n(p+1)
      end do
      
      ! get matrix elements and residuals for these indices
      do ilot=1,nlot
       i=lot(ilot,1); j=lot(ilot,2); k=lot(ilot,3); q=lot(ilot,4)
       b(ilot)=dmrgg_fun(i,j,k,q, p, l,m,rnk,vip,r,n, par)
       !write(*,'(4i4,a,4i4,a,e10.3)')i,j,k,q,' : ',ind(l:m),' : ',b(ilot)
       b(ilot)=b(ilot)-ddot(r(p),u(p)%p(i,j,1),r(p-1)*n(p),v(p+1)%p(1,k,q),1)
       !write(*,'(4i4,a,4i4,a,e10.3)')i,j,k,q,' : ',ind(l:m),' : ',b(ilot)
      end do

      ! find maximum pivot
      ilot=idamax(nlot,b,1)
      ii=lot(ilot,1); jj=lot(ilot,2); kk=lot(ilot,3); qq=lot(ilot,4)
      pivot=b(ilot)
      !write(*,'(a,4i4,a,e10.3)') 'pre-pivot ',ii,jj,kk,qq,' : ',pivot

      ! compute row / col through pivot
      dir=1
      if(dir.eq.1)then
       haverow=.true.
       forall(k=1:n(p+1),q=1:r(p+1))arow1(k,q)=dmrgg_fun(ii,jj,k,q, p, l,m,rnk,vip,r,n, par)
       call dcopy(n(p+1)*r(p+1),arow1,1,brow1,1)
       call dgemv('t',r(p),n(p+1)*r(p+1),-1.d0,v(p+1)%p,r(p),u(p)%p(ii,jj,1),r(p-1)*n(p),1.d0,brow1,1)
       kq=idamax(n(p+1)*r(p+1),brow1,1)
       qq=(kq-1)/n(p+1)+1; kk=kq-(qq-1)*n(p+1)
       !if(dabs(pivot).gt.(1.d0+small)*dabs(brow1(kk,qq)))then;write(*,*)subnam,': pivots wrong (row): ',pivot,brow1(kk,qq);stop;endif
       pivot=brow1(kk,qq)
      else if(dir.eq.2)then
       havecol=.true.
       forall(i=1:r(p-1),j=1:n(p))acol1(i,j)=dmrgg_fun(i,j,kk,qq, p, l,m,rnk,vip,r,n, par)
       call dcopy(r(p-1)*n(p),acol1,1,bcol1,1)
       call dgemv('n',r(p-1)*n(p),r(p),-1.d0,u(p)%p,r(p-1)*n(p),v(p+1)%p(1,kk,qq),1,1.d0,bcol1,1)
       ij=idamax(r(p-1)*n(p),bcol1,1)
       jj=(ij-1)/r(p-1)+1; ii=ij-(jj-1)*r(p-1)
       !if(dabs(pivot).gt.(1.d0+small)*dabs(bcol1(ii,jj)))then;write(*,*)subnam,': pivots wrong (col): ',pivot,bcol1(ii,jj);stop;endif
       pivot=bcol1(ii,jj)
      else
       write(*,*)subnam,': illegal dir: ',dir;stop
      end if
      
      deallocate(lot,b)
     case default
      write(*,*) subnam,': unknown pivoting: ',pivoting; stop
     end select
     !write(*,'(a,4i4,a,e10.3)') '    pivot ',ii,jj,kk,qq,' : ',pivot
     
     newvip(:,p)=(/-1,-1,-1,-1/)
     upd(p)=.false.
     rr(l-1:m)=arg%r(l-1:m)

     if(dabs(pivot).gt.1.e-14)then ! ACCEPT
      ! put new pivot in sets
      vip(p,lft,r(p)+1,1)=ii;vip(p,lft,r(p)+1,2)=jj
      vip(p,rgt,r(p)+1,2)=kk;vip(p,rgt,r(p)+1,1)=qq
      newvip(:,p)=(/ ii,jj,kk,qq /)
      upd(p)=.true.

      ! allocate memory for blocks etc
      allocate(bcol(r(p-1),n(p),r(p)+1),brow(r(p)+1,n(p+1),r(p+1)), stat=info)
      if(info.ne.0)then;write(*,*)subnam,': not enough memory for brow bcol';stop;endif
      call dcopy(r(p-1)*n(p)*r(p),arg%u(p)%p,1,bcol,1)
      !forall(i=1:r(p),j=1:n(p+1),k=1:r(p+1))brow(i,j,k)=arg%u(p+1)%p(i,j,k)
      do k=1,r(p+1);do j=1,n(p+1);call dcopy(r(p),arg%u(p+1)%p(1,j,k),1,brow(1,j,k),1);enddo;enddo

      !write(*,*) '(col)'
      ! compute new column TODO: -r_k^2
      if(havecol)then
       forall(i=1:r(p-1),j=1:n(p))bcol(i,j,r(p)+1)=acol1(i,j)
      else 
       forall(i=1:r(p-1),j=1:n(p))bcol(i,j,r(p)+1)=dmrgg_fun(i,j,kk,qq, p, l,m,rnk,vip,r,n, par)
       call dcopy(r(p-1)*n(p),bcol(1,1,r(p)+1),1,acol1,1)
!      call dcopy(r(p-1)*n(p),arg%u(p)%p(1,1,r(p)+1),1,acol1,1)
      endif

      !write(*,*) '(row)'
      ! compute new row TODO: -r_k^2
      if(haverow)then
       forall(k=1:n(p+1),q=1:r(p+1))brow(r(p)+1,k,q)=arow1(k,q)
      else 
       forall(k=1:n(p+1),q=1:r(p+1))brow(r(p)+1,k,q)=dmrgg_fun(ii,jj,k,q, p, l,m,rnk,vip,r,n, par)
       call dcopy(n(p+1)*r(p+1),brow(r(p)+1,1,1),r(p)+1,arow1,1)
!      call dcopy(n(p+1)*r(p+1),arg%u(p+1)%p(r(p)+1,1,1),r(p)+1,arow1,1)
      endif
     
      !!write(*,*) '(mat)'
      !! compute new submatrix TODO: -r_k^2
      !do t=1,r(p)+1
      ! do s=1,r(p)+1
      !  i=vip(p,lft,s,1);j=vip(p,lft,s,2)
      !  k=vip(p,rgt,t,2);q=vip(p,rgt,t,1)
      !  mat(s+(t-1)*(r(p)+1),p)=dmrgg_fun(i,j,k,q, p, l,m,rnk,vip,r,n, par)
      ! end do
      !end do
      
      ! update the inverse submatrix
      !call d2_lug(r(p)+1,mat(1,p),inv(1,p))
      pos=r(p)**2
      inv(pos+1,p)=1.d0/pivot
      call dcopy(r(p),v(p+1)%p(1,kk,qq),1,inv(pos+2,p),1)
      call dcopy(r(p),u(p)%p(ii,jj,1),r(p-1)*n(p),inv(pos+2+r(p),p),1)

      !write(*,*) '(arg)'
      ! recover TT-format for blocks
      deallocate(arg%u(p)%p, arg%u(p+1)%p)
      allocate(arg%u(p)%p(r(p-1),n(p),r(p)+1),arg%u(p+1)%p(r(p)+1,n(p+1),r(p+1)),stat=info)
      if(info.ne.0)then;write(*,*)subnam,': not enough memory for new blocks';stop;endif
      call dcopy(r(p-1)*n(p)*(r(p)+1),bcol,1,arg%u(p)%p,1)
      call dcopy((r(p)+1)*n(p+1)*r(p+1),brow,1,arg%u(p+1)%p,1)
      
      !write(*,*) '(uv)'
      ! recompute factors u,v
      call dcopy(r(p-1)*n(p)*r(p),u(p)%p,1,bcol,1)
      !forall(i=1:r(p),j=1:n(p+1),k=1:r(p+1))brow(i,j,k)=v(p+1)%p(i,j,k)
      do k=1,r(p+1);do j=1,n(p+1); call dcopy(r(p),v(p+1)%p(1,j,k),1,brow(1,j,k),1);enddo;enddo
      call d2_luac(r(p-1)*n(p),r(p)+1,  inv(1,p),arg%u(p)%p,  bcol,from=r(p)+1)
      call d2_luar(r(p)+1,n(p+1)*r(p+1),inv(1,p),arg%u(p+1)%p,brow,from=r(p)+1)
      deallocate(u(p)%p, v(p+1)%p)
      allocate(u(p)%p(r(p-1),n(p),r(p)+1),v(p+1)%p(r(p)+1,n(p+1),r(p+1)),stat=info)
      if(info.ne.0)then;write(*,*)subnam,': not enough memory for new factors';stop;endif
      call dcopy(r(p-1)*n(p)*(r(p)+1),bcol,1,u(p)%p,1)
      call dcopy((r(p)+1)*n(p+1)*r(p+1),brow,1,v(p+1)%p,1)
     
!      if(p.gt.l)then
!       !write(*,*) '(v)'
!       call dcopy(r(p-1)*n(p)*r(p),v(p)%p,1,bcol,1)
!       !call d2_luar(r(p-1),n(p)*(r(p)+1),inv(1,p-1),arg%u(p)%p,bcol)
!       call dcopy(r(p-1)*n(p),arg%u(p)%p(1,1,r(p)+1),1,acol1,1)
!       call d2_luar(r(p-1),n(p),inv(1,p-1),acol1,bcol1)
!       call dcopy(r(p-1)*n(p),bcol1,1,bcol(1,1,r(p)+1),1)
!       deallocate(v(p)%p)
!       allocate(v(p)%p(r(p-1),n(p),r(p)+1), stat=info)
!       if(info.ne.0)then;write(*,*)subnam,': not enough memory for v factor';stop;endif
!       call dcopy(r(p-1)*n(p)*(r(p)+1),bcol,1,v(p)%p,1)
!      end if
!      if(p.lt.m-1)then
!       !write(*,*) '(u)'
!       !forall(i=1:r(p),j=1:n(p+1),k=1:r(p+1))brow(i,j,k)=u(p+1)%p(i,j,k)
!       do k=1,r(p+1);do j=1,n(p+1); call dcopy(r(p),u(p+1)%p(1,j,k),1,brow(1,j,k),1);enddo;enddo
!       !call d2_luac((r(p)+1)*n(p+1),r(p+1),inv(1,p+1),arg%u(p+1)%p,brow)
!       call dcopy(n(p+1)*r(p+1),arg%u(p+1)%p(r(p)+1,1,1),r(p)+1,arow1,1)
!       call d2_luac(n(p+1),r(p+1),inv(1,p+1),arow1,brow1)
!       call dcopy(n(p+1)*r(p+1),brow1,1,brow(r(p)+1,1,1),r(p)+1)
!       deallocate(u(p+1)%p)
!       allocate(u(p+1)%p(r(p)+1,n(p+1),r(p+1)),stat=info)
!       if(info.ne.0)then;write(*,*)subnam,': not enough memory for u factor';stop;endif
!       call dcopy((r(p)+1)*n(p+1)*r(p+1),brow,1,u(p+1)%p,1) 
!      end if
      
      ! increase the rank and move on
      arg%r(p)=arg%r(p)+1
      deallocate(bcol,brow)
     end if
     deallocate(brow1,bcol1)

    ! BROADCAST
    
    ! share r and vip data
    call mpi_allgather(arg%r(p),1,MPI_INTEGER,arg%r(l),1,MPI_INTEGER,MPI_COMM_WORLD,info)
    if(info.ne.0)then;write(*,*)'mpi: allgather r fail: ',info;stop;endif
    
    call mpi_allgather(upd(p),1,MPI_LOGICAL,upd(l),1,MPI_LOGICAL,MPI_COMM_WORLD,info)
    if(info.ne.0)then;write(*,*)'mpi: allgather upd fail: ',info;stop;endif
    
    call mpi_allgather(newvip(1,p),4,MPI_INTEGER,newvip(1,l),4,MPI_INTEGER,MPI_COMM_WORLD,info)
    if(info.ne.0)then;write(*,*)'mpi: allgather vip fail: ',info;stop;endif
    do pp=l,m-1
     if(upd(pp))then
      vip(pp,lft,r(pp),1)=newvip(1,pp);vip(pp,lft,r(pp),2)=newvip(2,pp)
      vip(pp,rgt,r(pp),2)=newvip(3,pp);vip(pp,rgt,r(pp),1)=newvip(4,pp)
     end if
    end do

    ! SEND RECV
    ! share arg with neighbours
    allocate(bcol(rr(p-1),n(p),r(p)),brow(r(p),n(p+1),rr(p+1)),bcol1(n(p),rr(p)),brow1(rr(p),n(p+1)),stat=info)
    if(info.ne.0)then;write(*,*)'cannot allocate bcol/brow';stop;endif

    if(p.gt.l)then
     if(upd(p))then
      call mpi_send(acol1,rr(p-1)*n(p),MPI_DOUBLE_PRECISION,me-1,lft,MPI_COMM_WORLD,info)
      if(info.ne.0)then;write(*,*)'mpi: send(1) fail: ',info,me;stop;endif
     end if
    end if
    if(p.lt.m-1)then
     if(upd(p+1))then
      call mpi_recv(brow1,rr(p)*n(p+1),MPI_DOUBLE_PRECISION,me+1,lft,MPI_COMM_WORLD,stat,info)
      if(info.ne.0)then;write(*,*)'mpi: recv(1) fail: ',info,me;stop;endif
      call dcopy(r(p)*n(p+1)*rr(p+1),arg%u(p+1)%p,1,brow,1)
      deallocate(arg%u(p+1)%p)
      allocate(arg%u(p+1)%p(r(p),n(p+1),r(p+1)),stat=info)
      if(info.ne.0)then;write(*,*)subnam,': not enough memory for new blocks';stop;endif
      call dcopy(r(p)*n(p+1)*rr(p+1),brow,1,arg%u(p+1)%p,1)
      do j=1,n(p+1)
       call dcopy(rr(p),brow1(1,j),1,arg%u(p+1)%p(1,j,r(p+1)),1)
      end do
      if(upd(p))then
       ii=vip(p,lft,r(p),1);jj=vip(p,lft,r(p),2)
       do k=1,n(p+1)
        arg%u(p+1)%p(r(p),k,r(p+1))=dmrgg_fun(ii,jj,k,r(p+1), p, l,m,rnk,vip,r,n, par)
       end do
      end if
      ! expand factor v (r(p),n(p+1),rr(p+1)) -> (r(p),n(p+1),r(p+1))
      deallocate(brow1)
      allocate(brow1(r(p),n(p+1)),stat=info)
      call dcopy(r(p)*n(p+1)*rr(p+1),v(p+1)%p,1,brow,1)
      call dcopy(r(p)*n(p+1),arg%u(p+1)%p(1,1,r(p+1)),1,brow1,1)
      deallocate(v(p+1)%p)
      allocate(v(p+1)%p(r(p),n(p+1),r(p+1)),stat=info)
      if(info.ne.0)then;write(*,*)subnam,': not enough memory for new v';stop;endif
      call dcopy(r(p)*n(p+1)*rr(p+1),brow,1,v(p+1)%p,1)
      call d2_luar(r(p),n(p+1),inv(1,p),brow1,brow)
      call dcopy(r(p)*n(p+1),brow,1,v(p+1)%p(1,1,r(p+1)),1)
     end if
    end if
    
    if(p.lt.m-1)then
     if(upd(p))then
      call mpi_send(arow1,n(p+1)*rr(p+1),MPI_DOUBLE_PRECISION,me+1,rgt,MPI_COMM_WORLD,info)
      if(info.ne.0)then;write(*,*)'mpi: send(2) fail: ',info,me;stop;endif
     end if
    end if
    if(p.gt.l)then
     if(upd(p-1))then
      call mpi_recv(bcol1,n(p)*rr(p),MPI_DOUBLE_PRECISION,me-1,rgt,MPI_COMM_WORLD,stat,info)
      if(info.ne.0)then;write(*,*)'mpi: recv(2) fail: ',info,me;stop;endif
      call dcopy(rr(p-1)*n(p)*r(p),arg%u(p)%p,1,bcol,1)
      deallocate(arg%u(p)%p)
      allocate(arg%u(p)%p(r(p-1),n(p),r(p)),stat=info)
      if(info.ne.0)then;write(*,*)subnam,': not enough memory for new blocks';stop;endif
      do j=1,n(p); do k=1,r(p)
       call dcopy(rr(p-1),bcol(1,j,k),1,arg%u(p)%p(1,j,k),1)
      end do; end do
      call dcopy(n(p)*rr(p),bcol1,1,arg%u(p)%p(r(p-1),1,1),r(p-1))
      if(upd(p))then
       kk=vip(p,rgt,r(p),2);qq=vip(p,rgt,r(p),1)
       do j=1,n(p)
        arg%u(p)%p(r(p-1),j,r(p))=dmrgg_fun(r(p-1),j,kk,qq, p, l,m,rnk,vip,r,n, par)
       end do
      end if
      ! expand factor u (rr(p-1),n(p),r(p)) -> (r(p-1),n(p),r(p))
      deallocate(bcol1)
      allocate(bcol1(n(p),r(p)),stat=info)
      call dcopy(rr(p-1)*n(p)*r(p),u(p)%p,1,bcol,1)
      call dcopy(n(p)*r(p),arg%u(p)%p(r(p-1),1,1),r(p-1),bcol1,1)
      deallocate(u(p)%p)
      allocate(u(p)%p(r(p-1),n(p),r(p)),stat=info)
      if(info.ne.0)then;write(*,*)subnam,': not enough memory for new u';stop;endif
      do j=1,n(p); do k=1,r(p)
       call dcopy(rr(p-1),bcol(1,j,k),1,u(p)%p(1,j,k),1)
      end do; end do
      call d2_luac(n(p),r(p),inv(1,p),bcol1,bcol)
      call dcopy(n(p)*r(p),bcol,1,u(p)%p(r(p-1),1,1),r(p-1))
     end if
    end if

    deallocate(acol1,arow1,bcol1,brow1,bcol,brow)

    ! check u and v are correct
    !allocate(bcol(r(p-1),n(p),r(p)),brow(r(p),n(p+1),r(p+1)), stat=info)
    !if(info.ne.0)then;write(*,*)subnam,': not enough memory for brow bcol';stop;endif
    !call d2_luac(r(p-1)*n(p),r(p),  inv(1,p),arg%u(p)%p,  bcol)
    !call d2_luar(r(p),n(p+1)*r(p+1),inv(1,p),arg%u(p+1)%p,brow)
    !call daxpy(r(p-1)*n(p)*r(p),-1.d0,u(p)%p,1,bcol,1)
    !call daxpy(r(p)*n(p+1)*r(p+1),-1.d0,v(p+1)%p,1,brow,1)
    !write(*,'(a,i1,a,2e10.3)')'[',me,']',dnrm2(r(p-1)*n(p)*r(p),bcol,1)/dnrm2(r(p-1)*n(p)*r(p),u(p)%p,1),dnrm2(r(p)*n(p+1)*r(p+1),brow,1)/dnrm2(r(p)*n(p+1)*r(p+1),v(p+1)%p,1)
    !deallocate(bcol,brow)
    
    !end do !(loop over my cores)
  end do  !(loop over ranks)
  
  !write(*,'(a,i1,a,2f15.7)') '[',me,']: ',dnrm2(r(p-1)*n(p)*r(p),arg%u(p)%p,1),dnrm2(r(p)*n(p+1)*r(p+1),arg%u(p+1)%p,1)
  !write(*,'(a,i1,a,2f15.7)') '[',me,']: ',minval(arg%u(p)%p),minval(arg%u(p+1)%p)
  call mpi_barrier(MPI_COMM_WORLD,info)
  if(info.ne.0)then;write(*,*)subnam,': mpi barrier fail: ',info;stop;endif
  !write(*,'(a,i2,a,10i3)')'[',me,']: ',r(l-1:m)
    
  ! ASSEMBLE
  ! share inv matrices
  !if(me.eq.0)write(*,*)'share inv matrices...'
  if(p.gt.l)then
   call mpi_send(inv(1,p),r(p)**2,MPI_DOUBLE_PRECISION,me-1,lft,MPI_COMM_WORLD,info)
   if(info.ne.0)then;write(*,*)'mpi: send(inv) fail: ',info,me;stop;endif
  end if
  if(p.lt.m-1)then
   call mpi_recv(inv(1,p+1),r(p+1)**2,MPI_DOUBLE_PRECISION,me+1,lft,MPI_COMM_WORLD,stat,info)
   if(info.ne.0)then;write(*,*)'mpi: recv(inv) fail: ',info,me;stop;endif
  end if
  
  ! compute approximation
  !if(me.eq.0)write(*,*)'compute approximation...'
  if(p.eq.l)call dcopy(r(l-1)*n(l)*r(l),u(l)%p,1,arg%u(l)%p,1)
  if(p.eq.m-1)call dcopy(r(m-1)*n(m)*r(m),v(m)%p,1,arg%u(m)%p,1)
  if(p.lt.m-1)call d2_luac(r(p)*n(p+1),r(p+1),inv(1,p+1),v(p+1)%p,arg%u(p+1)%p)
  
  ! reallocate all inactive cores
  !if(me.eq.0)write(*,*)'reallocate cores...'
  do pp=l,m-1; deallocate(u(pp)%p); enddo
  do pp=l+1,m; deallocate(v(pp)%p); enddo
  do pp=l,m
   if(.not.(pp.eq.p .or. pp.eq.p+1))then
    deallocate(arg%u(pp)%p)
    allocate(arg%u(pp)%p(r(pp-1),n(pp),r(pp)),stat=info)
    if(info.ne.0)then;write(*,*)'cannot reallocate cores';stop;endif
   end if
  end do
  
  ! share cores
  !if(me.eq.0)write(*,*)'broadcast cores...'
  call mpi_bcast(arg%u(l)%p,r(l-1)*n(l)*r(l),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
  if(info.ne.0)then;write(*,*)subnam,': mpi bcast',l,'failed:',info;stop;endif
  do pp=l+1,m
   call mpi_bcast(arg%u(pp)%p,r(pp-1)*n(pp)*r(pp),MPI_DOUBLE_PRECISION,pp-l-1,MPI_COMM_WORLD,info)
   if(info.ne.0)then;write(*,*)subnam,': mpi bcast', pp,'failed:',info;stop;endif
  end do

  deallocate(inv,vip)
 end subroutine 
 

 subroutine dtt_accchk(nlot,arg,einf,efro,ainf,afro,par,pivot)
  implicit none
  include "mpif.h"
  integer,intent(in) :: nlot
  type(dtt),intent(in) :: arg
  double precision,intent(out) :: einf,efro,ainf,afro
  double precision,intent(in),optional :: par(*)
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
   aval=dfunc(m-l+1,ind,arg%n(l),par)
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
  
  efro=dsqrt(efro)
  afro=dsqrt(afro)
  if(present(pivot))pivot(1:m-l+1)=piv(1:m-l+1)
 end subroutine
 
 pure double precision function dmrgg_fun(i,j,k,q, p, l,m,rnk,vip,r,n, par) result(f)
  implicit none
  integer,intent(in) :: i,j,k,q, p, l,m,rnk
  integer,intent(in) :: vip(l-1:m,2,rnk,2),r(l-1:m),n(l:m)
  double precision,intent(in),optional :: par(*)
  integer :: t,s
  integer :: ind(tt_size)
  t=i; do s=p-1,1,-1; ind(s-l+1)=vip(s,lft,t,2); t=vip(s,lft,t,1); enddo
  ind(p  )=j
  ind(p+1)=k
  t=q; do s=p+1,m-1,+1; ind(s+1-l+1)=vip(s,rgt,t,2); t=vip(s,rgt,t,1); enddo
  f=dfunc(m,ind,n,par)
 end function 
 
 pure double precision function dfunc(m,ind,n,par) result(f)
 ! function to compute tensor entry
  implicit none
  integer,intent(in) :: m
  integer,intent(in) :: ind(m),n(m)
  double precision,intent(in),optional :: par(*)
  double precision :: r
  integer :: i,j
  if(present(par))then
   j=ind(m)-1
   do i=m-1,1,-1; j=j*n(i)+ind(i)-1; enddo
   j=j+1
   f=par(j)
   return
  end if 
  r=0.d0
  do i=1,m
   r=r+ind(i)**2
  end do
  f=1.d0/dsqrt(r)
  !f=1.d0/(sum(ind(1:m))-m+1)
 end function
 
end module
