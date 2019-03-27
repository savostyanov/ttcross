program main
 use tt_lib
 use ttaux_lib
 use dmrgg_lib
 use time_lib
 use quad_lib
 use mat_lib
 implicit none
 include "mpif.h"
 type(dtt) :: tt,qq
 integer :: i,j,m,n,p,nx,r,piv,decay,a,info,nproc,me,ind(tt_size),nn(tt_size)
 integer(kind=8) :: neval
 double precision :: f,t1,t2,tcrs, einf,efro,ainf,afro, mem,logdet, val,tru
 double precision,allocatable::par(:),mat(:,:)
 integer,allocatable :: own(:)
 double precision,parameter :: bnd=20.d0
 double precision,parameter :: tpi=6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843722347638652235862102370961489247599254991347037715054497824558763660238982596673467248813132861720427898927904494743814043597218874055410784343525863535047693496369353388102640011362542905271216555715426855155792183472743574429368818024499068602930991707421015845593785178470840399122242580439217280688363196272595495426199210374144226999999967459560999021194634656321926371900489189106938166052850446165066893700705238623763420200062756775057731750664167628412343553382946071965069808575109374623191257277647075751875039155637155610643424536132260038557532223918184328403d0
 character(len=32) :: aa
 double precision,external :: dfunc_gauss_discr
 double precision,external :: dlamch

 ! Read params
 call getarg(1,aa); read(aa,'(i10)')m;   if(m.eq.0)m=8       ! dimension of the problem
 call getarg(2,aa); read(aa,'(i10)')n;   if(n.eq.0)n=25      ! stoch. mode size
 call getarg(3,aa); read(aa,'(i10)')r;   if(r.eq.0)r=20      ! max rank
 call getarg(4,aa); read(aa,'(i10)')piv;                     ! pivoting strategy

 call mpi_init(info)
 if(info.ne.0)then;write(*,*)'mpi: init fail: ',info;stop;endif
 call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
 if(info.ne.0)then;write(*,*)'mpi: comm_size fail: ',info;stop;endif
 call mpi_comm_rank(MPI_COMM_WORLD,me,info)
 if(info.ne.0)then;write(*,*)'mpi: comm_rank fail: ',info;stop;endif
 !write(*,'(a,i3,a,i3)')'mpi: I am ',me,' of ',nproc
 allocate(own(0:nproc))
 
 if(me.eq.0)then
  write(*,'(a)') 'Hi, this is TT cross interpolation integrating generalised Gaussian ...'
  write(*,'(3x,a,i10)') 'dimension:',m
  write(*,'(3x,a,i10)') 'quad.mode:',n
  write(*,'(3x,a,i10)') 'TT ranks :',r
  write(*,'(3x,a,i10)') 'pivoting :',piv
  write(*,'(3x,a,i10)') 'MPI procs:',nproc
!!$OMP PARALLEL
!  if (omp_get_thread_num().eq.0) then
!   write(*,'(3x,a,i10)')'OMP thrds:', omp_get_num_threads()
!  end if 
!!$OMP END PARALLEL
  write(*,'(3x,a,i10)') 'sizeof(d):',storage_size(1.d0)
  write(*,'(3x,a,e10.3)') 'epsilon  :',epsilon(1.d0)
 end if


 allocate(par(5*nx+2*n+16+m*m+2*m),mat(m,m), stat=info)
 if(info.ne.0)then;write(*,*)'cannot allocate par:',info;stop;endif
 par(1)=bnd
 a=10
 !call eye(mat)
 !call dscal(m*m,2.d0,mat,1)
 !call laplace(mat); call matinv(mat)
 forall(i=1:m,j=1:m)mat(i,j)=dexp(-dble(abs(i-j)))
 call dcopy(m*m,mat,1,par(a+1),1)
 call dpotrf('l',m,mat,m,info)
 if(info.ne.0)then;write(*,*)'dportf: info: ',info; stop; endif
 logdet=0.d0
 do i=1,m
  logdet=logdet+2*dlog(dabs(mat(i,i)))
 end do
 tru=dexp( dlog(tpi)*dble(m)/2 - logdet/2) 

 qq%l=1;qq%m=m;forall(i=1:m)qq%n(i)=n+2*i;qq%r=1;call alloc(qq);
 do i=1,m
  qq%u(i)%p=2*bnd/(qq%n(i)-1)
  qq%u(i)%p(1,1,1)=1*bnd/(qq%n(i)-1)
  qq%u(i)%p(1,qq%n(i),1)=1*bnd/(qq%n(i)-1)
 enddo
 t1=timef()
 tt%l=1;tt%m=m;forall(i=1:m)tt%n(i)=qq%n(i);tt%r=1;call alloc(tt)

 !distribute bonds (ranks) between procs
 own(0)=tt%l
 do p=1,nproc-1; own(p) = tt%l + int(dble(tt%m-tt%l)*dble(p)/nproc); enddo
 own(nproc) = tt%m
 !write(*,'(a,i3,a,8i5)')'[',me,']: own: ',own(0:nproc)
 
 call dtt_dmrgg(tt,dfunc_gauss_discr,par,maxrank=r,own=own,pivoting=piv,neval=neval,quad=qq,tru=tru)
 t2=timef()
 tcrs=t2-t1
 if(me.eq.0)write(*,'(a,i12,a,e12.4,a)') '...with',neval,' evaluations completed in ',tcrs,' sec.'

 val = dtt_quad(tt,qq,own)
 if(me.eq.0) then
  !write(*,*) dexp(logdet)
  write(*,'(a,e40.30)') 'computed value:', val
  write(*,'(a,e40.30)') 'analytic value:', tru
  write(*,'(a,f7.2)')   'correct digits:', -dlog(dabs(1.d0-val/tru))/dlog(10.d0)
  write(*,'(a)')'Good bye.'
 endif
 call dealloc(tt)
 call mpi_finalize(info)
 if(info.ne.0)then;write(*,*)'mpi: finalize fail: ',info;stop;endif
end program

double precision function dfunc_gauss_discr(m,ind,n,par) result(f)
 implicit none
 integer,intent(in) :: m
 integer,intent(in) :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 integer :: nodes,weights,matrix,i,a
 double precision :: x(m),ax(m)
 double precision :: arg,tol,term,bnd
 double precision,parameter :: e=2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011573834187930702154089149934884167509244761460668082264800168477411853742345442437107539077744992069551702761838606261331384583000752044933826560297606737113200709328709127443747047230696977209310141692836819025515108657463772111252389784425056953696770785449969967946864454905987931636889230098793127736178215424999229576351482208269895193668033182528869398496465105820939239829488793320362509443117301238197068416140397019837679320683282376464804295311802328782509819455815301756717361332069811250996181881593041690351598888519345807273866738589422879228499892086805825749279610484198444363463244968487560233624827041978623209002160990235304369941849146314093431738143640546253152096183690888707016768396424378140592714563549061303107208510383750510115747704171898610687396965521267154688957035035d0
 double precision,external :: ddot
 bnd=par(1)
 matrix=10
 forall(i=1:m)x(i)=bnd*dble(2*ind(i)-1-n(i))/(n(i)-1)
 call dgemv('n',m,m,1.d0,par(matrix+1),m,x,1,0.d0,ax,1)
 arg=ddot(m,x,1,ax,1)
 arg=-arg/2
 f=dexp(arg)
 
 ! apply weights
 !do i=1,m
 ! f=f*par(weights+ind(i))
 !end do
end function

