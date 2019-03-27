module quad_lib
 implicit none
 double precision,parameter,private :: tpi=6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843722347638652235862102370961489247599254991347037715054497824558763660238982596673467248813132861720427898927904494743814043597218874055410784343525863535047693496369353388102640011362542905271216555715426855155792183472743574429368818024499068602930991707421015845593785178470840399122242580439217280688363196272595495426199210374144226999999967459560999021194634656321926371900489189106938166052850446165066893700705238623763420200062756775057731750664167628412343553382946071965069808575109374623191257277647075751875039155637155610643424536132260038557532223918184328403d0
contains
 subroutine quad_rinv1(n,q)
  implicit none
  ! 1/t = \sum w_p exp(-a_p t^2), t >= 1
  ! w_p = 2 h c_p / sqrt(pi) (1+exp(-s_p)); a_p = log^2(1+exp(s_p))
  ! c_p = cosh t_p, s_p = sinh t_p, t_p = h*(p-n/2-1), p=1,n
  ! [Hackbusch, Khoromskij, 2006 part 1, (5.3)]
  integer,intent(inout) :: n
  double precision :: q(2,n)
  character(len=*),parameter :: subnam='quad_rinv1'
  integer :: nq,i,m
  double precision :: s,c,es,loghuge,wq,aq,ht,t
  
  t=0.d0
  loghuge=dlog(huge(t))
  nq=(n-3)/2;m=1
  ht=dlog(tpi*nq)/nq
  q=0.d0
  do i=-nq,nq
   t=dble(i)*ht; s=dsinh(t); c=dcosh(t)
   if(s.lt.-loghuge)then
    !q(1,1)=q(1,1)+2*c*ht/dsqrt(tpi/2)
    !q(1,1)=0.d0
   ! write(*,*) i, ' small'
   else if(s.gt.loghuge)then
   ! write(*,*) i, ' big'
   else
    m=m+1
    q(1,m)=2*c*ht/(dsqrt(tpi/2)*(1.d0+dexp(-s)))
    q(2,m)=(dlog(1.d0+dexp(s)))**2
   end if 
   !write(*,'(2e25.15)') q(:,m)
  end do
  n=m
 end subroutine

 subroutine testquad_rinv(nq,q,a,b,n,err)
  implicit none
  integer,intent(in) :: nq
  double precision,intent(in) :: q(2,nq)
  double precision,intent(in) :: a,b
  integer,intent(in) :: n
  double precision,intent(out),optional :: err
  character(len=*),parameter :: subnam='testquad_rinv'
  double precision :: t,x,y,lt,lx,ly,val,err1,errm
  integer :: i,j
  if(a.le.0 .or. b.le.0)then; write(*,*)subnam,': illegal interval: ',a,b;stop;endif
  if(b.lt.1)then;write(*,*)subnam,': too few step numbers: ',n;stop;endif
  x=dmin1(a,b);y=dmax1(a,b);lx=dlog(x);ly=dlog(y); errm=0.d0
  open(10,file='_testquad.rinv',access='sequential')
  do i=0,n-1
   lt=lx+dble(i)*(ly-lx)/(n-1)
   t=dexp(lt)
   val=0.d0
   do j=1,nq
    val=val+q(1,j)*dexp(-q(2,j)*t*t)
   end do
   err1=t*dabs(1.d0/t-val)
   errm=dmax1(err1,errm)
   write(10,'(e12.5,3x,2f25.15,3x,e10.2)')t,1.d0/t,val,err1
  end do
  close(10)
  if(present(err))err=errm
 end subroutine


 subroutine lgwt(n, x, w)
  implicit none
  integer, intent(in) :: n
  double precision, intent(out) :: x(n), w(n)
  !double precision,parameter :: small=3.e-15  !dble
  !double precision,parameter :: small=1.e-33   !quad
  double precision :: small ! adaptive
  integer :: i, j, m
  double precision::  p1, p2, p3, pp, z, z1
  small=5*epsilon(1.d0)

  m = (n + 1) / 2
  do i = 1,m
    z = dcos( (tpi * (4*i-1)) / (8*n+4) )
100 p1 = 1.0d0
    p2 = 0.0d0
    do j = 1, n
    p3 = p2
    p2 = p1
    p1 = ((2*j-1) * z * p2 - (j-1)*p3) / j
    enddo

    pp = n*(z*p1-p2)/(z*z-1)
    z1 = z
    z = z1 - p1/pp             ! newton's method

    if (dabs(z-z1) .gt. small) goto  100

      x(i)     =  -z
      x(n+1-i) =   z
      w(i)     = 2.d0/((1-z*z)*pp*pp)
      w(n+1-i) = w(i)
    end do     ! i loop
 end subroutine

end module
