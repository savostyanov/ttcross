module time_lib
use omp_lib
!include 'mpif.h'
integer,private,save :: time_cold=-1, time_per=0
contains
 double precision function timef( )
 integer :: c,r,m
 real :: t(2)
 timef=etime(t)
! call system_clock(count=c,count_rate=r,count_max=m)
! if(c.lt.cold)then
!  time_per=time_per+1
!  write(*,*)'TIME_LIB: increase time_per to: ',time_per,' maxcount: ',m
! end if 
! timef=dble(c+time_per*m)/r
! time_cold=c
 
 !timef=omp_get_wtime()

 !timef=0.d0
 !time=mpi_wtime()
 end function
end module 
