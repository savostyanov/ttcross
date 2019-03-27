module time_lib
!use omp_lib
integer,private,save :: time_cold=-1, time_per=0
contains
 double precision function timef( )
 implicit none
 !include 'mpif.h'
 integer :: c,r,m
 real*4 :: t(2),time
 real*8 :: mpitime
 real*8,external :: mpi_wtime
 !call etime(t,time)
 !timef=1.d0*time
 !write(*,*)t,time,timef
! call system_clock(count=c,count_rate=r,count_max=m)
! if(c.lt.cold)then
!  time_per=time_per+1
!  write(*,*)'TIME_LIB: increase time_per to: ',time_per,' maxcount: ',m
! end if 
! timef=dble(c+time_per*m)/r
! time_cold=c
 !timef=1.d0*omp_get_wtime()

 !timef=0.d0
 timef=1.d0*mpi_wtime()
 end function
end module 
