program main
 use tt_lib
 use ttio_lib
 use dmrgg_lib
 use time_lib
 implicit none
 include "mpif.h"
 type(dtt) :: tt
 integer :: m,n,r,piv,info,nproc,me
 double precision :: t1,t2,tcrs, einf,efro,ainf,afro, mem
 integer,parameter :: NCHK=2**20

 call mpi_init(info)
 if(info.ne.0)then;write(*,*)'mpi: init fail: ',info;stop;endif

 call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
 if(info.ne.0)then;write(*,*)'mpi: comm_size fail: ',info;stop;endif
 call mpi_comm_rank(MPI_COMM_WORLD,me,info)
 if(info.ne.0)then;write(*,*)'mpi: comm_rank fail: ',info;stop;endif

 write(*,'(a,i3,a,i3)')'mpi: I am ',me,' of ',nproc
 
 m=nproc+1 ! dimension of the problem
 n=2**7 ! mode size of the array
 r=16   ! TT rank bounds of the approximation
 piv=1  ! pivoting `strength'. 
        ! piv=-1 means full pivoting
        ! piv=p, 1<=p<=3 means random pivoting on n*r^p entries
 
 if(me.eq.0)then
  write(*,'(a)') 'Hi, this is TT cross interpolation ...'
  write(*,'(3x,a,i10)') 'dimension: ',m
  write(*,'(3x,a,i10)') 'mode size: ',n
  write(*,'(3x,a,i10)') 'TT ranks : ',r
  write(*,'(3x,a,i10)') 'pivoting : ',piv
 end if

 t1=timef()
 tt%l=1;tt%m=m;tt%n=n;tt%r=1;call alloc(tt)
 call dtt_dmrgg(tt, r, pivpar=piv)
 t2=timef()
 tcrs=t2-t1
 
 if(me.eq.0)then
  write(*,'(a,e12.4,a)') '... completed in ',tcrs,' sec.'
  write(*,'(a,i10,a)') 'Accuracy check on random set of ',NCHK,' entries ...'
 endif
 
 t1=timef()
 call dtt_accchk(NCHK,tt,einf,efro,ainf,afro)
 t2=timef()
 if(me.eq.0)then
  write(*,'(a,e12.4,a)') '... completed in ',t2-t1,' sec.'
  write(*,'(a,e12.4)') 'Relative accuracy estimate in infinity  norm: ',einf/ainf
  write(*,'(a,e12.4)') 'Relative accuracy estimate in Frobenius norm: ',efro/afro
 endif

! mem=dtt_mb(tt)
! if(mem.lt.2**10)then
!  write(*,'(a)')'Write the interpolant to file tt.tt ...'
!  call dtt_write(tt,'tt.tt',info)
!  if(info.ne.0)then;write(*,*)'cannot save tt';stop;endif
!  write(*,'(a,e12.4,a)')'... file saved (',mem,' MB).'
! end if
! call dealloc(tt)
 if(me.eq.0)then
  write(*,'(a)')'Good bye.'
 endif
 call mpi_finalize(info)
 if(info.ne.0)then;write(*,*)'mpi: finalize fail: ',info;stop;endif
end program
