program main
 use tt_lib
 use ttio_lib
 use dmrgg_lib
 use time_lib
 implicit none
 type(dtt) :: tt
 integer :: m,n,r,piv,info
 double precision :: t1,t2,tcrs, einf,efro,ainf,afro, mem
 integer,parameter :: NCHK=2**10
 
 m=2**5 ! dimension of the problem
 n=2**5 ! mode size of the array
 r=16   ! TT rank bounds of the approximation
 piv=1  ! pivoting `strength'. 
        ! piv=-1 means full pivoting
        ! piv=p, 1<=p<=3 means random pivoting on n*r^p entries

 write(*,'(a)') 'Hi, this is TT cross interpolation ...'
 write(*,'(3x,a,i10)') 'dimension: ',m
 write(*,'(3x,a,i10)') 'mode size: ',n
 write(*,'(3x,a,i10)') 'TT ranks : ',r
 write(*,'(3x,a,i10)') 'pivoting : ',piv

 t1=timef()
 tt%l=1;tt%m=m;tt%n=n;tt%r=1;call alloc(tt)
 call dtt_dmrgg(tt, r, pivpar=piv)
 t2=timef()
 tcrs=t2-t1
 
 write(*,'(a,e12.4,a)') '... completed in ',tcrs,' sec.'
 
 write(*,'(a,i10,a)') 'Accuracy check on random set of ',NCHK,' entries ...'
 
 t1=timef()
 call dtt_accchk(2**10,tt,einf,efro,ainf,afro)
 t2=timef()
 write(*,'(a,e12.4,a)') '... completed in ',t2-t1,' sec.'
 write(*,'(a,e12.4)') 'Relative accuracy estimate in infinity  norm: ',einf/ainf
 write(*,'(a,e12.4)') 'Relative accuracy estimate in Frobenius norm: ',efro/afro
 
 mem=dtt_mb(tt)
 if(mem.lt.2**10)then
  write(*,'(a)')'Write the interpolant to file tt.tt ...'
  call dtt_write(tt,'tt.tt',info)
  if(info.ne.0)then;write(*,*)'cannot save tt';stop;endif
  write(*,'(a,e12.4,a)')'... file saved (',mem,' MB).'
 end if
 call dealloc(tt)
 write(*,'(a)')'Good bye.'
end program
