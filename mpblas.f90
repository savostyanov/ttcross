!module mpblas_lib
!use mpmodule

!contains


 type(mp_real) function mpnrm2(n,mpx,incx) result(f)
  use mpmodule
  implicit none
  integer,intent(in)    :: n,incx
  type(mp_real)         :: mpx(*)
  type(mp_real)         :: tmp
  integer :: i
  f="0.d0"
!$OMP PARALLEL PRIVATE(i,tmp)
  tmp="0.d0"
!$OMP DO
  do i=0,n-1
   f=f+mpx(i*incx+1)**2
  end do
!$OMP END DO
!$OMP CRITICAL
  f=f+tmp
!$OMP END CRITICAL
!$OMP END PARALLEL
  f=sqrt(f)
 end function

 type(mp_real) function mpdot(n,mpx,incx,mpy,incy) result(f)
  use mpmodule
  implicit none
  integer,intent(in)    :: n,incx,incy
  type(mp_real)         :: mpx(*),mpy(*)
  type(mp_real)         :: tmp
  integer :: i
  f="0.d0"
!$OMP PARALLEL PRIVATE(i,tmp)
 tmp="0.d0"
!$OMP DO
  do i=0,n-1
   tmp=tmp+mpx(i*incx+1)*mpy(i*incy+1)
  end do
!$OMP END DO
!$OMP CRITICAL
  f=f+tmp
!$OMP END CRITICAL
!$OMP END PARALLEL
 end function

 integer function impamax(n,mpx,incx) result(ii)
  use mpmodule
  implicit none
  integer,intent(in)    :: n,incx
  type(mp_real)         :: mpx(*)
  integer :: i,itmp
  type(mp_real) :: dmax,tmp
  ii = 0; if (n.lt.1 .or. incx.le.0) return
  ii = 1; if (n.eq.1) return
  ii = 1
  dmax = abs(mpx(1))
!$OMP PARALLEL PRIVATE(i,tmp,itmp)
  tmp=dmax
  itmp=ii
!$OMP DO
  do i = 1,n-1
   if(abs(mpx(i*incx+1)).gt.tmp)then
    itmp = i+1
    tmp = abs(mpx(i*incx+1))
   end if
  end do
!$OMP END DO
!$OMP CRITICAL
  if(tmp.gt.dmax)then
   dmax=tmp
   ii=itmp
  end if
!$OMP END CRITICAL
!$OMP END PARALLEL
 end function

 subroutine mpcopy(n,mpx,incx,mpy,incy)
  use mpmodule
  implicit none
  integer,intent(in)    :: n,incx,incy
  type(mp_real)         :: mpx(*),mpy(*)
  integer :: i
!$OMP PARALLEL DO PRIVATE(i)
  do i=0,n-1
   mpy(i*incy+1)=mpx(i*incx+1)
  end do
!$OMP END PARALLEL DO
 end subroutine

 subroutine mpscal(n,mpa,mpx,incx)
  use mpmodule
  implicit none
  integer,intent(in)    :: n,incx
  type(mp_real)         :: mpa,mpx(*)
  integer :: i
!$OMP PARALLEL DO PRIVATE(i)
  do i=0,n-1
   mpx(i*incx+1)=mpa*mpx(i*incx+1)
  end do
!$OMP END PARALLEL DO
 end subroutine

 subroutine mpgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
  use mpmodule
  implicit none
  character,intent(in)  :: trans
  integer,intent(in)    :: m,n,lda,incx,incy
  type(mp_real)         :: alpha,beta, a(lda,*),x(*),y(*), zero
  integer :: i,j
  zero = "0"
  if(trans.eq.'n' .or. trans.eq.'N')then
!$OMP PARALLEL PRIVATE(i,j)
   if (beta.eq.zero) then
!$OMP DO
    do i=1,m
     y((i-1)*incy+1)=zero
    end do
!$OMP END DO
   else
!$OMP DO
    do i=1,m
     y((i-1)*incy+1)=beta*y((i-1)*incy+1)
    end do
!$OMP END DO
   end if
!$OMP DO
   do i=1,m
    do j=1,n
     y((i-1)*incy+1)=y((i-1)*incy+1)+alpha*a(i,j)*x((j-1)*incx+1)
    end do
   end do
!$OMP END DO
!$OMP END PARALLEL

  else if(trans.eq.'t' .or. trans.eq.'T')then
!$OMP PARALLEL PRIVATE(i,j)
   if (beta.eq.zero) then
!$OMP DO
    do j=1,n
     y((j-1)*incy+1)=zero
    end do
!$OMP END DO
   else
!$OMP DO
    do j=1,n
     y((j-1)*incy+1)=beta*y((j-1)*incy+1)
    end do
!$OMP END DO
   end if
!$OMP DO
   do j=1,n
    do i=1,m
     y((j-1)*incy+1)=y((j-1)*incy+1)+alpha*a(i,j)*x((i-1)*incx+1)
    end do
   end do
!$OMP END DO
!$OMP END PARALLEL
  else
   write(*,*)'mpgemv: unknown trans: ',trans
   stop
  end if
 end subroutine

 subroutine mpgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  use mpmodule
  implicit none
  character,intent(in)  :: transa,transb
  integer,intent(in)    :: m,n,k,lda,ldb,ldc
  type(mp_real)         :: alpha,beta, a(lda,*),b(ldb,*),c(ldc,*), zero
  type(mp_real)         :: tmp
  integer :: i,j,l
  logical :: nota,notb
  if(.not.(transa.eq.'n'.or.transa.eq.'N'.or.transa.eq.'t'.or.transa.eq.'T'.or.transa.eq.'c'.or.transa.eq.'C'))then
   write(*,*)'mpgemm: unknown transa: ',transa; stop
  end if
  if(.not.(transb.eq.'n'.or.transb.eq.'N'.or.transb.eq.'t'.or.transb.eq.'T'.or.transb.eq.'c'.or.transb.eq.'C'))then
   write(*,*)'mpgemm: unknown transb: ',transb; stop
  end if
  nota=(transa.eq.'n' .or. transa.eq.'N')
  notb=(transb.eq.'n' .or. transb.eq.'N')

  zero = "0"

!$OMP PARALLEL PRIVATE(i,j,l,tmp)

  if (beta.eq.zero) then
!$OMP DO
    do i=1,m
     do j=1,n
      c(i,j)=zero
     end do
    end do
!$OMP END DO
  else
!$OMP DO
    do i=1,m
     do j=1,n
      c(i,j)=beta*c(i,j)
     end do
    end do
!$OMP END DO
  end if

  if(notb)then
   if(nota)then ! C = alpha*A*B + beta*C
!$OMP DO
    do i=1,m
     do j=1,n
      tmp=zero
      do l=1,k
       tmp=tmp+a(i,l)*b(l,j)
      end do
      c(i,j)=c(i,j)+alpha*tmp
     end do
    end do
!$OMP END DO
   else         ! C = alpha*A'*B + beta*C
!$OMP DO
    do i=1,m
     do j=1,n
      tmp=zero
      do l=1,k
       tmp=tmp+a(l,i)*b(l,j)
      end do
      c(i,j)=c(i,j)+alpha*tmp
     end do
    end do
!$OMP END DO
   end if
  else
   if(nota)then ! C = alpha*A*B' + beta*C
!$OMP DO
    do i=1,m
     do j=1,n
      tmp=zero
      do l=1,k
       tmp=tmp+a(i,l)*b(j,l)
      end do
      c(i,j)=c(i,j)+alpha*tmp
     end do
    end do
!$OMP END DO
   else         ! C = alpha*A'*B' + beta*C
!$OMP DO
    do i=1,m
     do j=1,n
      tmp=zero
      do l=1,k
       tmp=tmp+a(l,i)*b(l,j)
      end do
      c(i,j)=c(i,j)+alpha*tmp
     end do
    end do
!$OMP END DO
   end if
  endif
!$OMP END PARALLEL
 end subroutine
!end module
