module zero_lib
 implicit none
 
 interface zeros
  module procedure d1_zero,d2_zero,d3_zero
  module procedure z1_zero,z2_zero,z3_zero
 end interface 
 
contains
 
 subroutine d1_zero(arg)
  double precision,intent(out) :: arg(:)
  double precision,parameter :: zero=0.d0
  integer :: n,i
  n=size(arg)
!$OMP PARALLEL DO SHARED(arg)
  do i=1,n
   arg(i)=zero
  end do
!$OMP END PARALLEL DO
 end subroutine
 subroutine d2_zero(arg)
  double precision,intent(out) :: arg(:,:)
  double precision,parameter :: zero=0.d0
  integer :: m,n,i,j
  m=size(arg,1);n=size(arg,2)
!$OMP PARALLEL DO SHARED(arg) PRIVATE(i)
  do j=1,n
   do i=1,m
    arg(i,j)=zero
   end do 
  end do
!$OMP END PARALLEL DO
 end subroutine
 subroutine d3_zero(arg)
  double precision,intent(out) :: arg(:,:,:)
  double precision,parameter :: zero=0.d0
  integer :: m,n,p,i,j,k
  m=size(arg,1);n=size(arg,2);p=size(arg,3)
!$OMP PARALLEL DO SHARED(arg) PRIVATE(i,j)
  do k=1,p
   do j=1,n
    do i=1,m
     arg(i,j,k)=zero
    end do 
   end do
  end do
!$OMP END PARALLEL DO
 end subroutine

 subroutine z1_zero(arg)
  double complex,intent(out) :: arg(:)
  double complex,parameter :: zero=(0.d0,0.d0)
  integer :: n,i
  n=size(arg)
!$OMP PARALLEL DO SHARED(arg)
  do i=1,n
   arg(i)=zero
  end do
!$OMP END PARALLEL DO
 end subroutine
 subroutine z2_zero(arg)
  double complex,intent(out) :: arg(:,:)
  double complex,parameter :: zero=(0.d0,0.d0)
  integer :: m,n,i,j
  m=size(arg,1);n=size(arg,2)
!$OMP PARALLEL DO SHARED(arg) PRIVATE(i)
  do j=1,n
   do i=1,m
    arg(i,j)=zero
   end do 
  end do
!$OMP END PARALLEL DO
 end subroutine
 subroutine z3_zero(arg)
  double complex,intent(out) :: arg(:,:,:)
  double complex,parameter :: zero=(0.d0,0.d0)
  integer :: m,n,p,i,j,k
  m=size(arg,1);n=size(arg,2);p=size(arg,3)
!$OMP PARALLEL DO SHARED(arg) PRIVATE(i,j)
  do k=1,p
   do j=1,n
    do i=1,m
     arg(i,j,k)=zero
    end do 
   end do
  end do
!$OMP END PARALLEL DO
 end subroutine

end module
