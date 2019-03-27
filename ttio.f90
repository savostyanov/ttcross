module ttio_lib
 use tt_lib
 implicit none

 integer,parameter,private :: un=51, rec1=1
 character(len=*),parameter,private :: frm='unformatted', acc='stream'

 integer(kind=4),private  :: ver(2)=(/ 1, 0 /)

 type,private :: tthead
  sequence
  character(len=8) :: txt='TT      '
  integer(kind=4)  :: ver(2)=(/ 1, 0 /)
  integer(kind=4)  :: inf(4)=(/tt_size, 0, 0, 0/)
  character(len=64) :: comment
  integer(kind=4) :: i(8)
 end type

 interface read
  module procedure dtt_read,ztt_read
 end interface
 interface write
  module procedure dtt_write,ztt_write
 end interface

contains

! WRITE
 subroutine dtt_write(arg,fnam,info)
  implicit none
  type(dtt),intent(in) :: arg
  character(len=*),intent(in) :: fnam
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='dtt_write'
  type(tthead) :: head
  integer :: io,u,i,j,k,b,sz
  integer(kind=4) :: l,m,n(tt_size),r(0:tt_size)
  real(kind=8),allocatable :: x(:)
  logical :: ex,op

  if(present(info))info=-11
  inquire (file=fnam, exist=ex, opened=op)
  if(op)then
   write(*,*)subnam,': file is open, trying to close: ',fnam
   inquire(file=fnam, number=u)
   write(*,*)subnam,': establish unit: ',u
   close(unit=u,status='keep')
   write(*,*)subnam,': closed ok'
  end if

  u=un; op=.true.
  do  while(op)
   inquire(unit=u,opened=op)
   if(op)then
    write(*,*)subnam,': unit ',u,' is busy, trying next '
    u=u+1
   end if
  end do

  sz=memory(arg)
  if(sz.le.0)then;write(*,*)subnam,': tt structure has invalid size: ',sz;endif
  allocate(x(max(sz,1)),stat=i)
  if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif

  l=arg%l; m=arg%m; n=arg%n; r=arg%r; sz=0
  do b=l,m
   call dcopy(r(b-1)*n(b)*r(b),arg%u(b)%p,1,x(sz+1),1)
   sz=sz+r(b-1)*n(b)*r(b)
  end do

  open(unit=u,file=fnam,form=frm,access=acc,action='write',position='rewind',status='replace',err=101,iostat=io)

  head%i(1)=l
  head%i(2)=m
  write(u,err=111,iostat=io) head
  write(u,err=112,iostat=io) l,m
  write(u,err=113,iostat=io) n(l:m),r(l-1:m)
  write(u,err=114,iostat=io) x
  close(u,err=121,iostat=io)
  if(present(info))info=0
  deallocate(x,stat=i)
  if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif
  return

101 continue
  write(*,*) subnam,': error opening file: ',io
  if(present(info))info=io
  return
111 continue
  write(*,*) subnam,': error writing header: ',io
  if(present(info))info=io
  return
112 continue
  write(*,*) subnam,': error writing lm: ',io
  if(present(info))info=io
  return
113 continue
  write(*,*) subnam,': error writing nr: ',io
  if(present(info))info=io
  return
114 continue
  write(*,*) subnam,': error writing cores: ',io
  if(present(info))info=io
  return
121 continue
  write(*,*) subnam,': error closing file: ',io
  if(present(info))info=io
  return
 end subroutine
 subroutine ztt_write(arg,fnam,info)
  implicit none
  type(ztt),intent(in) :: arg
  character(len=*),intent(in) :: fnam
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='ztt_write'
  type(tthead) :: head
  integer :: io,u,i,j,k,b,sz
  integer(kind=4) :: l,m,n(tt_size),r(0:tt_size)
  complex(kind=8),allocatable :: x(:)
  logical :: ex,op

  if(present(info))info=-11
  inquire (file=fnam, exist=ex, opened=op)
  if(op)then
   write(*,*)subnam,': file is open, trying to close: ',fnam
   inquire(file=fnam, number=u)
   write(*,*)subnam,': establish unit: ',u
   close(unit=u,status='keep')
   write(*,*)subnam,': closed ok'
  end if

  u=un; op=.true.
  do  while(op)
   inquire(unit=u,opened=op)
   if(op)then
    write(*,*)subnam,': unit ',u,' is busy, trying next '
    u=u+1
   end if
  end do

  sz=memory(arg)
  if(sz.le.0)then;write(*,*)subnam,': tt structure has invalid size: ',sz;endif
  allocate(x(max(sz,1)),stat=i)
  if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif

  l=arg%l; m=arg%m; n=arg%n; r=arg%r; sz=0
  do b=l,m
   call zcopy(r(b-1)*n(b)*r(b),arg%u(b)%p,1,x(sz+1),1)
   sz=sz+r(b-1)*n(b)*r(b)
  end do

  open(unit=u,file=fnam,form=frm,access=acc,action='write',position='rewind',status='replace',err=101,iostat=io)

  head%i(1)=l
  head%i(2)=m
  head%inf(2)=1 ! Complex!
  write(u,err=111,iostat=io) head
  write(u,err=112,iostat=io) l,m
  write(u,err=113,iostat=io) n(l:m),r(l-1:m)
  write(u,err=114,iostat=io) x
  close(u,err=121,iostat=io)
  if(present(info))info=0
  deallocate(x,stat=i)
  if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif
  return


101 continue
  write(*,*) subnam,': error opening file: ',io
  if(present(info))info=io
  return
111 continue
  write(*,*) subnam,': error writing header: ',io
  if(present(info))info=io
  return
112 continue
  write(*,*) subnam,': error writing lm: ',io
  if(present(info))info=io
  return
113 continue
  write(*,*) subnam,': error writing nr: ',io
  if(present(info))info=io
  return
114 continue
  write(*,*) subnam,': error writing cores: ',io
  if(present(info))info=io
  return
121 continue
  write(*,*) subnam,': error closing file: ',io
  if(present(info))info=io
  return
 end subroutine


! READ
 subroutine dtt_read(arg,fnam,info)
  implicit none
  type(dtt),intent(inout) :: arg
  character(len=*),intent(in) :: fnam
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='dtt_read'
  type(tthead) :: head
  integer :: io,u,i,j,k,b,sz
  integer(kind=4) :: l,m,n(tt_size),r(0:tt_size)
  real(kind=8),allocatable :: x(:)
  logical :: ex,op

  if(present(info))info=-11
  inquire (file=fnam, exist=ex, opened=op)
  if(.not.ex)then
   write(*,*)subnam,': file not exist: ',fnam
   if(present(info))info=-1
   return
  endif
  if(op)then
   write(*,*)subnam,': file is open, trying to close: ',fnam
   inquire(file=fnam, number=u)
   write(*,*)subnam,': establish unit: ',u
   close(unit=u,status='keep')
   write(*,*)subnam,': closed ok'
  end if

  u=un; op=.true.
  do  while(op)
   inquire(unit=u,opened=op)
   if(op)then
    write(*,*)subnam,': unit ',u,' is busy, trying next '
    u=u+1
   end if
  end do

  open(unit=u,file=fnam,form=frm,access=acc,action='read',position='rewind',status='old',err=101,iostat=io)
  read(u,err=111,iostat=io) head

  if(head%txt(1:2).ne.'TT')then
   write(*,*)subnam,': not TT header in file: ',fnam
   if(present(info))info=-2
   return
  end if
  if(head%ver(1).ne.ver(1))then
   write(*,*)subnam,': not correct version of TT file: ',head%ver
   if(present(info))info=-3
   return
  end if

  read(u,err=112,iostat=io) l,m
  arg%l=l; arg%m=m
  if(l.lt.0)then
   write(*,*)subnam,': read strange l,m: ',l,m
  end if

  read(u,err=113,iostat=io) n(l:m),r(l-1:m)
  arg%n(l:m)=n(l:m);arg%r(l-1:m)=r(l-1:m); call alloc(arg)

  sz=memory(arg)
  if(sz.le.0)then;write(*,*)subnam,': tt structure has invalid size: ',sz;endif
  allocate(x(max(sz,1)),stat=i)
  if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif

  read(u,err=114,iostat=io) x
  sz=0
  do b=l,m
   call dcopy(r(b-1)*n(b)*r(b),x(sz+1),1,arg%u(b)%p,1)
   sz=sz+r(b-1)*n(b)*r(b)
  end do

  close(u,err=121,iostat=io)
  if(present(info))info=0
  deallocate(x,stat=i)
  if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif
  return

101 continue
  write(*,*) subnam,': error opening file: ',io
  if(present(info))info=io
  return
111 continue
  write(*,*) subnam,': error reading header: ',io
  if(present(info))info=io
  return
112 continue
  write(*,*) subnam,': error reading lm: ',io
  if(present(info))info=io
  return
113 continue
  write(*,*) subnam,': error reading nr: ',io
  if(present(info))info=io
  return
114 continue
  write(*,*) subnam,': error writing cores: ',io
  if(present(info))info=io
  return
121 continue
  write(*,*) subnam,': error closing file: ',io
  if(present(info))info=io
  return
 end subroutine
 subroutine ztt_read(arg,fnam,info)
  implicit none
  type(ztt),intent(inout) :: arg
  character(len=*),intent(in) :: fnam
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='ztt_read'
  type(tthead) :: head
  integer :: io,u,i,j,k,b,sz
  integer(kind=4) :: l,m,n(tt_size),r(0:tt_size)
  complex(kind=8),allocatable :: x(:)
  logical :: ex,op

  if(present(info))info=-11
  inquire (file=fnam, exist=ex, opened=op)
  if(.not.ex)then
   write(*,*)subnam,': file not exist: ',fnam
   if(present(info))info=-1
   return
  endif
  if(op)then
   write(*,*)subnam,': file is open, trying to close: ',fnam
   inquire(file=fnam, number=u)
   write(*,*)subnam,': establish unit: ',u
   close(unit=u,status='keep')
   write(*,*)subnam,': closed ok'
  end if

  u=un; op=.true.
  do  while(op)
   inquire(unit=u,opened=op)
   if(op)then
    write(*,*)subnam,': unit ',u,' is busy, trying next '
    u=u+1
   end if
  end do

  open(unit=u,file=fnam,form=frm,access=acc,action='read',position='rewind',status='old',err=101,iostat=io)
  read(u,err=111,iostat=io) head

  if(head%txt(1:2).ne.'TT')then
   write(*,*)subnam,': not TT header in file: ',fnam
   if(present(info))info=-2
   return
  end if
  if(head%ver(1).ne.ver(1))then
   write(*,*)subnam,': not correct version of TT file: ',head%ver
   if(present(info))info=-3
   return
  end if

  read(u,err=112,iostat=io) l,m
  arg%l=l; arg%m=m
  if(l.lt.0)then
   write(*,*)subnam,': read strange l,m: ',l,m
  end if

  read(u,err=113,iostat=io) n(l:m),r(l-1:m)
  arg%n(l:m)=n(l:m);arg%r(l-1:m)=r(l-1:m); call alloc(arg)

  sz=memory(arg)
  if(sz.le.0)then;write(*,*)subnam,': tt structure has invalid size: ',sz;endif
  allocate(x(max(sz,1)),stat=i)
  if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif

  read(u,err=114,iostat=io) x
  sz=0
  do b=l,m
   call zcopy(r(b-1)*n(b)*r(b),x(sz+1),1,arg%u(b)%p,1)
   sz=sz+r(b-1)*n(b)*r(b)
  end do

  close(u,err=121,iostat=io)
  if(present(info))info=0
  deallocate(x,stat=i)
  if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif
  return

101 continue
  write(*,*) subnam,': error opening file: ',io
  if(present(info))info=io
  return
111 continue
  write(*,*) subnam,': error reading header: ',io
  if(present(info))info=io
  return
112 continue
  write(*,*) subnam,': error reading lm: ',io
  if(present(info))info=io
  return
113 continue
  write(*,*) subnam,': error reading nr: ',io
  if(present(info))info=io
  return
114 continue
  write(*,*) subnam,': error writing cores: ',io
  if(present(info))info=io
  return
121 continue
  write(*,*) subnam,': error closing file: ',io
  if(present(info))info=io
  return
 end subroutine



end module
