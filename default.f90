module default_lib
 implicit none
 interface default
  module procedure default_i, default_d, default_z, default_a
 end interface
contains
 integer function default_i(def,opt) result(fun)
  integer,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
 double precision function default_d(def,opt) result(fun)
  double precision,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
 double complex function default_z(def,opt) result(fun)
  double complex,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
 character function default_a(def,opt) result(fun)
  character,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
end module
