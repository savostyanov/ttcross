module ptype_lib
 implicit none
 
 type,public :: pointd3
  double precision,dimension(:,:,:),pointer :: p=>null()
 end type
 type,public :: pointd2
  double precision,dimension(:,:),pointer :: p=>null()
 end type
 type,public :: pointd
  double precision,dimension(:),pointer :: p=>null()
 end type
 
 type,public :: pointi
  integer,dimension(:),pointer :: p=>null()
 end type
 type,public :: pointi2
  integer,dimension(:,:),pointer :: p=>null()
 end type
 type,public :: pointi3
  integer,dimension(:,:,:),pointer :: p=>null()
 end type
 
 type,public :: pointz3
  double complex,dimension(:,:,:),pointer :: p=>null()
 end type
 type,public :: pointz2
  double complex,dimension(:,:),pointer :: p=>null()
 end type
 type,public :: pointz
  double complex,dimension(:),pointer :: p=>null()
 end type

end module
