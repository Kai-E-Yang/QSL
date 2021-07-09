module mod_param
implicit none
  save
  integer, parameter :: i4=SELECTED_INT_KIND(R=9) 
  integer, parameter :: i8=SELECTED_INT_KIND(R=18) 
  integer, parameter :: r8=SELECTED_REAL_KIND(P=15,R=300)
  character*100 :: par
  character*100 :: indataformat
  character*100 :: BfieldName, OutFileName

  integer :: dimx,dimy,dimz
  integer :: refine_dimx,refine_dimy,refine_dimz
  integer :: x_start,y_start,z_start
  integer :: x_end,y_end,z_end
  integer :: nlevel
  integer :: nthreads
  integer :: isn=10000
  integer :: iterMax

  real(kind=r8) :: delta_s
  real(kind=r8),dimension(:,:,:),allocatable :: Bx,By,Bz
  real,dimension(:),allocatable :: x,y,z
  real :: x_min,x_max,y_min,y_max,z_min,z_max 

  real(kind=r8),dimension(:,:,:,:),allocatable :: cal_data

  real(kind=r8) :: BoundaryEps
  real(kind=r8) :: eps_B

  namelist /filename_par/ BfieldName, OutFileName, indataformat
  namelist /cal_par/ nthreads,dimx,dimy,dimz,x_start,y_start, &
    z_start,x_end,y_end,z_end,nlevel,delta_s,x_min,x_max,y_min,y_max,z_min,z_max
    
end module mod_param
