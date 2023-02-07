module mod_param
implicit none
  save
  integer, parameter :: i4=SELECTED_INT_KIND(R=9) 
  integer, parameter :: i8=SELECTED_INT_KIND(R=18) 
  integer, parameter :: r8=SELECTED_REAL_KIND(P=15,R=300)
  real,    parameter :: Pi=real(3.1415926,kind=r8)

  character*100 :: par
  character*100 :: indataformat
  character*100 :: BfieldName, OutFileName

  integer :: dim_th,dim_ph,dim_ra
  integer :: refine_dim_th,refine_dim_ph,refine_dim_ra
  logical :: outputvtk

  real(kind=r8) :: th_start,ph_start,ra_start
  real(kind=r8) :: th_end,ph_end,ra_end
  real(kind=r8) :: min_th,min_ph,min_ra
  real(kind=r8) :: max_th,max_ph,max_ra
  real(kind=r8) :: d_th,d_ph,d_ra

  integer :: nlevel
  integer :: nthreads
  integer :: isn=10000
  integer :: iterMax

  real(kind=r8),dimension(:,:,:),allocatable :: Bth,Bph,Bra
  real(kind=r8),dimension(:,:,:),allocatable :: Bx,By,Bz
  real(kind=r8),dimension(:,:,:,:),allocatable :: cal_data

  real(kind=r8),dimension(:),allocatable :: th,ph,ra
  real(kind=r8),dimension(:),allocatable :: refine_th
  real(kind=r8),dimension(:),allocatable :: refine_ph
  real(kind=r8),dimension(:),allocatable :: refine_ra

  real(kind=r8) :: delta_s
  real(kind=r8) :: BoundaryEps
  real(kind=r8) :: eps_B

  namelist /filename_par/ BfieldName, OutFileName, indataformat, outputvtk
  namelist /cal_par/ nthreads,dim_ra,dim_th,dim_ph,th_start,ph_start, &
    ra_start,ph_end,ra_end,th_end,nlevel,delta_s
    
end module mod_param