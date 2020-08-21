module mod_solver
use mod_param
use mod_operator
use omp_lib
implicit none
contains
!-------- the main loop for calculation --------
  subroutine do_cal()
    integer(kind=selected_int_kind(8)) :: TID
    integer,dimension(:),allocatable :: mark
    integer :: sum0,flag(5),k
    integer,dimension(8) :: time_begin,time_end
    real(kind=r8) :: time_delta(8)
    
    !------ Set parallel threads ------
    if(nthreads == 0) then
      nthreads = OMP_GET_MAX_THREADS()
    end if
    
    call OMP_SET_NUM_THREADS(nthreads)
    allocate(mark(0:(nthreads-1)))
    
    mark = 0
    sum0 = 0
    do k=1,5
       flag(k) = floor(0.2*k*refine_dimx)
    end do
    
    write(*,'(A)')'| Start the computation !'
    call date_and_time(VALUES=time_begin)

    if (refine_dimz.eq.1) then
!$omp parallel private(TID, sum0, k)
!$omp do schedule(dynamic,1)

    do k=1,refine_dimx
       call cal_plane_z(k)
       TID = OMP_GET_THREAD_NUM()
       mark(TID) = mark(TID) + 1
       sum0 = sum(mark)
       if(sum0 .eq. flag(1) .or. sum0 .eq. flag(2) .or. sum0 &
         & .eq. flag(3) .or. sum0 .eq. flag(4) .or. sum0 .eq.  &
         & flag(5)) then
         call date_and_time(VALUES=time_end)
         time_delta = 1.0d0*(time_end - time_begin)
         write(*,'(A,f6.1,A,A,f8.2)')'| Calculation progress: '&
         & ,100.*sum(mark)/(refine_dimx*1.0),'%','; costs time (min): '&
         & ,(time_delta(7)/3600.0 + time_delta(6)/60. &
         & + time_delta(5) + time_delta(3)*24.0)*60
       end if
    
    end do
!$omp end do
!$omp end parallel
    end if

    if (refine_dimy.eq.1) then
!$omp parallel private(TID, sum0, k)
!$omp do schedule(dynamic,1)

    do k=1,refine_dimx
       call cal_plane_y(k)
       TID = OMP_GET_THREAD_NUM()
       mark(TID) = mark(TID) + 1
       sum0 = sum(mark)
       if(sum0 .eq. flag(1) .or. sum0 .eq. flag(2) .or. sum0 &
         & .eq. flag(3) .or. sum0 .eq. flag(4) .or. sum0 .eq.  &
         & flag(5)) then
         call date_and_time(VALUES=time_end)
         time_delta = 1.0d0*(time_end - time_begin)
         write(*,'(A,f6.1,A,A,f8.2)')'| Calculation progress: '&
         & ,100.*sum(mark)/(refine_dimx*1.0),'%','; costs time (min): '&
         & ,(time_delta(7)/3600.0 + time_delta(6)/60. &
         & + time_delta(5) + time_delta(3)*24.0)*60
       end if
    
    end do
!$omp end do
!$omp end parallel
    end if

    if (refine_dimx.eq.1) then
!$omp parallel private(TID, sum0, k)
!$omp do schedule(dynamic,1)

    do k=1,refine_dimy
       call cal_plane_x(k)
       TID = OMP_GET_THREAD_NUM()
       mark(TID) = mark(TID) + 1
       sum0 = sum(mark)
       if(sum0 .eq. flag(1) .or. sum0 .eq. flag(2) .or. sum0 &
         & .eq. flag(3) .or. sum0 .eq. flag(4) .or. sum0 .eq.  &
         & flag(5)) then
         call date_and_time(VALUES=time_end)
         time_delta = 1.0d0*(time_end - time_begin)
         write(*,'(A,f6.1,A,A,f8.2)')'| Calculation progress: '&
         & ,100.*sum(mark)/(refine_dimy*1.0),'%','; costs time (min): '&
         & ,(time_delta(7)/3600.0 + time_delta(6)/60. &
         & + time_delta(5) + time_delta(3)*24.0)*60
       end if
    
    end do
!$omp end do
!$omp end parallel
    end if

    if(refine_dimx.ne.1 .and. refine_dimy.ne.1 .and. refine_dimz.ne.1) then
!$omp parallel private(TID, sum0, k)
!$omp do schedule(dynamic,1)

    do k=1,refine_dimx
       call cal_plane(k)
       TID = OMP_GET_THREAD_NUM()
       mark(TID) = mark(TID) + 1
       sum0 = sum(mark)
       if(sum0 .eq. flag(1) .or. sum0 .eq. flag(2) .or. sum0 &
         & .eq. flag(3) .or. sum0 .eq. flag(4) .or. sum0 .eq.  &
         & flag(5)) then
         call date_and_time(VALUES=time_end)
         time_delta = 1.0d0*(time_end - time_begin)
         write(*,'(A,f6.1,A,A,f8.2)')'| Calculation progress: '&
         & ,100.*sum(mark)/(refine_dimx*1.0),'%','; costs time (min): '&
         & ,(time_delta(7)/3600.0 + time_delta(6)/60. &
         & + time_delta(5) + time_delta(3)*24.0)*60
       end if
    
    end do
!$omp end do
!$omp end parallel
    end if
    call date_and_time(VALUES=time_end)
    write(*,'(A)')'| Computation finished !'
  end subroutine do_cal
! ------- do parallel cal on each layer ------
  subroutine cal_plane(k)
    integer,intent(in)::k
    integer::i,j
    real(kind=r8) :: q_slice(refine_dimy,refine_dimz)
    real(kind=r8) :: end_slice(refine_dimy,refine_dimz)
    real(kind=r8) :: length_slice(refine_dimy,refine_dimz)
    
    real(kind=r8) :: xindex,yindex,zindex
    
    xindex = x_start+(k-1)/real(nlevel,kind=r8)
    
    do j=1,refine_dimz
    do i=1,refine_dimy
       yindex = y_start+(i-1)/real(nlevel,kind=r8)
       zindex = z_start+(j-1)/real(nlevel,kind=r8)
      call cal_point(xindex,yindex,zindex,q_slice(i,j),&
            &end_slice(i,j),length_slice(i,j))
    end do
    end do
    cal_data(k,:,:,1) = q_slice
    cal_data(k,:,:,2) = length_slice
    cal_data(k,:,:,3) = end_slice
  end subroutine cal_plane
! ------- do parallel cal on z layer ------
  subroutine cal_plane_x(k)
    integer,intent(in)::k
    integer::i,j
    real(kind=r8) :: q_slice(refine_dimz)
    real(kind=r8) :: end_slice(refine_dimz)
    real(kind=r8) :: length_slice(refine_dimz)
    
    real(kind=r8) :: xindex,yindex,zindex
    
    yindex = y_start+(k-1)/real(nlevel,kind=r8)
    xindex = x_start

    do i=1,refine_dimz
       zindex = z_start+(i-1)/real(nlevel,kind=r8)
      call cal_point(xindex,yindex,zindex,q_slice(i),&
            &end_slice(i),length_slice(i))
    end do
    cal_data(1,k,:,1) = q_slice
    cal_data(1,k,:,2) = length_slice
    cal_data(1,k,:,3) = end_slice
  end subroutine cal_plane_x
! ------- do parallel cal on z layer ------
  subroutine cal_plane_y(k)
    integer,intent(in)::k
    integer::i,j
    real(kind=r8) :: q_slice(refine_dimz)
    real(kind=r8) :: end_slice(refine_dimz)
    real(kind=r8) :: length_slice(refine_dimz)
    
    real(kind=r8) :: xindex,yindex,zindex
    
    xindex = x_start+(k-1)/real(nlevel,kind=r8)
    yindex = y_start

    do i=1,refine_dimz
       zindex = z_start+(i-1)/real(nlevel,kind=r8)
      call cal_point(xindex,yindex,zindex,q_slice(i),&
            &end_slice(i),length_slice(i))
    end do
    cal_data(k,1,:,1) = q_slice
    cal_data(k,1,:,2) = length_slice
    cal_data(k,1,:,3) = end_slice
  end subroutine cal_plane_y
! ------- do parallel cal on z layer ------
  subroutine cal_plane_z(k)
    integer,intent(in)::k
    integer::i,j
    real(kind=r8) :: q_slice(refine_dimy)
    real(kind=r8) :: end_slice(refine_dimy)
    real(kind=r8) :: length_slice(refine_dimy)
    
    real(kind=r8) :: xindex,yindex,zindex
    
    xindex = x_start+(k-1)/real(nlevel,kind=r8)
    zindex = z_start

    do i=1,refine_dimy
       yindex = y_start+(i-1)/real(nlevel,kind=r8)
      call cal_point(xindex,yindex,zindex,q_slice(i),&
            &end_slice(i),length_slice(i))
    end do
    cal_data(k,:,1,1) = q_slice
    cal_data(k,:,1,2) = length_slice
    cal_data(k,:,1,3) = end_slice
  end subroutine cal_plane_z
! ------ do calculation at each point ------
  subroutine cal_point(PosiX,PosiY,PosiZ,SquashingQ,end_mark,length)
    real(kind=r8),intent(in) :: PosiX
    real(kind=r8),intent(in) :: PosiY
    real(kind=r8),intent(in) :: PosiZ
    real(kind=r8),intent(out) :: SquashingQ
    real(kind=r8),intent(out) :: end_mark
    real(kind=r8) :: DetMatrix
    real(kind=r8) :: NormMatrix
    real(kind=r8) :: Uf(3)
    real(kind=r8) :: Ub(3)
    real(kind=r8) :: Vf(3)
    real(kind=r8) :: Vb(3)
    real(kind=r8) :: Uf_tmp(3)
    real(kind=r8) :: Ub_tmp(3)
    real(kind=r8) :: Vf_tmp(3)
    real(kind=r8) :: Vb_tmp(3)
    real(kind=r8) :: NormVectorF(3)
    real(kind=r8) :: NormVectorB(3)
    real(kind=r8) :: length
    ! linep is the present position
    ! linef is the forward end of the integral line
    ! linef is the backward end of the integral line
    ! linep,f,b contains the value in the order of 
    ! 1,2,3: x,y,z coordinate, 
    ! 4,5,6, 7,8,9: x,y,z components of U, V vactors,
    ! 10,11,12,13: xyz components of the tangent vector and the norm of B
    ! at the end of the line.
    real(kind=r8) :: LineP(13)
    real(kind=r8) :: LineF(13)
    real(kind=r8) :: LineB(13)
    
    LineP(1:3)=(/PosiX,PosiY,PosiZ/)
    call initializeUV(LineP)
    call fieldline(LineP,LineF,LineB,length)
    ! define the normal vectors at the two targed planes.
    
    NormVectorF = LineF(10:12)
    NormVectorB = LineB(10:12)
    
    ! ------------ get the mark flag for the two end ---------
    ! meaning of the value of end_mark
    ! 1: close field line
    ! 2: open field line, positive end roots in bottom
    ! 3: open field line, negative end roots in bottom
    ! 4: field line without end roots in the bottom boundary
    end_mark = real(0,kind=r8)
    if(abs(LineF(3)) .lt. BoundaryEps .and. abs(LineB(3)) &
      & .lt. BoundaryEps) then
      end_mark = real(1,kind=r8)
    end if
    
    if(abs(LineB(3)) .lt. BoundaryEps .and. abs(LineF(3)) &
      & .gt. BoundaryEps) then
      end_mark = real(2,kind=r8)
    end if
    
    if(abs(LineF(3)) .lt. BoundaryEps .and. abs(LineB(3)) &
      & .gt. BoundaryEps) then
      end_mark = real(3,kind=r8)
    end if
    
    if(abs(LineF(3)) .gt. BoundaryEps .and. abs(LineB(3)) &
      & .gt. BoundaryEps) then
      end_mark = real(4,kind=r8)
    end if
    
    ! get the determinant of the Jacobian matrix
    DetMatrix = LineP(13)**2/(LineF(13)*LineB(13))
    
    ! get the projected U and V
    Uf_tmp = LineF(4:6)
    Vf_tmp = LineF(7:9)
    
    Ub_tmp = LineB(4:6)
    Vb_tmp = LineB(7:9)
    
    Uf = Uf_tmp - dot_product(Uf_tmp,NormVectorF)*NormVectorF
    Vf = Vf_tmp - dot_product(Vf_tmp,NormVectorF)*NormVectorF
    
    Ub = Ub_tmp - dot_product(Ub_tmp,NormVectorB)*NormVectorB 
    Vb = Vb_tmp - dot_product(Vb_tmp,NormVectorB)*NormVectorB 
    
    ! get the norm of the Jacobian matrix
    NormMatrix = dot_product(Uf,Uf)*dot_product(Vb,Vb) &
             & + dot_product(Ub,Ub)*dot_product(Vf,Vf) &
             & - 2.0*dot_product(Uf,Vf)*dot_product(Ub,Vb)
    
    ! value the squashing factor Q
    SquashingQ = log10(NormMatrix) - log10(DetMatrix)
    
    ! if the value of log10(Q) is less than log10(2), it will be revalued as log10(2)     
    
    if (SquashingQ .lt. log10(2.0d0)) then
       SquashingQ = real(log10(2.0d0), kind=r8)
    end if
    if (LineP(12) .lt. 0.0) then
       SquashingQ = real(-1, kind=r8)*SquashingQ
    end if
  end subroutine cal_point
! -------- calculate the field line ---------
  subroutine fieldline(LineP,LineF,LineB,length)
    real(kind=r8)::LineP(13)
    real(kind=r8)::LineF(13)
    real(kind=r8)::LineB(13)
    real(kind=r8)::sig,lengthF,lengthB,length
    
    sig = real(1,kind=r8)
    call integralLine(LineP,LineF,sig,lengthF)    
    sig = real(-1,kind=r8)
    call integralLine(LineP,LineB,sig,lengthB)
    call diffLine(LineP(1:3),LineP(10:13))
    length = abs(lengthB) + abs(lengthF)
  end subroutine fieldline
! -------- do the line integral --------
  subroutine integralLine(LineP,LineI,sig,s_end)
    integer,parameter :: neqn = 9
    integer::n_step,flag
    real(kind=r8),intent(in)::sig
    real(kind=r8)::LineP(13)
    real(kind=r8)::LineI(13)
    real(kind=r8)::tmp(9)
    real(kind=r8)::dtmp(4)
    real(kind=r8)::s_start,s_end,ds
    
    s_start = real(0,kind=r8)
    s_end = real(0,kind=r8)
    
    tmp = LineP(1:9)
    ds = real(delta_s,kind=r8)*sig
    n_step = 0
    flag = 1
    
    do while (n_step .le. isn .and. flag .eq. 1)
       s_end = s_start + ds
       call rk4(rhs,neqn,s_start,s_end,tmp,flag)
    
       s_start = s_end
       n_step  = n_step+1
    end do
    
    LineI(1:9) = tmp
    call diffLine(tmp(1:3),dtmp)
    LineI(10:13) = dtmp

  end subroutine integralLine
end module mod_solver