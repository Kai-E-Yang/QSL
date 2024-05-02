module mod_operator
use mod_param
implicit none
contains

  function cross(a, b)
    real(kind = r8),dimension(3) :: cross
    real(kind = r8),dimension(3),intent(in) :: a, b
    
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross
! ------- function check whether the point outside the computational domain
  function outside_boundary(posi)
    logical:: outside_boundary
    real(kind=r8),intent(in) :: posi(3)
    real(kind=r8)::posi_rtp(3)
    
    call car2polar(posi,posi_rtp)
    
    if(posi_rtp(1) .ge. min_ra .and. posi_rtp(1) .le. max_ra) then
      outside_boundary = .false.
    else
      outside_boundary = .true.
    end if

  end function outside_boundary
! ------ function give the 3 components of a vector in Cartesian from Spherical coordinate
  function rtp2xyz(posi, vec)
    real(kind=r8) :: rtp2xyz(3)
    real(kind=r8),intent(in) :: posi(3)
    real(kind=r8),intent(in) :: vec(3)
    real(kind=r8) :: et(3),ep(3),er(3)
    
    er = posi/sqrt(dot_product(posi,posi))
    ep = (/-1.0*posi(2),posi(1),real(0,kind=r8)/)
    ep = ep/sqrt(dot_product(ep,ep))
    et = cross(ep, er)
    
    rtp2xyz(1) = vec(1)*et(1) + vec(1)*ep(1) + vec(1)*er(1)
    rtp2xyz(2) = vec(2)*et(2) + vec(2)*ep(2) + vec(2)*er(2)
    rtp2xyz(3) = vec(3)*et(3) + vec(3)*ep(3) + vec(3)*er(3)

  end function rtp2xyz

  subroutine car2polar(vec,vec_rtp)
    real(kind=r8),intent(in)::vec(3)
    real(kind=r8),intent(out)::vec_rtp(3)
    ! input is the position in catisian coordinate x,y,z
    ! output is the spherical coordinate radius,theta,phi
    vec_rtp(1) = sqrt(dot_product(vec,vec))
    vec_rtp(2) = acos(vec(3)/vec_rtp(1))
    
    if(vec(1)**2 + vec(2)**2 .gt. eps_B ) then
      if(vec(2) .ge. 0) then
        vec_rtp(3) = acos(vec(1)/sqrt(vec(1)**2 + vec(2)**2))
      else
        vec_rtp(3) = 2*Pi - acos(vec(1)/sqrt(vec(1)**2 + vec(2)**2))
      end if
    else
      vec_rtp(3) = 0.0_r8
    end if
  end subroutine car2polar

  function vec_rtp2car(rtp,B_r,B_th,B_ph)
    real(kind=r8) :: vec_rtp2car(3)
    real(kind=r8),intent(in) :: rtp(3)
    real(kind=r8),intent(in) :: B_r,B_th,B_ph
    real(kind=r8) :: r_i,t_i,p_i
    
    r_i = rtp(1)
    t_i = rtp(2)
    p_i = rtp(3)
    
    vec_rtp2car(1) = B_r*sin(t_i)*cos(p_i) + B_th*cos(t_i)*cos(p_i) - B_ph*sin(p_i)
    vec_rtp2car(2) = B_r*sin(t_i)*sin(p_i) + B_th*cos(t_i)*sin(p_i) + B_ph*cos(p_i)
    vec_rtp2car(3) = B_r*cos(t_i)          - B_th*sin(t_i)

  end function vec_rtp2car
! ---------------------------------------------------------------
! the first dimension of xv indicates the x,y,z components of B
! the location in the second dimension of xv 
! 
!               g       g
!              /       /
!             /       /
!    g- - - -8- - - -7- - - -g
!           /|      /|   
!          / |     / |  
! g- - - -5- - - -6  |
!   g- - -|- 4 - -|- 3- - - -g
!         | /     | /
!         |/      |/     
! g- - - -1- - - -2- - - -g
!         |       |
!         |       |
!         |       |
!         g- - - -g
! 
! 
! ------ calculate the Bxyz on 8 vertex -----------
  ! subroutine corner (nnn,xv)
  !   integer,intent(in)::nnn(3)
  !   integer::n1,n2,n3
  !   integer::n1p,n2p,n3p
  !   real(kind = r8),intent(out)::xv(3,8)
  !   n1 = nnn(1)
  !   n2 = nnn(2)
  !   n3 = nnn(3)
    
  !   n1p = nnn(1)+1
  !   n2p = nnn(2)+1
  !   n3p = nnn(3)+1

  !   xv(1,1) = Bx(n1,n2,n3)
  !   xv(1,2) = Bx(n1p,n2,n3)
  !   xv(1,3) = Bx(n1p,n2p,n3)
  !   xv(1,4) = Bx(n1,n2p,n3)
  !   xv(1,5) = Bx(n1,n2,n3p)
  !   xv(1,6) = Bx(n1p,n2,n3p)
  !   xv(1,7) = Bx(n1p,n2p,n3p)
  !   xv(1,8) = Bx(n1,n2p,n3p)
    
  !   xv(2,1) = By(n1,n2,n3)
  !   xv(2,2) = By(n1p,n2,n3)
  !   xv(2,3) = By(n1p,n2p,n3)
  !   xv(2,4) = By(n1,n2p,n3)
  !   xv(2,5) = By(n1,n2,n3p)
  !   xv(2,6) = By(n1p,n2,n3p)
  !   xv(2,7) = By(n1p,n2p,n3p)
  !   xv(2,8) = By(n1,n2p,n3p)
    
  !   xv(3,1) = Bz(n1,n2,n3)
  !   xv(3,2) = Bz(n1p,n2,n3)
  !   xv(3,3) = Bz(n1p,n2p,n3)
  !   xv(3,4) = Bz(n1,n2p,n3)
  !   xv(3,5) = Bz(n1,n2,n3p)
  !   xv(3,6) = Bz(n1p,n2,n3p)
  !   xv(3,7) = Bz(n1p,n2p,n3p)
  !   xv(3,8) = Bz(n1,n2p,n3p)
  ! end subroutine corner
  subroutine corner(nnn,xv)
      integer, intent(in) :: nnn(3)
      real(kind = r8), intent(out) :: xv(3,8)
  
      integer :: n1p, n2p, n3p
  
      ! Incremented indices
      n1p = nnn(1) + 1
      n2p = nnn(2) + 1
      n3p = nnn(3) + 1
  
      ! Set values for xv based on indices
      xv(1,:) = (/Bx(nnn(1), nnn(2), nnn(3)), Bx(n1p, nnn(2), nnn(3)), &
                 Bx(n1p, n2p, nnn(3)), Bx(nnn(1), n2p, nnn(3)), &
                 Bx(nnn(1), nnn(2), n3p), Bx(n1p, nnn(2), n3p), &
                 Bx(n1p, n2p, n3p), Bx(nnn(1), n2p, n3p)/)
  
      xv(2,:) = (/By(nnn(1), nnn(2), nnn(3)), By(n1p, nnn(2), nnn(3)), &
                 By(n1p, n2p, nnn(3)), By(nnn(1), n2p, nnn(3)), &
                 By(nnn(1), nnn(2), n3p), By(n1p, nnn(2), n3p), &
                 By(n1p, n2p, n3p), By(nnn(1), n2p, n3p)/)
  
      xv(3,:) = (/Bz(nnn(1), nnn(2), nnn(3)), Bz(n1p, nnn(2), nnn(3)), &
                 Bz(n1p, n2p, nnn(3)), Bz(nnn(1), n2p, nnn(3)), &
                 Bz(nnn(1), nnn(2), n3p), Bz(n1p, nnn(2), n3p), &
                 Bz(n1p, n2p, n3p), Bz(nnn(1), n2p, n3p)/)
  end subroutine corner

  subroutine corner_B(nnn,xv)
    integer, intent(in) :: nnn(7,3)
    real(kind = r8), intent(out) :: xv(7,3,8)
    integer :: n1p(7), n2p(7), n3p(7),idx
  
   ! Incremented indices
    n1p = nnn(:,1) + 1
    n2p = nnn(:,2) + 1
    n3p = nnn(:,3) + 1
    do idx=1,7
      ! Set values for xv based on indices
      xv(idx,1,:) = (/Bx(nnn(idx,1), nnn(idx,2), nnn(idx,3)), Bx(n1p(idx), nnn(idx,2), nnn(idx,3)), &
                 Bx(n1p(idx), n2p(idx), nnn(idx,3)), Bx(nnn(idx,1), n2p(idx), nnn(idx,3)), &
                 Bx(nnn(idx,1), nnn(idx,2), n3p(idx)), Bx(n1p(idx), nnn(idx,2), n3p(idx)), &
                 Bx(n1p(idx), n2p(idx), n3p(idx)), Bx(nnn(idx,1), n2p(idx), n3p(idx))/)
      xv(idx,2,:) = (/By(nnn(idx,1), nnn(idx,2), nnn(idx,3)), By(n1p(idx), nnn(idx,2), nnn(idx,3)), &
                 By(n1p(idx), n2p(idx), nnn(idx,3)), By(nnn(idx,1), n2p(idx), nnn(idx,3)), &
                 By(nnn(idx,1), nnn(idx,2), n3p(idx)), By(n1p(idx), nnn(idx,2), n3p(idx)), &
                 By(n1p(idx), n2p(idx), n3p(idx)), By(nnn(idx,1), n2p(idx), n3p(idx))/)
      xv(idx,3,:) = (/Bz(nnn(idx,1), nnn(idx,2), nnn(idx,3)), Bz(n1p(idx), nnn(idx,2), nnn(idx,3)), &
                 Bz(n1p(idx), n2p(idx), nnn(idx,3)), Bz(nnn(idx,1), n2p(idx), nnn(idx,3)), &
                 Bz(nnn(idx,1), nnn(idx,2), n3p(idx)), Bz(n1p(idx), nnn(idx,2), n3p(idx)), &
                 Bz(n1p(idx), n2p(idx), n3p(idx)), Bz(nnn(idx,1), n2p(idx), n3p(idx))/)
    end do
  end subroutine corner_B

  ! subroutine corner_dB (nnn,xv)
  !   integer,intent(in)::nnn(3)
  !   integer::n1,n2,n3
  !   integer::n1p,n2p,n3p
  !   integer::idx
  !   real(kind = r8),intent(out)::xv(9,8)
  !   n1 = nnn(1)
  !   n2 = nnn(2)
  !   n3 = nnn(3)
    
  !   n1p = nnn(1)+1
  !   n2p = nnn(2)+1
  !   n3p = nnn(3)+1

  !   do idx=1,9
  !     xv(idx,1) = dBxyz(n1,n2,n3,idx)
  !     xv(idx,2) = dBxyz(n1p,n2,n3,idx)
  !     xv(idx,3) = dBxyz(n1p,n2p,n3,idx)
  !     xv(idx,4) = dBxyz(n1,n2p,n3,idx)
  !     xv(idx,5) = dBxyz(n1,n2,n3p,idx)
  !     xv(idx,6) = dBxyz(n1p,n2,n3p,idx)
  !     xv(idx,7) = dBxyz(n1p,n2p,n3p,idx)
  !     xv(idx,8) = dBxyz(n1,n2p,n3p,idx)
  !   end do
  ! end subroutine corner_dB

! --------- get the tri-linear interpolation ----------
    ! the location of each index 
    !     8- - - -7
    !    /|      /|   
    !   / |     / |  
    !  5- - - -6  |
    !  |  4 - -|- 3   
    !  | /     | /
    !  |/      |/     
    !  1- - - -2
    ! 
    ! o indicates the location of dxyz in the unit square
    ! (0,1)---------(1,1)
    !     |   |     |
    !     |   |     |
    !     |---o-----|
    !     |   |     |
    ! (0,0)---------(1,0)
    ! the left bottom one is dxyz1
    ! the right up one is dxyz2=1-dxyz
  subroutine xitp (interp,xv,dxyz1)
    real(kind = r8),intent(out)::interp(3)
    real(kind = r8)::dxyz2(3),weight(8)
    real(kind = r8)::xv(3,8)
    real(kind = r8),intent(in)::dxyz1(3)

    dxyz2 = real(1,kind=r8) - dxyz1
    
    weight(1) = dxyz2(1)*dxyz2(2)*dxyz2(3)
    weight(2) = dxyz1(1)*dxyz2(2)*dxyz2(3)
    weight(3) = dxyz1(1)*dxyz1(2)*dxyz2(3)
    weight(4) = dxyz2(1)*dxyz1(2)*dxyz2(3)
    
    weight(5) = dxyz2(1)*dxyz2(2)*dxyz1(3)
    weight(6) = dxyz1(1)*dxyz2(2)*dxyz1(3)
    weight(7) = dxyz1(1)*dxyz1(2)*dxyz1(3)
    weight(8) = dxyz2(1)*dxyz1(2)*dxyz1(3)
    
    interp(1) = dot_product(xv(1,:),weight)
    interp(2) = dot_product(xv(2,:),weight)
    interp(3) = dot_product(xv(3,:),weight)

  end subroutine xitp

  subroutine xitp_dB (interp,xv,dxyz1)
    real(kind = r8),intent(out)::interp(9)
    real(kind = r8)::dxyz2(3),weight(8)
    real(kind = r8)::xv(9,8)
    real(kind = r8),intent(in)::dxyz1(3)
    integer :: idx

    dxyz2 = 1.0_r8 - dxyz1
    
    weight(1) = dxyz2(1)*dxyz2(2)*dxyz2(3)
    weight(2) = dxyz1(1)*dxyz2(2)*dxyz2(3)
    weight(3) = dxyz1(1)*dxyz1(2)*dxyz2(3)
    weight(4) = dxyz2(1)*dxyz1(2)*dxyz2(3)
    
    weight(5) = dxyz2(1)*dxyz2(2)*dxyz1(3)
    weight(6) = dxyz1(1)*dxyz2(2)*dxyz1(3)
    weight(7) = dxyz1(1)*dxyz1(2)*dxyz1(3)
    weight(8) = dxyz2(1)*dxyz1(2)*dxyz1(3)
    
    do idx=1,9
      interp(idx) = dot_product(xv(idx,:),weight)
    enddo
  end subroutine xitp_dB

  subroutine xitp_B (interp,xv,dxyz1)
    real(kind = r8),intent(out)::interp(7,3)
    real(kind = r8)::dxyz2(7,3),weight(7,8)
    real(kind = r8)::xv(7,3,8)
    real(kind = r8),intent(in)::dxyz1(7,3)
    integer :: idx

    dxyz2 = 1.0_r8 - dxyz1

    do idx=1,7
      weight(idx,1) = dxyz2(idx,1)*dxyz2(idx,2)*dxyz2(idx,3)
      weight(idx,2) = dxyz1(idx,1)*dxyz2(idx,2)*dxyz2(idx,3)
      weight(idx,3) = dxyz1(idx,1)*dxyz1(idx,2)*dxyz2(idx,3)
      weight(idx,4) = dxyz2(idx,1)*dxyz1(idx,2)*dxyz2(idx,3)
      
      weight(idx,5) = dxyz2(idx,1)*dxyz2(idx,2)*dxyz1(idx,3)
      weight(idx,6) = dxyz1(idx,1)*dxyz2(idx,2)*dxyz1(idx,3)
      weight(idx,7) = dxyz1(idx,1)*dxyz1(idx,2)*dxyz1(idx,3)
      weight(idx,8) = dxyz2(idx,1)*dxyz1(idx,2)*dxyz1(idx,3)

      interp(idx,1) = dot_product(xv(idx,1,:),weight(idx,:))
      interp(idx,2) = dot_product(xv(idx,2,:),weight(idx,:))
      interp(idx,3) = dot_product(xv(idx,3,:),weight(idx,:))
    enddo
    
  end subroutine xitp_B

! ------ give the equation of ODE ------
  subroutine diffLine(Posi,Tangent)
    real(kind=r8),intent(in)::Posi(3)
    real(kind=r8),intent(out)::Tangent(4)
    real(kind=r8)::posi_rtp(3)
    real(kind=r8)::drtp(3)
    real(kind=r8)::xv(3,8)
    real(kind=r8)::bvec(3)
    real(kind=r8)::binter
    integer      ::nnn(3)

    call car2polar(Posi,posi_rtp)
    call find_index(posi_rtp,nnn,drtp)
    call corner(nnn,xv)
    call xitp(bvec,xv,drtp)
    binter = sqrt(dot_product(bvec,bvec))    
    if (binter .lt. eps_B) then
       Tangent = (/0.0_r8,0.0_r8,0.0_r8,0.0_r8/)
    else
       Tangent(1:3) = bvec/binter
       Tangent(4)   = binter
    end if
  end subroutine diffLine

  ! subroutine diffB(Posi,dBxOut,dByOut,dBzOut,Tangent)
  !   real(kind = r8),intent(in)::Posi(3)
  !   real(kind = r8),intent(out)::dBxOut(3),dByOut(3),dBzOut(3)
  !   real(kind=r8),intent(out)::Tangent(4)
  !   real(kind=r8)::posi_rtp(3)
  !   real(kind=r8)::drtp(3)
  !   real(kind=r8)::xv(3,8)
  !   real(kind=r8)::bvec(3)
  !   real(kind=r8)::xv_db(9,8)
  !   real(kind=r8)::dbvec(9)
  !   real(kind=r8)::binter
  !   integer      ::nnn(3)
    
  !   call car2polar(Posi,posi_rtp)
  !   call find_index(posi_rtp,nnn,drtp)

  !   call corner_dB(nnn,xv_db)
  !   call xitp_dB(dbvec,xv_db,drtp)
  !   dBxOut = dbvec(1:3)
  !   dByOut = dbvec(4:6)
  !   dBzOut = dbvec(7:9)

  !   call corner(nnn,xv)
  !   call xitp(bvec,xv,drtp)
  !   binter = sqrt(dot_product(bvec,bvec))    
  !   if (binter .lt. eps_B) then
  !      Tangent = (/0.0_r8,0.0_r8,0.0_r8,0.0_r8/)
  !   else
  !      Tangent(1:3) = bvec/binter
  !      Tangent(4)   = binter
  !   end if

  ! end subroutine diffB

  ! subroutine cal_diffB(Posi,dBxOut,dByOut,dBzOut)
  !   real(kind = r8)::Posi(3)
  !   real(kind = r8)::dBxOut(3),dByOut(3),dBzOut(3)
  !   real(kind = r8)::Tan1(4),Tan2(4),Tan3(4)
  !   real(kind = r8)::Tan4(4),Tan5(4),Tan6(4)
  !   ! real(kind = r8)::dx(3),dy(3),dz(3)
  !   ! real(kind = r8)::dxyz
    
  !   ! dxyz = 1e-4
    
  !   ! dx = (/dxyz,0.0_r8,0.0_r8/)
  !   ! dy = (/0.0_r8,dxyz,0.0_r8/)
  !   ! dz = (/0.0_r8,0.0_r8,dxyz/)
    
  !   call diffLine(Posi+dx,Tan1)
  !   call diffLine(Posi-dx,Tan2)
  !   call diffLine(Posi+dy,Tan3)
  !   call diffLine(Posi-dy,Tan4)
  !   call diffLine(Posi+dz,Tan5)
  !   call diffLine(Posi-dz,Tan6)
        
  !   dBxOut(1) = Tan1(1) - Tan2(1)
  !   dBxOut(2) = Tan3(1) - Tan4(1)
  !   dBxOut(3) = Tan5(1) - Tan6(1)
    
  !   dByOut(1) = Tan1(2) - Tan2(2)
  !   dByOut(2) = Tan3(2) - Tan4(2)
  !   dByOut(3) = Tan5(2) - Tan6(2)
    
  !   dBzOut(1) = Tan1(3) - Tan2(3)
  !   dBzOut(2) = Tan3(3) - Tan4(3)
  !   dBzOut(3) = Tan5(3) - Tan6(3)
    
  !   dBxOut = dBxOut/(dxyz*2.0_r8)
  !   dByOut = dByOut/(dxyz*2.0_r8)
  !   dBzOut = dBzOut/(dxyz*2.0_r8)
  ! end subroutine cal_diffB

  subroutine diffLine_A(nnn,drtp,Tangent)
    integer,intent(in)::nnn(3)
    real(kind=r8),intent(in)::drtp(3)
    real(kind=r8),intent(out)::Tangent(4)
    real(kind=r8)::xv(3,8)
    real(kind=r8)::bvec(3)
    real(kind=r8)::binter

    call corner(nnn,xv)
    call xitp(bvec,xv,drtp)
    binter = sqrt(dot_product(bvec,bvec))    
    if (binter .lt. eps_B) then
       Tangent = (/0.0_r8,0.0_r8,0.0_r8,0.0_r8/)
    else
       Tangent(1:3) = bvec/binter
       Tangent(4)   = binter
    end if
  end subroutine diffLine_A

  subroutine diffLine_B(nnn,drtp,Tangent)
    integer,intent(in)::nnn(7,3)
    real(kind=r8),intent(in)::drtp(7,3)
    real(kind=r8),intent(out)::Tangent(7,4)
    real(kind=r8)::xv(7,3,8)
    real(kind=r8)::bvec(7,3)
    real(kind=r8)::binter(7)
    integer :: idx
    call corner_B(nnn,xv)
    call xitp_B(bvec,xv,drtp)
    binter=sqrt(bvec(:,1)**2+bvec(:,2)**2+bvec(:,3)**2)
    do idx=1,7
      if (binter(idx) .lt. eps_B) then
         Tangent(idx,:) = (/0.0_r8,0.0_r8,0.0_r8,0.0_r8/)
      else
         Tangent(idx,1:3) = bvec(idx,:)/binter(idx)
         Tangent(idx,4)   = binter(idx)
      end if
    enddo
  end subroutine diffLine_B

  subroutine cal_diffLine_diffB(Posi,dBxOut,dByOut,dBzOut,Tan0,nnn)
    real(kind = r8)::Posi(3)
    real(kind = r8)::dBxOut(3),dByOut(3),dBzOut(3)
    real(kind = r8)::Tan0(4)
    real(kind = r8)::Tan1(4),Tan2(4),Tan3(4)
    real(kind = r8)::Tan4(4),Tan5(4),Tan6(4)
    real(kind=r8)::posi_rtp(3),posi_rtp_tmp(3)
    real(kind=r8)::drtp(3),drtp_px(3),drtp_nx(3),drtp_py(3),drtp_ny(3),drtp_pz(3),drtp_nz(3)
    integer      ::nnn(3),n_px(3),n_nx(3),n_py(3),n_ny(3),n_pz(3),n_nz(3)

    call car2polar(Posi,posi_rtp)
    call find_index(posi_rtp,nnn,drtp)

    n_px=nnn
    n_nx=nnn
    n_py=nnn
    n_ny=nnn
    n_pz=nnn
    n_nz=nnn

    call car2polar(Posi+dx,posi_rtp_tmp)
    ! call remap_index(posi_rtp_tmp,nnn,n_px,drtp_px)
    call find_index(posi_rtp_tmp,n_px,drtp_px)

    call car2polar(Posi-dx,posi_rtp_tmp)
    ! call remap_index(posi_rtp_tmp,nnn,n_nx,drtp_nx)
    call find_index(posi_rtp_tmp,n_nx,drtp_nx)

    call car2polar(Posi+dy,posi_rtp_tmp)
    ! call remap_index(posi_rtp_tmp,nnn,n_py,drtp_py)
    call find_index(posi_rtp_tmp,n_py,drtp_py)

    call car2polar(Posi-dy,posi_rtp_tmp)
    ! call remap_index(posi_rtp_tmp,nnn,n_ny,drtp_ny)
    call find_index(posi_rtp_tmp,n_ny,drtp_ny)

    call car2polar(Posi+dz,posi_rtp_tmp)
    ! call remap_index(posi_rtp_tmp,nnn,n_pz,drtp_pz)
    call find_index(posi_rtp_tmp,n_pz,drtp_pz)

    call car2polar(Posi-dz,posi_rtp_tmp)
    ! call remap_index(posi_rtp_tmp,nnn,n_nz,drtp_nz)
    call find_index(posi_rtp_tmp,n_nz,drtp_nz)

    call diffLine_A(nnn,drtp,Tan0)
    call diffLine_A(n_px,drtp_px,Tan1)
    call diffLine_A(n_nx,drtp_nx,Tan2)
    call diffLine_A(n_py,drtp_py,Tan3)
    call diffLine_A(n_ny,drtp_ny,Tan4)
    call diffLine_A(n_pz,drtp_pz,Tan5)
    call diffLine_A(n_nz,drtp_nz,Tan6)

    dBxOut(1) = Tan1(1) - Tan2(1)
    dBxOut(2) = Tan3(1) - Tan4(1)
    dBxOut(3) = Tan5(1) - Tan6(1)
    
    dByOut(1) = Tan1(2) - Tan2(2)
    dByOut(2) = Tan3(2) - Tan4(2)
    dByOut(3) = Tan5(2) - Tan6(2)
    
    dBzOut(1) = Tan1(3) - Tan2(3)
    dBzOut(2) = Tan3(3) - Tan4(3)
    dBzOut(3) = Tan5(3) - Tan6(3)

    dBxOut = dBxOut/(dxyz*2.0_r8)
    dByOut = dByOut/(dxyz*2.0_r8)
    dBzOut = dBzOut/(dxyz*2.0_r8)
  end subroutine cal_diffLine_diffB

  subroutine cal_diffLine_diffB_B(Posi,dBxOut,dByOut,dBzOut,Tan0,nnn0)
    real(kind = r8)::Posi(3)
    real(kind = r8)::dBxOut(3),dByOut(3),dBzOut(3)
    real(kind = r8)::Tan0(4),Tangent(7,4),drtp(7,3)
    real(kind=r8)::posi_rtp(3),posi_rtp_tmp(3)
    integer      ::nnn(7,3),nnn0(3),idx

    do idx=1,7
      call car2polar(Posi+dxyz_vec(idx,:),posi_rtp)
      call find_index(posi_rtp,nnn(idx,:),drtp(idx,:))
    end do

    ! call car2polar(Posi,posi_rtp)
    ! call find_index(posi_rtp,nnn(1,:),drtp(1,:))

    ! call car2polar(Posi+dx,posi_rtp_tmp)
    ! call find_index(posi_rtp_tmp,nnn(2,:),drtp(2,:))

    ! call car2polar(Posi-dx,posi_rtp_tmp)
    ! call find_index(posi_rtp_tmp,nnn(3,:),drtp(3,:))

    ! call car2polar(Posi+dy,posi_rtp_tmp)
    ! call find_index(posi_rtp_tmp,nnn(4,:),drtp(4,:))

    ! call car2polar(Posi-dy,posi_rtp_tmp)
    ! call find_index(posi_rtp_tmp,nnn(5,:),drtp(5,:))

    ! call car2polar(Posi+dz,posi_rtp_tmp)
    ! call find_index(posi_rtp_tmp,nnn(6,:),drtp(6,:))

    ! call car2polar(Posi-dz,posi_rtp_tmp)
    ! call find_index(posi_rtp_tmp,nnn(7,:),drtp(7,:))

    call diffLine_B(nnn,drtp,Tangent)
    nnn0=nnn(1,:)

    dBxOut(1) = Tangent(2,1) - Tangent(3,1)
    dBxOut(2) = Tangent(4,1) - Tangent(5,1)
    dBxOut(3) = Tangent(6,1) - Tangent(7,1)
    
    dByOut(1) = Tangent(2,2) - Tangent(3,2)
    dByOut(2) = Tangent(4,2) - Tangent(5,2)
    dByOut(3) = Tangent(6,2) - Tangent(7,2)
    
    dBzOut(1) = Tangent(2,3) - Tangent(3,3)
    dBzOut(2) = Tangent(4,3) - Tangent(5,3)
    dBzOut(3) = Tangent(6,3) - Tangent(7,3)

    dBxOut = dBxOut/(dxyz*2.0_r8)
    dByOut = dByOut/(dxyz*2.0_r8)
    dBzOut = dBzOut/(dxyz*2.0_r8)
    Tan0 = Tangent(1,:)
  end subroutine cal_diffLine_diffB_B
  subroutine remap_index(posi_rtp_in,nnn,n_in,drtp)
    real(kind=r8),intent(in)::posi_rtp_in(3)
    real(kind=r8),intent(out)::drtp(3)
    integer,intent(in) :: nnn(3)
    integer,intent(out) :: n_in(3)

    if (posi_rtp_in(1).lt.ra(nnn(1))) then
      n_in(1)=nnn(1)-1
    endif
    if (posi_rtp_in(1).gt.ra(nnn(1)+1)) then
      n_in(1)=nnn(1)+1
    endif
    drtp(1)= (posi_rtp_in(1) - ra(n_in(1))) / (ra(n_in(1) + 1) - ra(n_in(1)))

    if (posi_rtp_in(2).lt.th(nnn(2))) then
      n_in(2)=nnn(2)-1
    endif
    if (posi_rtp_in(2).gt.th(nnn(2)+1)) then
      n_in(2)=nnn(2)+1
    endif
    drtp(2)= (posi_rtp_in(2) - th(n_in(2))) / (th(n_in(2) + 1) - th(n_in(2)))

    if (nnn(3) .eq. 1 .and. posi_rtp_in(3).gt. (min_ph+max_ph)*0.5) then
        n_in(3)=dim_ph-1
        drtp(3)= (posi_rtp_in(3) - ph(n_in(3))) / (ph(n_in(3) + 1) - ph(n_in(3)))
    else if (nnn(3) .eq. dim_ph-1 .and. posi_rtp_in(3).lt. (min_ph+max_ph)*0.5) then
        n_in(3)=0
        drtp(3)= (posi_rtp_in(3) - ph(n_in(3))) / (ph(n_in(3) + 1) - ph(n_in(3)))
    else
      if (posi_rtp_in(3).lt. min_ph) then
        n_in(3)=0
        drtp(3) = 0.0_r8
      else if (posi_rtp_in(3).ge. max_ph) then
        n_in(3)=dim_ph
        drtp(3) = 0.0_r8
      else
        if (posi_rtp_in(3).lt.ph(nnn(3))) then
          n_in(3)=nnn(3)-1
        endif
        if (posi_rtp_in(3).gt.ph(nnn(3)+1)) then
          n_in(3)=nnn(3)+1
        endif
        drtp(3)= (posi_rtp_in(3) - ph(n_in(3))) / (ph(n_in(3) + 1) - ph(n_in(3)))
      end if
    end if
  end subroutine remap_index

  subroutine obtain_ds(Posi,dstep)
    real(kind=r8),intent(in)::Posi(3)
    real(kind=r8),intent(out)::dstep
    real(kind=r8)::posi_rtp(3)
    real(kind=r8)::drtp(3)
    real(kind=r8)::ds1,ds2,ds3,ds4,ds5,ds6,ds7
    integer      ::nnn(3)

    call car2polar(Posi,posi_rtp)
    call find_index(posi_rtp,nnn,drtp)

    ds1=ra(nnn(1)+1)-ra(nnn(1))
    ds2=ra(nnn(1))*(th(nnn(2)+1)-th(nnn(2)))
    ds3=ra(nnn(1)+1)*(th(nnn(2)+1)-th(nnn(2)))

    if (nnn(2) .eq. 1) then
      ds6=ra(nnn(1))*dsin(th(nnn(2)+1))*(ph(nnn(3)+1)-ph(nnn(3)))
      ds7=ra(nnn(1)+1)*dsin(th(nnn(2)+1))*(ph(nnn(3)+1)-ph(nnn(3)))
      dstep=minval((/ds1,ds2,ds3,ds6,ds7/))
    else if (nnn(2) .eq. dim_th-1) then
      ds4=ra(nnn(1))*dsin(th(nnn(2)))*(ph(nnn(3)+1)-ph(nnn(3)))
      ds5=ra(nnn(1)+1)*dsin(th(nnn(2)))*(ph(nnn(3)+1)-ph(nnn(3)))
      dstep=minval((/ds1,ds2,ds3,ds4,ds5/))
    else
      ds4=ra(nnn(1))*dsin(th(nnn(2)))*(ph(nnn(3)+1)-ph(nnn(3)))
      ds5=ra(nnn(1)+1)*dsin(th(nnn(2)))*(ph(nnn(3)+1)-ph(nnn(3)))
      ds6=ra(nnn(1))*dsin(th(nnn(2)+1))*(ph(nnn(3)+1)-ph(nnn(3)))
      ds7=ra(nnn(1)+1)*dsin(th(nnn(2)+1))*(ph(nnn(3)+1)-ph(nnn(3)))
      dstep=minval((/ds1,ds2,ds3,ds4,ds5,ds6,ds7/))
    end if

  end subroutine obtain_ds

  subroutine obtain_ds_n(nnn,dstep)
    integer,intent(in)::nnn(3)
    real(kind=r8),intent(out)::dstep
    real(kind=r8)::ds1,ds2,ds3,ds4,ds5,ds6,ds7

    ds1=ra(nnn(1)+1)-ra(nnn(1))
    ds2=ra(nnn(1))*(th(nnn(2)+1)-th(nnn(2)))
    ds3=ra(nnn(1)+1)*(th(nnn(2)+1)-th(nnn(2)))

    if (nnn(2) .eq. 1) then
      ds6=ra(nnn(1))*dsin(th(nnn(2)+1))*(ph(nnn(3)+1)-ph(nnn(3)))
      ds7=ra(nnn(1)+1)*dsin(th(nnn(2)+1))*(ph(nnn(3)+1)-ph(nnn(3)))
      dstep=minval((/ds1,ds2,ds3,ds6,ds7/))
    else if (nnn(2) .eq. dim_th-1) then
      ds4=ra(nnn(1))*dsin(th(nnn(2)))*(ph(nnn(3)+1)-ph(nnn(3)))
      ds5=ra(nnn(1)+1)*dsin(th(nnn(2)))*(ph(nnn(3)+1)-ph(nnn(3)))
      dstep=minval((/ds1,ds2,ds3,ds4,ds5/))
    else
      ds4=ra(nnn(1))*dsin(th(nnn(2)))*(ph(nnn(3)+1)-ph(nnn(3)))
      ds5=ra(nnn(1)+1)*dsin(th(nnn(2)))*(ph(nnn(3)+1)-ph(nnn(3)))
      ds6=ra(nnn(1))*dsin(th(nnn(2)+1))*(ph(nnn(3)+1)-ph(nnn(3)))
      ds7=ra(nnn(1)+1)*dsin(th(nnn(2)+1))*(ph(nnn(3)+1)-ph(nnn(3)))
      dstep=minval((/ds1,ds2,ds3,ds4,ds5,ds6,ds7/))
    end if

  end subroutine obtain_ds_n
! ---- the ODE solver runge-kutta 4 order method -----
  subroutine rk4 (f, neqn, t, t_out, y, flag, nidx)
    integer ( kind = i4 ) neqn
    external f
    real(kind = r8) :: t,t_out,dt
    real(kind = r8) :: y(neqn)
    real(kind = r8) :: y_tmp(neqn)
    real(kind = r8) :: k1(neqn)
    real(kind = r8) :: k2(neqn)
    real(kind = r8) :: k3(neqn)
    real(kind = r8) :: k4(neqn)
    real(kind = r8) :: yp(neqn)
    real(kind = r8) :: Tangent(4)
    integer :: flag
    integer :: iterNum
    integer :: nidx(3)
        
    iterNum = 0
    flag = 1
    dt = t_out - t
    y_tmp = y
    
    call f(t,y,yp, nidx)
    k1 = dt*yp
    
    call f(t+0.5*dt,y + 0.5*k1,yp, nidx)
    k2 = dt*yp
    
    call f(t+0.5*dt,y+0.5*k2,yp, nidx)
    k3 = dt*yp
    
    call f(t+dt,y+k3,yp, nidx)
    k4 = dt*yp
    
    y_tmp =  y &
        & +( k1 &
        & +  k2*2.0d0 &
        & +  k3*2.0d0 &
        & +  k4 )/6.0d0
    
    do while(outside_boundary(y_tmp(1:3)) .and. iterNum .le. iterMax)
    
       iterNum = iterNum + 1
       dt = dt * 0.5
       
       call f(t,y,yp, nidx)
       k1 = dt*yp
    
       call f(t+0.5*dt,y + 0.5*k1,yp, nidx)
       k2 = dt*yp
    
       call f(t+0.5*dt,y+0.5*k2,yp, nidx)
       k3 = dt*yp
    
       call f(t+dt,y+k3,yp, nidx)
       k4 = dt*yp
    
       y_tmp =  y &
           & +( k1 &
           & +  k2*2.0d0 &
           & +  k3*2.0d0 &
           & +  k4 )/6.0d0
    end do
    
    y = y_tmp
    call f(t,y,yp, nidx)
    
    ! call diffLine(y(1:3),Tangent)
    
    ! if(iterNum .ge. iterMax .or. Tangent(4) .le. eps_B) then
    !   flag = 0
    ! end if

    if(iterNum .ge. iterMax) then
      flag = 0
    end if
    t_out = t + dt
  end subroutine rk4
! ----------- RK3 code ----------
  subroutine rk3 (f, neqn, t, t_out, y, flag, nidx)
      integer(kind=i4) :: neqn
      external f
      real(kind=r8) :: t, t_out, dt
      real(kind=r8) :: y(neqn)
      real(kind=r8) :: y_tmp(neqn)
      real(kind=r8) :: k1(neqn), k2(neqn), k3(neqn)
      real(kind=r8) :: yp(neqn)
      integer :: flag, iterNum
      integer :: nidx(3)
  
      iterNum = 0
      flag = 1
      dt = t_out - t
      y_tmp = y
      
      ! Stage 1
      call f(t, y, yp, nidx)
      k1 = dt * yp
      
      ! Stage 2
      call f(t + 0.5 * dt, y + 0.5 * k1, yp, nidx)
      k2 = dt * yp
      
      ! Stage 3
      call f(t + dt, y - k1 + 2.0 * k2, yp, nidx)
      k3 = dt * yp
      
      ! Update the solution
      y_tmp = y + (k1 / 6.0 + 2.0 * k2 / 3.0 + k3 / 6.0)
      
      ! Handle boundary conditions and potential adjustment of the time step
      do while(outside_boundary(y_tmp(1:3)) .and. iterNum .le. iterMax)
          iterNum = iterNum + 1
          dt = dt * 0.5
          
          call f(t, y, yp, nidx)
          k1 = dt * yp
          
          call f(t + 0.5 * dt, y + 0.5 * k1, yp, nidx)
          k2 = dt * yp
          
          call f(t + dt, y - k1 + 2.0 * k2, yp, nidx)
          k3 = dt * yp
          
          y_tmp = y + (k1 / 6.0 + 2.0 * k2 / 3.0 + k3 / 6.0)
      end do
      
      y = y_tmp
      call f(t, y, yp, nidx)
      
      if (iterNum >= iterMax) then
          flag = 0
      endif
      
      t_out = t + dt
  end subroutine rk3
! ----------- RK2 code ----------
  subroutine rk2 (f, neqn, t, t_out, y, flag, nidx)
      integer(kind=i4) :: neqn
      external f
      real(kind=r8) :: t, t_out, dt
      real(kind=r8) :: y(neqn)
      real(kind=r8) :: y_tmp(neqn)
      real(kind=r8) :: k1(neqn), k2(neqn)
      real(kind=r8) :: yp(neqn)
      integer :: flag, iterNum
      integer :: nidx(3)
  
      iterNum = 0
      flag = 1
      dt = t_out - t
      y_tmp = y
      
      ! Compute k1 at the initial point
      call f(t, y, yp, nidx)
      k1 = dt * yp
      
      ! Compute k2 at the end point using the slope found from k1
      call f(t + dt*0.5, y + k1*0.5, yp, nidx)
      k2 = dt * yp
      
      ! Update the solution using the average of k1 and k2
      y_tmp = y + k2
      
      ! Check for boundary conditions and adjust time step if needed
      do while(outside_boundary(y_tmp(1:3)) .and. iterNum .le. iterMax)
          iterNum = iterNum + 1
          dt = dt * 0.5
          
          call f(t, y, yp, nidx)
          k1 = dt * yp
          
          call f(t + dt*0.5, y + k1*0.5, yp, nidx)
          k2 = dt * yp
          
          y_tmp = y + k2
      end do
      
      y = y_tmp
      call f(t, y, yp, nidx)
      
      if (iterNum >= iterMax) then
          flag = 0
      endif
      
      t_out = t + dt
  end subroutine rk2
  subroutine rk2_ralston (f, neqn, t, t_out, y, flag, nidx)
      integer(kind=i4) :: neqn
      external f
      real(kind=r8) :: t, t_out, dt
      real(kind=r8) :: y(neqn)
      real(kind=r8) :: y_tmp(neqn)
      real(kind=r8) :: k1(neqn), k2(neqn)
      real(kind=r8) :: yp(neqn)
      integer :: flag, iterNum
      integer :: nidx(3)
  
      iterNum = 0
      flag = 1
      dt = t_out - t
      y_tmp = y
      
      ! Compute k1 at the initial point
      call f(t, y, yp, nidx)
      k1 = dt * yp
      
      ! Compute k2 at the end point using the slope found from k1
      call f(t + dt*2.0/3.0, y + k1*2.0/3.0, yp, nidx)
      k2 = dt * yp
      
      ! Update the solution using the average of k1 and k2
      y_tmp = y + (k1*0.25 + k2*0.75)
      
      ! Check for boundary conditions and adjust time step if needed
      do while(outside_boundary(y_tmp(1:3)) .and. iterNum .le. iterMax)
          iterNum = iterNum + 1
          dt = dt * 0.5
          
          call f(t, y, yp, nidx)
          k1 = dt * yp
          
          call f(t + dt*2.0/3.0, y + k1*2.0/3.0, yp, nidx)
          k2 = dt * yp
          
          y_tmp = y + (k1*0.25 + k2*0.75)
      end do
      
      y = y_tmp
      call f(t, y, yp, nidx)
      
      if (iterNum >= iterMax) then
          flag = 0
      endif
      
      t_out = t + dt
  end subroutine rk2_ralston

  subroutine rk2_heun (f, neqn, t, t_out, y, flag, nidx)
      integer(kind=i4) :: neqn
      external f
      real(kind=r8) :: t, t_out, dt
      real(kind=r8) :: y(neqn)
      real(kind=r8) :: y_tmp(neqn)
      real(kind=r8) :: k1(neqn), k2(neqn)
      real(kind=r8) :: yp(neqn)
      integer :: flag, iterNum
      integer :: nidx(3)
  
      iterNum = 0
      flag = 1
      dt = t_out - t
      y_tmp = y
      
      ! Compute k1 at the initial point
      call f(t, y, yp, nidx)
      k1 = dt * yp
      
      ! Compute k2 at the end point using the slope found from k1
      call f(t + dt, y + k1, yp, nidx)
      k2 = dt * yp
      
      ! Update the solution using the average of k1 and k2
      y_tmp = y + 0.5 * (k1 + k2)
      
      ! Check for boundary conditions and adjust time step if needed
      do while(outside_boundary(y_tmp(1:3)) .and. iterNum .le. iterMax)
          iterNum = iterNum + 1
          dt = dt * 0.5
          
          call f(t, y, yp, nidx)
          k1 = dt * yp
          
          call f(t + dt, y + k1, yp, nidx)
          k2 = dt * yp
          
          y_tmp = y + 0.5 * (k1 + k2)
      end do
      
      y = y_tmp
      call f(t, y, yp, nidx)
      
      if (iterNum >= iterMax) then
          flag = 0
      endif
      
      t_out = t + dt
  end subroutine rk2_heun

! ----------- R.H.S. of the ODE ----------
  subroutine rhs(t,y,yp, nidx)
    integer,parameter :: neqn = 9
    real(kind = r8),intent(in)::t
    real(kind = r8)::y(neqn)
    real(kind = r8)::yp(neqn)
    real(kind = r8)::TangentB(4)
    real(kind = r8)::dBx(3),dBy(3),dBz(3)
    integer, intent(out):: nidx(3)
    ! call diffB(y(1:3),dBx,dBy,dBz,TangentB)

    ! call diffLine(y(1:3),TangentB)
    ! call cal_diffB(y(1:3),dBx,dBy,dBz)

    call cal_diffLine_diffB_B(y(1:3),dBx,dBy,dBz,TangentB,nidx)
    ! call cal_diffLine_diffB(y(1:3),dBx,dBy,dBz,TangentB,nidx)

    yp(1:3) = TangentB(1:3)
    yp(4) = dot_product(y(4:6),dBx)
    yp(5) = dot_product(y(4:6),dBy)
    yp(6) = dot_product(y(4:6),dBz)
    
    yp(7) = dot_product(y(7:9),dBx)
    yp(8) = dot_product(y(7:9),dBy)
    yp(9) = dot_product(y(7:9),dBz)
  end subroutine rhs

  subroutine initializeUV(LineP)
    real(kind=r8),intent(inout)::LineP(13)
    real(kind=r8)::Posi(3)
    real(kind=r8)::Tangent(4)
    real(kind=r8)::UVector(3)
    real(kind=r8)::VVector(3)
    real(kind=r8)::norm
    
    Posi = LineP(1:3)
    call diffLine(Posi,Tangent)

    if (abs(Tangent(3)).lt.eps_N) then
      if(abs(Tangent(1)).lt.eps_N) then
        VVector(1)=real(1,kind=r8)
        VVector(2)=-Tangent(1)/Tangent(2)
        VVector(3)=real(0,kind=r8)
      else
        VVector(1)=-Tangent(2)/Tangent(1)
        VVector(2)=real(1,kind=r8)
        VVector(3)=real(0,kind=r8)
      end if
    else
      VVector(1)=real(0,kind=r8)
      VVector(2)=real(1,kind=r8)
      VVector(3)=-Tangent(2)/Tangent(3)
    end if

    UVector=cross(Tangent(1:3),VVector)

    norm = sqrt(dot_product(VVector,VVector))
    LineP(7:9) = VVector/norm
        
    norm = sqrt(dot_product(UVector,UVector))
    LineP(4:6) = UVector / norm
    
  end subroutine initializeUV

! --------- search for the index along each direction -------
  subroutine find_index(posi, nnn, drtp)
      real(kind=r8), intent(in) :: posi(3)
      real(kind=r8), intent(out) :: drtp(3)
      integer, intent(out) :: nnn(3)

      if (uniform_ra) then
        call find_index_all_interpolation(posi, nnn, drtp)
      else if (uniform_th .and.uniform_ph) then
        call find_index_th_ph_interpolation_ra_binary(posi, nnn, drtp)
      else
        call find_index_all_binary(posi, nnn, drtp)
      end if

  end subroutine find_index

  ! interpolation searching
  subroutine find_index_all_interpolation(posi, nnn, drtp)
      real(kind=r8), intent(in) :: posi(3)
      real(kind=r8), intent(out) :: drtp(3)
      integer, intent(out) :: nnn(3)
      integer :: lo, hi, mid
  
      if (posi(1) < min_ra) then
          nnn(1) = 0
          drtp(1) = 0.0_r8
      else if (posi(1) >= max_ra) then
          nnn(1) = dim_ra
          drtp(1) = 0.0_r8
      else
          mid = floor((posi(1)-min_ra)/d_ra+1)
          nnn(1)=mid
          drtp(1) = (posi(1) - ra(mid)) / (ra(mid + 1) - ra(mid))
      end if

      if (posi(2) < min_th) then
          nnn(2) = 0
          drtp(2) = 0.0_r8
      else if (posi(2) >= max_th) then
          nnn(2) = dim_th
          drtp(2) = 0.0_r8
      else
          mid = floor((posi(2)-min_th)/d_th+1)
          nnn(2)=mid
          drtp(2) = (posi(2) - th(mid)) / (th(mid + 1) - th(mid))
      end if

      if (posi(3) < min_ph) then
          nnn(3) = 0
          drtp(3) = 0.0_r8
      else if (posi(3) >= max_ph) then
          nnn(3) = dim_ph
          drtp(3) = 0.0_r8
      else
          mid = floor((posi(3)-min_ph)/d_ph+1)
          nnn(3)=mid
          drtp(3) = (posi(3) - ph(mid)) / (ph(mid + 1) - ph(mid))
      end if

  end subroutine find_index_all_interpolation

  ! binary searching
  subroutine find_index_all_binary(posi, nnn, drtp)
      real(kind=r8), intent(in) :: posi(3)
      real(kind=r8), intent(out) :: drtp(3)
      integer, intent(out) :: nnn(3)
      integer :: lo, hi, mid
  
      ! Binary search for the first dimension
      if (posi(1) < min_ra) then
          nnn(1) = 0
          drtp(1) = 0.0_r8
      else if (posi(1) >= max_ra) then
          nnn(1) = dim_ra
          drtp(1) = 0.0_r8
      else
          lo = 1
          hi = dim_ra
          do while (lo <= hi)
              mid = lo + (hi - lo) / 2
              if (ra(mid) > posi(1)) then
                  hi = mid - 1
              else if (ra(mid + 1) <= posi(1)) then
                  lo = mid + 1
              else
                  nnn(1) = mid
                  drtp(1) = (posi(1) - ra(mid)) / (ra(mid + 1) - ra(mid))
                  exit
              end if
          end do
      end if
  
      ! Repeat binary search logic for the second and third dimensions
      ! This is for the second dimension (th array)
      if (posi(2) < min_th) then
          nnn(2) = 0
          drtp(2) = 0.0_r8
      else if (posi(2) >= max_th) then
          nnn(2) = dim_th
          drtp(2) = 0.0_r8
      else
          lo = 1
          hi = dim_th
          do while (lo <= hi)
              mid = lo + (hi - lo) / 2
              if (th(mid) > posi(2)) then
                  hi = mid - 1
              else if (th(mid + 1) <= posi(2)) then
                  lo = mid + 1
              else
                  nnn(2) = mid
                  drtp(2) = (posi(2) - th(mid)) / (th(mid + 1) - th(mid))
                  exit
              end if
          end do
      end if
  
      ! This is for the third dimension (ph array)
      if (posi(3) < min_ph) then
          nnn(3) = 0
          drtp(3) = 0.0_r8
      else if (posi(3) >= max_ph) then
          nnn(3) = dim_ph
          drtp(3) = 0.0_r8
      else
          lo = 1
          hi = dim_ph
          do while (lo <= hi)
              mid = lo + (hi - lo) / 2
              if (ph(mid) > posi(3)) then
                  hi = mid - 1
              else if (ph(mid + 1) <= posi(3)) then
                  lo = mid + 1
              else
                  nnn(3) = mid
                  drtp(3) = (posi(3) - ph(mid)) / (ph(mid + 1) - ph(mid))
                  exit
              end if
          end do
      end if
  
  end subroutine find_index_all_binary

  ! interpolation in theta and phi, binary in radius
  subroutine find_index_th_ph_interpolation_ra_binary(posi, nnn, drtp)
      real(kind=r8), intent(in) :: posi(3)
      real(kind=r8), intent(out) :: drtp(3)
      integer, intent(out) :: nnn(3)
      integer :: lo, hi, mid
  
      ! Binary search for the first dimension
      if (posi(1) < min_ra) then
          nnn(1) = 0
          drtp(1) = 0.0_r8
      else if (posi(1) >= max_ra) then
          nnn(1) = dim_ra
          drtp(1) = 0.0_r8
      else
          lo = 1
          hi = dim_ra
          do while (lo <= hi)
              mid = lo + (hi - lo) / 2
              if (ra(mid) > posi(1)) then
                  hi = mid - 1
              else if (ra(mid + 1) <= posi(1)) then
                  lo = mid + 1
              else
                  nnn(1) = mid
                  drtp(1) = (posi(1) - ra(mid)) / (ra(mid + 1) - ra(mid))
                  exit
              end if
          end do
      end if
  
      if (posi(2) < min_th) then
          nnn(2) = 0
          drtp(2) = 0.0_r8
      else if (posi(2) >= max_th) then
          nnn(2) = dim_th
          drtp(2) = 0.0_r8
      else
          mid = floor((posi(2)-min_th)/d_th+1)
          nnn(2)=mid
          drtp(2) = (posi(2) - th(mid)) / (th(mid + 1) - th(mid))
      end if

      if (posi(3) < min_ph) then
          nnn(3) = 0
          drtp(3) = 0.0_r8
      else if (posi(3) >= max_ph) then
          nnn(3) = dim_ph
          drtp(3) = 0.0_r8
      else
          mid = floor((posi(3)-min_ph)/d_ph+1)
          nnn(3)=mid
          drtp(3) = (posi(3) - ph(mid)) / (ph(mid + 1) - ph(mid))
      end if
  
  end subroutine find_index_th_ph_interpolation_ra_binary

  ! subroutine find_index_old(posi,nnn,drtp)
  !   real(kind=r8),intent(in) ::posi(3)
  !   real(kind=r8),intent(out)::drtp(3)
  !   integer,intent(out) :: nnn(3)
  !   integer :: i

  !   if (posi(1) .lt. min_ra) then
  !      nnn(1) = 0
  !      drtp(1) = real(0,kind=r8)
  !   else if (posi(1) .ge. max_ra) then
  !      nnn(1) = dim_ra
  !      drtp(1) = real(0,kind=r8)
  !   else
  !      do i=1,dim_ra
  !         if(ra(i) .le. posi(1) .and. ra(i+1) .gt. posi(1) ) then
  !           nnn(1) = i
  !           drtp(1) = (posi(1) - ra(i))/( ra(i+1) - ra(i) )
  !           exit
  !         end if
  !      end do
  !   end if
    
  !   if (posi(2) .lt. min_th) then
  !      nnn(2) = 0
  !      drtp(2) = real(0,kind=r8)
  !   else if (posi(2) .ge. max_th) then
  !      nnn(2) = dim_th
  !      drtp(2) = real(0,kind=r8)
  !   else
  !      do i=1,dim_th
  !         if(th(i) .le. posi(2) .and. th(i+1) .gt. posi(2) ) then
  !           nnn(2) = i
  !           drtp(2) = (posi(2) - th(i))/( th(i+1) - th(i) )
  !           exit
  !         end if
  !      end do
  !   end if

  !   if (posi(3) .lt. min_ph) then
  !      nnn(3) = 0
  !      drtp(3) = real(0,kind=r8)
  !   else if (posi(3) .ge. max_ph) then
  !      nnn(3) = dim_ph
  !      drtp(3) = real(0,kind=r8)
  !   else
  !      do i=1,dim_ph
  !         if(ph(i) .le. posi(3) .and. ph(i+1) .gt. posi(3) ) then
  !           nnn(3) = i
  !           drtp(3) = (posi(3) - ph(i))/( ph(i+1) - ph(i) )
  !           exit
  !         end if
  !      end do
  !   end if
  ! end subroutine find_index_old

end module mod_operator
