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
      vec_rtp(3) = real(0,kind=r8)
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
  subroutine corner (nnn,xv)
    integer,intent(in)::nnn(3)
    integer::n1,n2,n3
    integer::n1p,n2p,n3p
    real(kind = r8),intent(out)::xv(3,8)
    n1 = nnn(1)
    n2 = nnn(2)
    n3 = nnn(3)
    
    n1p = nnn(1)+1
    n2p = nnn(2)+1
    n3p = nnn(3)+1

    xv(1,1) = Bx(n1,n2,n3)
    xv(1,2) = Bx(n1p,n2,n3)
    xv(1,3) = Bx(n1p,n2p,n3)
    xv(1,4) = Bx(n1,n2p,n3)
    xv(1,5) = Bx(n1,n2,n3p)
    xv(1,6) = Bx(n1p,n2,n3p)
    xv(1,7) = Bx(n1p,n2p,n3p)
    xv(1,8) = Bx(n1,n2p,n3p)
    
    xv(2,1) = By(n1,n2,n3)
    xv(2,2) = By(n1p,n2,n3)
    xv(2,3) = By(n1p,n2p,n3)
    xv(2,4) = By(n1,n2p,n3)
    xv(2,5) = By(n1,n2,n3p)
    xv(2,6) = By(n1p,n2,n3p)
    xv(2,7) = By(n1p,n2p,n3p)
    xv(2,8) = By(n1,n2p,n3p)
    
    xv(3,1) = Bz(n1,n2,n3)
    xv(3,2) = Bz(n1p,n2,n3)
    xv(3,3) = Bz(n1p,n2p,n3)
    xv(3,4) = Bz(n1,n2p,n3)
    xv(3,5) = Bz(n1,n2,n3p)
    xv(3,6) = Bz(n1p,n2,n3p)
    xv(3,7) = Bz(n1p,n2p,n3p)
    xv(3,8) = Bz(n1,n2p,n3p)
  end subroutine corner
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
       Tangent = real(0,kind=r8)*(/1,1,1,1/)
    else
       Tangent(1:3) = bvec/binter
       Tangent(4)   = binter
    end if
  end subroutine diffLine

  subroutine diffB(Posi,dBx,dBy,dBz)
    real(kind = r8)::Posi(3)
    real(kind = r8)::dBx(3),dBy(3),dBz(3)
    real(kind = r8)::Tan1(4),Tan2(4),Tan3(4)
    real(kind = r8)::Tan4(4),Tan5(4),Tan6(4)
    real(kind = r8)::dx(3),dy(3),dz(3)
    real(kind = r8)::dxyz
    
    dxyz = real(0.001,kind=r8)
    
    dx = real(0,kind=r8)*(/1,1,1/)
    dy = real(0,kind=r8)*(/1,1,1/)
    dz = real(0,kind=r8)*(/1,1,1/)
    
    dx(1) = dxyz
    dy(2) = dxyz
    dz(3) = dxyz
    
    call diffLine(Posi+dx,Tan1)
    call diffLine(Posi-dx,Tan2)
    call diffLine(Posi+dy,Tan3)
    call diffLine(Posi-dy,Tan4)
    call diffLine(Posi+dz,Tan5)
    call diffLine(Posi-dz,Tan6)
        
    dBx(1) = Tan1(1) - Tan2(1)
    dBx(2) = Tan3(1) - Tan4(1)
    dBx(3) = Tan5(1) - Tan6(1)
    
    dBy(1) = Tan1(2) - Tan2(2)
    dBy(2) = Tan3(2) - Tan4(2)
    dBy(3) = Tan5(2) - Tan6(2)
    
    dBz(1) = Tan1(3) - Tan2(3)
    dBz(2) = Tan3(3) - Tan4(3)
    dBz(3) = Tan5(3) - Tan6(3)
    
    dBx = dBx/(dxyz*real(2,kind=r8))
    dBy = dBy/(dxyz*real(2,kind=r8))
    dBz = dBz/(dxyz*real(2,kind=r8))
  end subroutine diffB

! ---- the ODE solver runge-kutta 4 order method -----
  subroutine rk4 (f, neqn, t, t_out, y, flag)
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
        
    iterNum = 0
    flag = 1
    dt = t_out - t
    y_tmp = y
    
    call f(t,y,yp)
    k1 = dt*yp
    
    call f(t+0.5*dt,y + 0.5*k1,yp)
    k2 = dt*yp
    
    call f(t+0.5*dt,y+0.5*k2,yp)
    k3 = dt*yp
    
    call f(t+dt,y+k3,yp)
    k4 = dt*yp
    
    y_tmp =  y &
        & +( k1 &
        & +  k2*real(2,kind=r8) &
        & +  k3*real(2,kind=r8) &
        & +  k4 )/6.0d0
    
    do while(outside_boundary(y_tmp(1:3)) .and. iterNum .le. iterMax)
    
       iterNum = iterNum + 1
       dt = dt * real(0.5,kind=r8)
       
       call f(t,y,yp)
       k1 = dt*yp
    
       call f(t+0.5*dt,y + 0.5*k1,yp)
       k2 = dt*yp
    
       call f(t+0.5*dt,y+0.5*k2,yp)
       k3 = dt*yp
    
       call f(t+dt,y+k3,yp)
       k4 = dt*yp
    
       y_tmp =  y &
           & +( k1 &
           & +  k2*real(2,kind=r8) &
           & +  k3*real(2,kind=r8) &
           & +  k4 )/6.0d0
    end do
    
    y = y_tmp
    call f(t,y,yp)
    
    call diffLine(y(1:3),Tangent)
    
    if(iterNum .ge. iterMax .or. Tangent(4) .le. eps_B) then
      flag = 0
    end if
    
    t_out = t + dt
  end subroutine rk4
! ----------- R.H.S. of the ODE ----------
  subroutine rhs(t,y,yp)
    integer,parameter :: neqn = 9
    real(kind = r8)::t
    real(kind = r8)::y(neqn)
    real(kind = r8)::yp(neqn)
    real(kind = r8)::TangentB(4)
    real(kind = r8)::dBx(3),dBy(3),dBz(3)
    
    call diffLine(y(1:3),TangentB)
    call diffB(y(1:3),dBx,dBy,dBz)
    
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
    VVector(1) = Tangent(2)
    VVector(2) = -1.0*Tangent(1)
    VVector(3) = real(0,kind=r8)
    VVector = VVector - dot_product(VVector,Tangent(1:3))*Tangent(1:3)
    norm = sqrt(dot_product(VVector,VVector))
    VVector = VVector/norm
    
    UVector = -1.0*cross(VVector,Tangent(1:3))
    
    norm = sqrt(dot_product(UVector,UVector))
    UVector = UVector / norm
    
    LineP(4:6) = UVector
    LineP(7:9) = VVector

  end subroutine initializeUV

! --------- search for the index along each direction -------
  subroutine find_index(posi,nnn,drtp)
    real(kind=r8),intent(in) ::posi(3)
    real(kind=r8),intent(out)::drtp(3)
    integer,intent(out) :: nnn(3)
    integer :: i
    
    if (posi(1) .lt. min_ra) then
       nnn(1) = 0
       drtp(1) = real(0,kind=r8)
    else if (posi(1) .gt. max_ra) then
       nnn(1) = dim_ra
       drtp(1) = real(0,kind=r8)
    else
       do i=1,dim_ra-1
          if(ra(i) .le. posi(1) .and. ra(i+1) .gt. posi(1) ) then
            nnn(1) = i
            drtp(1) = (posi(1) - ra(i))/( ra(i+1) - ra(i) )
            exit
          end if
       end do
    end if
    
    if (posi(2) .lt. min_th) then
       nnn(2) = 0
       drtp(2) = real(0,kind=r8)
    else if (posi(2) .gt. max_th) then
       nnn(2) = dim_th
       drtp(2) = real(0,kind=r8)
    else
       do i=1,dim_th-1
          if(th(i) .gt. posi(2) .and. th(i+1) .le. posi(2) ) then
            nnn(2) = i
            drtp(2) = (posi(2) - th(i))/( th(i+1) - th(i) )
            exit
          end if
       end do
    end if

    if (posi(3) .lt. min_ph) then
       nnn(3) = 0
       drtp(3) = real(0,kind=r8)
    else if (posi(3) .gt. max_ph) then
       nnn(3) = dim_ph
       drtp(3) = real(0,kind=r8)
    else
       do i=1,dim_ph-1
          if(ph(i) .le. posi(3) .and. ph(i+1) .gt. posi(3) ) then
            nnn(3) = i
            drtp(3) = (posi(3) - ph(i))/( ph(i+1) - ph(i) )
            exit
          end if
       end do
    end if

  end subroutine find_index

end module mod_operator
