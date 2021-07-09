module mod_io
use mod_param
use mod_operator
implicit none
contains
  subroutine read_parameter()
    implicit none
    ! set default value    
    BfieldName   = 'bfield.binary'
    OutFileName  = 'result.binary'
    indataformat = 'binary'

    nthreads    = 1
    dim_th        = 100
    dim_ph        = 100
    dim_ra        = 100
    
    th_start     = 100
    ph_start     = 100
    ra_start     = 100
    
    th_end       = 100
    ph_end       = 100
    ra_end       = 100
    
    nlevel      = 1
    delta_s     = 0.1d0
    ! read parameters   
    write(*,'(A)')'| = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ='   
    write(*,'(A)')'| Reading parameters from file: '//trim(par)
    open(Unit=1,file=par,STATUS='OLD',ACTION='READ')
    read(1,nml=filename_par,end=3)
3   rewind(1)
    read(1,nml=cal_par,end=5)
5   close(1)

  end subroutine read_parameter

  subroutine allocate_var()
    allocate(Bth(0:dim_ra+1,0:dim_th+1,0:dim_ph+1))
    allocate(Bph(0:dim_ra+1,0:dim_th+1,0:dim_ph+1))
    allocate(Bra(0:dim_ra+1,0:dim_th+1,0:dim_ph+1))
    
    allocate(Bx(0:dim_ra+1,0:dim_th+1,0:dim_ph+1))
    allocate(By(0:dim_ra+1,0:dim_th+1,0:dim_ph+1))
    allocate(Bz(0:dim_ra+1,0:dim_th+1,0:dim_ph+1))
    
    allocate(th(0:dim_th+1))
    allocate(ph(0:dim_ph+1))
    allocate(ra(0:dim_ra+1))
  end subroutine allocate_var

  subroutine post_allocate()
    integer :: i
    refine_dim_th = nlevel*floor((th_end - th_start)/d_th)+1
    refine_dim_ph = nlevel*floor((ph_end - ph_start)/d_ph)+1
    refine_dim_ra = nlevel*floor((ra_end - ra_start)/d_ra)+1

    allocate(refine_ra(refine_dim_ra))
    allocate(refine_ph(refine_dim_ph))
    allocate(refine_th(refine_dim_th))

    allocate(cal_data(refine_dim_ra,refine_dim_th,refine_dim_ph,3))

    refine_th = th_end+(th_start-th_end)&
            & * (/(real(i-1,kind=r8),i=1,refine_dim_th)/)&
            & / (refine_dim_th-1)
    
    refine_ph = ph_start+(ph_end-ph_start)&
            & * (/(real(i-1,kind=r8),i=1,refine_dim_ph)/)&
            & / (refine_dim_ph-1)
    
    if(refine_dim_ra .eq. 1) then 
      refine_ra = (/ra_start/)
    else
      refine_ra = ra_start+(ra_end-ra_start)&
              & * (/(real(i-1,kind=r8),i=1,refine_dim_ra)/)&
              & / (refine_dim_ra-1)
    end if

    delta_s = ra(1)*d_th*real(delta_s,kind=r8)
    BoundaryEps = real(1,kind=r8)+1e-2*delta_s
    iterMax = 20
    eps_B = epsilon(real(1,kind=r8))

  end subroutine post_allocate

  subroutine deallocate_var()
    implicit none  
    deallocate(Bth,Bph,Bra)
    deallocate(Bx,By,Bz)
    deallocate(th,ph,ra)
    deallocate(refine_ra)
    deallocate(refine_ph)
    deallocate(refine_th)
    deallocate(cal_data)
  end subroutine deallocate_var

  subroutine show_information()
    write(*,'(A)')'|--------------------------------------------------------'
    write(*,'(A)')'| Mission: Calculate the Squashing Factor Q'
    write(*,'(A)')'| The method based on Scott Roger B. et al. 2017 ApJ 848, 117.'
    write(*,'(A)')'| 3D B field data comes from: '
    write(*,'(A)')'|                    '//trim(BfieldName)
    write(*,'(A)')'| The result is saved in file: '
    write(*,'(A)')'|                    '//trim(OutFileName)
    write(*,'(A)')'| The array of the result contains three value'
    write(*,'(A)')'| The first dimension is the sign(Bra)*log10(Q)'
    write(*,'(A)')'| The second dimension is length of the field lines'
    write(*,'(A)')'| The third dimension is mark of the boundary'
    write(*,'(A)')'| The meaning of this value is:'
    write(*,'(A)')'| '
    write(*,'(A)')'| 1: close field line with both ends root in bottom boundary'
    write(*,'(A)')'| 2: open field line, positive end roots in bottom'
    write(*,'(A)')'| 3: open field line, negative end roots in bottom'
    write(*,'(A)')'| 4: field line without end roots in the bottom boundary'
    write(*,'(A)')'| '
    write(*,'(A)')'|--------------------------------------------------------'
    write(*,'(A)')'| The NDIM of the Magnetic data is: '
    write(*,'(A,I3)')'|                  ra dimension: ',dim_ra
    write(*,'(A,I3)')'|                  th dimension: ',dim_th
    write(*,'(A,I3)')'|                  ph dimension: ',dim_ph
    write(*,'(A)')'| '
    write(*,'(A)')'| The calculated domain is: '
    write(*,'(A,g10.4,A,g10.4)')'|             ra ',ra_start,' to ',ra_end
    write(*,'(A,g10.4,A,g10.4)')'|             th ',th_start,' to ',th_end
    write(*,'(A,g10.4,A,g10.4)')'|             ph ',ph_start,' to ',ph_end
    write(*,'(A,I4)')'| The refine level is: ',nlevel
    write(*,'(A)')'| The NDIM of the result is: '
    write(*,'(A,I8)')'|                  ra dimension: ',refine_dim_ra
    write(*,'(A,I8)')'|                  th dimension: ',refine_dim_th
    write(*,'(A,I8)')'|                  ph dimension: ',refine_dim_ph
    write(*,'(A,I4,A)')'| This mission use ',nthreads,' threads'
    write(*,'(A)')'|--------------------------------------------------------'
  end subroutine show_information

  subroutine read_blk()
    logical :: alive
    integer(kind=i4) :: i,j,k
    integer(kind=i8) :: ntotal
    character*1024   :: header
    double precision :: ttime

    inquire(file=BfieldName,exist=alive)
    if(alive) then 
      write(*,*)'OUTPUT_B: loading data and preparing input files for QSL3D'
      OPEN(UNIT=3,STATUS='OLD',ACTION='READ',FILE=BfieldName,POSITION='REWIND',FORM='UNFORMATTED')
      read(3) header
      read(3) ntotal,dim_ra,dim_th,dim_ph
      read(3) ttime
      do k=1,dim_ph
      do j=1,dim_th
      do i=1,dim_ra
         read(3) Bth(i,j,k),Bph(i,j,k),Bra(i,j,k)
      end do
      end do
      end do 
      close(3)
    else
      write(*,'(A)')'| File '//BfieldName//' does not exist'
      stop
    end if
  end subroutine read_blk

  subroutine read_binary()
    logical          :: alive

    inquire(file=BfieldName,exist=alive)
    if(alive) then 
      write(*,'(A)')'| Loading magnetic field data...'
      open(3,File=BfieldName,Access="stream",Form = "unformatted" )
      read(3) ra(1:dim_ra)
      read(3) th(1:dim_th)
      read(3) ph(1:dim_ph)
      read(3) Bra(1:dim_ra,1:dim_th,1:dim_ph)
      read(3) Bth(1:dim_ra,1:dim_th,1:dim_ph)
      read(3) Bph(1:dim_ra,1:dim_th,1:dim_ph)
      close(3)
      write(*,'(A)')'| Data sucessfully loaded !'
    else
      write(*,'(A)')'| File '//BfieldName//' does not exist'
      stop
    end if
  end subroutine read_binary

  subroutine read_data()
    real(kind=r8) :: vec_tmp(3)
    real(kind=r8) :: rtp(3)
    integer       :: i,j,k

    select case(indataformat)
    case('blk')
        call read_blk()
    case('binary')
        call read_binary()
    case default
        write(*,*)'no data input'
        stop
    end select

    ! --------- get the value of B in the ghost cell ------------
    ! surface in ghost shell
    Bth(0,1:dim_th,1:dim_ph) = Bth(1,1:dim_th,1:dim_ph)
    Bph(0,1:dim_th,1:dim_ph) = Bph(1,1:dim_th,1:dim_ph)
    Bra(0,1:dim_th,1:dim_ph) = Bra(1,1:dim_th,1:dim_ph)
      
    Bth(dim_ra+1,1:dim_th,1:dim_ph) = Bth(dim_ra,1:dim_th,1:dim_ph)
    Bph(dim_ra+1,1:dim_th,1:dim_ph) = Bph(dim_ra,1:dim_th,1:dim_ph)
    Bra(dim_ra+1,1:dim_th,1:dim_ph) = Bra(dim_ra,1:dim_th,1:dim_ph)
      
    Bth(1:dim_ra,0,1:dim_ph) = Bth(1:dim_ra,1,1:dim_ph)
    Bph(1:dim_ra,0,1:dim_ph) = Bph(1:dim_ra,1,1:dim_ph)
    Bra(1:dim_ra,0,1:dim_ph) = Bra(1:dim_ra,1,1:dim_ph)
      
    Bth(1:dim_ra,dim_th+1,1:dim_ph) = Bth(1:dim_ra,dim_th,1:dim_ph)
    Bph(1:dim_ra,dim_th+1,1:dim_ph) = Bph(1:dim_ra,dim_th,1:dim_ph)
    Bra(1:dim_ra,dim_th+1,1:dim_ph) = Bra(1:dim_ra,dim_th,1:dim_ph)
      
    Bth(1:dim_ra,1:dim_th,0) = Bth(1:dim_ra,1:dim_th,1)
    Bph(1:dim_ra,1:dim_th,0) = Bph(1:dim_ra,1:dim_th,1)
    Bra(1:dim_ra,1:dim_th,0) = Bra(1:dim_ra,1:dim_th,1)
      
    Bth(1:dim_ra,1:dim_th,dim_ph+1) = Bth(1:dim_ra,1:dim_th,dim_ph)
    Bph(1:dim_ra,1:dim_th,dim_ph+1) = Bph(1:dim_ra,1:dim_th,dim_ph)
    Bra(1:dim_ra,1:dim_th,dim_ph+1) = Bra(1:dim_ra,1:dim_th,dim_ph)
    
    ! lines in ghost shell
    Bth(0,0,1:dim_ph) = Bth(1,1,1:dim_ph)
    Bph(0,0,1:dim_ph) = Bph(1,1,1:dim_ph)
    Bra(0,0,1:dim_ph) = Bra(1,1,1:dim_ph)
      
    Bth(dim_ra+1,0,1:dim_ph) = Bth(dim_ra,1,1:dim_ph)
    Bph(dim_ra+1,0,1:dim_ph) = Bph(dim_ra,1,1:dim_ph)
    Bra(dim_ra+1,0,1:dim_ph) = Bra(dim_ra,1,1:dim_ph)
      
    Bth(0,dim_th+1,1:dim_ph) = Bth(1,dim_th,1:dim_ph)
    Bph(0,dim_th+1,1:dim_ph) = Bph(1,dim_th,1:dim_ph)
    Bra(0,dim_th+1,1:dim_ph) = Bra(1,dim_th,1:dim_ph)
      
    Bth(dim_ra+1,dim_th+1,1:dim_ph) = Bth(dim_ra,dim_th,1:dim_ph)
    Bph(dim_ra+1,dim_th+1,1:dim_ph) = Bph(dim_ra,dim_th,1:dim_ph)
    Bra(dim_ra+1,dim_th+1,1:dim_ph) = Bra(dim_ra,dim_th,1:dim_ph)
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    Bth(1:dim_ra,0,0) = Bth(1:dim_ra,1,1)
    Bph(1:dim_ra,0,0) = Bph(1:dim_ra,1,1)
    Bra(1:dim_ra,0,0) = Bra(1:dim_ra,1,1)
      
    Bth(1:dim_ra,dim_th+1,0) = Bth(1:dim_ra,dim_th,1)
    Bph(1:dim_ra,dim_th+1,0) = Bph(1:dim_ra,dim_th,1)
    Bra(1:dim_ra,dim_th+1,0) = Bra(1:dim_ra,dim_th,1)
      
    Bth(1:dim_ra,dim_th+1,dim_ph+1) = Bth(1:dim_ra,dim_th,dim_ph)
    Bph(1:dim_ra,dim_th+1,dim_ph+1) = Bph(1:dim_ra,dim_th,dim_ph)
    Bra(1:dim_ra,dim_th+1,dim_ph+1) = Bra(1:dim_ra,dim_th,dim_ph)
      
    Bth(1:dim_ra,0,dim_ph+1) = Bth(1:dim_ra,1,dim_ph)
    Bph(1:dim_ra,0,dim_ph+1) = Bph(1:dim_ra,1,dim_ph)
    Bra(1:dim_ra,0,dim_ph+1) = Bra(1:dim_ra,1,dim_ph)
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    Bth(0,1:dim_th,0) = Bth(1,1:dim_th,1)
    Bph(0,1:dim_th,0) = Bph(1,1:dim_th,1)
    Bra(0,1:dim_th,0) = Bra(1,1:dim_th,1)
      
    Bth(dim_ra+1,1:dim_th,0) = Bth(dim_ra,1:dim_th,1)
    Bph(dim_ra+1,1:dim_th,0) = Bph(dim_ra,1:dim_th,1)
    Bra(dim_ra+1,1:dim_th,0) = Bra(dim_ra,1:dim_th,1)
      
    Bth(0,1:dim_th,dim_ph+1) = Bth(1,1:dim_th,dim_ph)
    Bph(0,1:dim_th,dim_ph+1) = Bph(1,1:dim_th,dim_ph)
    Bra(0,1:dim_th,dim_ph+1) = Bra(1,1:dim_th,dim_ph)
      
    Bth(dim_ra+1,1:dim_th,dim_ph+1) = Bth(dim_ra,1:dim_th,dim_ph)
    Bph(dim_ra+1,1:dim_th,dim_ph+1) = Bph(dim_ra,1:dim_th,dim_ph)
    Bra(dim_ra+1,1:dim_th,dim_ph+1) = Bra(dim_ra,1:dim_th,dim_ph)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    ! points in ghost shell
    Bth(0,0,0) = Bth(1,1,1)
    Bph(0,0,0) = Bph(1,1,1)
    Bra(0,0,0) = Bra(1,1,1)
      
    Bth(dim_ra+1,0,0) = Bth(dim_ra,1,1)
    Bph(dim_ra+1,0,0) = Bph(dim_ra,1,1)
    Bra(dim_ra+1,0,0) = Bra(dim_ra,1,1)
      
    Bth(0,dim_th+1,0) = Bth(1,dim_th,1)
    Bph(0,dim_th+1,0) = Bph(1,dim_th,1)
    Bra(0,dim_th+1,0) = Bra(1,dim_th,1)
      
    Bth(0,0,dim_ph+1) = Bth(1,1,dim_ph)
    Bph(0,0,dim_ph+1) = Bph(1,1,dim_ph)
    Bra(0,0,dim_ph+1) = Bra(1,1,dim_ph)
      
    Bth(dim_ra+1,dim_th+1,0) = Bth(dim_ra,dim_th,1)
    Bph(dim_ra+1,dim_th+1,0) = Bph(dim_ra,dim_th,1)
    Bra(dim_ra+1,dim_th+1,0) = Bra(dim_ra,dim_th,1)
      
    Bth(dim_ra+1,0,dim_ph+1) = Bth(dim_ra,1,dim_ph)
    Bph(dim_ra+1,0,dim_ph+1) = Bph(dim_ra,1,dim_ph)
    Bra(dim_ra+1,0,dim_ph+1) = Bra(dim_ra,1,dim_ph)
      
    Bth(0,dim_th+1,dim_ph+1) = Bth(1,dim_th,dim_ph)
    Bph(0,dim_th+1,dim_ph+1) = Bph(1,dim_th,dim_ph)
    Bra(0,dim_th+1,dim_ph+1) = Bra(1,dim_th,dim_ph)
      
    Bth(dim_ra+1,dim_th+1,dim_ph+1) = Bth(dim_ra,dim_th,dim_ph)
    Bph(dim_ra+1,dim_th+1,dim_ph+1) = Bph(dim_ra,dim_th,dim_ph)
    Bra(dim_ra+1,dim_th+1,dim_ph+1) = Bra(dim_ra,dim_th,dim_ph)

    ra(0)            = ra(1) 
    th(0)             = th(1)
    ph(0)               = ph(1)
    ra(dim_ra+1) = ra(dim_ra)
    th(dim_th+1)   = th(dim_th)
    ph(dim_ph+1)       = ph(dim_ph)
    
    do k=0,dim_ph+1
    do j=0,dim_th+1
    do i=0,dim_ra+1
       rtp       = (/ra(i),th(j),ph(k)/)
       vec_tmp   = vec_rtp2car(rtp,Bra(i,j,k),Bth(i,j,k),Bph(i,j,k))
       Bx(i,j,k) = vec_tmp(1)
       By(i,j,k) = vec_tmp(2)
       Bz(i,j,k) = vec_tmp(3)
    end do
    end do
    end do
    
    max_th = maxval(th)
    max_ph = maxval(ph)
    max_ra = maxval(ra)

    min_th = minval(th)
    min_ph = minval(ph)
    min_ra = minval(ra)

    d_th  = (max_th - min_th)/real(dim_th - 1, kind=r8)
    d_ph  = (max_ph - min_ph)/real(dim_ph - 1, kind=r8)
    d_ra  = (max_ra - min_ra)/real(dim_ra - 1, kind=r8)

    call post_allocate()
  end subroutine read_data

  subroutine write_data()
    write(*,'(A)')'| Writing the result ...'
    open(4,File=OutFileName,Access="stream",STATUS="REPLACE",&
    &Form = "unformatted" )
    write(4) cal_data
    close(4)
  end subroutine write_data

end module mod_io
