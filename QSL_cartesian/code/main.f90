!=============================================================
! This program is written in Fortran
! @Author: kaiyang 
! @Date:   2019-01-11 13:06:36
! @Email:  yangkaijilin@gmail.com
! @Last Modified by:   kaiyang
! @Last Modified time: 2019-11-15 20:57:48
! @Version: 
! @Licensing: Copy right from Kai E. Yang
! @Purpose: calculate the twist number by using the parallel 
!           current intrgral.
!=============================================================
program main
  use mod_param
  use mod_io
  use mod_solver
  use mod_operator
  implicit none

  call get_command_argument(1,par)
  call read_parameter()
  call allocate_var()
  call read_data()
  call show_information()
  write(*,'(A)')'| The calculation start!'
  call do_cal()
  write(*,'(A)')'| The calculation finished!'
  call write_data()
  call deallocate_var()

end program main
