# QSL
 
The program main.f90 is the main program.
This program is used for calculating the squashing factor Q.
The method is proposed by Pariat & Demoulin 2012.
This version is usually used by myself.
Version v2.0
History:
	K.Y. 2016@NJU
	K.Y. Jun/2017@MSU change the compiler from ifort to gfortran.
	K.Y. Mar/2017@MSU change the reading file from the initial formatted to unformatted.
	K.Y. Jun/2018@NJU change the most of names of the subroutines to make it look better.
	K.Y. Oct/2018@NJU use Scott's method in ApJ 2017
	K.Y. Nov/2018@NJU add the namelist format for reading the parameters.
   
The steps for this program is:

1	creat the magnetic field data with the program output_b.pro

2	> make -f makefile

3	copy the executable file "qsl" and the parameter file "par" to the document 'data' where the magnetic field data locates.

4 	change the input parameter in "par" by hand to select the subregion, the level of grid refinement, the name of the output file, or the integral step. The defoult level of grid refinement is 4, the default computational region is the whole region, the defoult value of the integral step is 0.25.

-----------------------------------------
the format of the par file is namelist in fortran, e.g.,

&filename
	BfieldName = '*****.binary'  ! file name for the Bxyz companent 
	OutFileName = '*****.binary' ! file name for the out put file
/
&cal_par
	nthreads = 1                 ! thread number used for parallel computation
	dimx = 50                    ! dimension for x
	dimy = 50                    ! dimension for y
	dimz = 50                    ! dimension for z
	x_start = 1                  ! start point in x (it starts from 1 to dim X)
	y_start = 1                  ! start point in y (it starts from 1 to dim Y)
	z_start = 1                  ! start point in z (it starts from 1 to dim Z)
	x_end = 50                   ! end point in x (it starts from 1 to dim X)
	y_end = 50                   ! end point in y (it starts from 1 to dim Y)
	z_end = 1                    ! end point in z (it starts from 1 to dim Z)
	nlevel = 10                  ! grid level
	delta_s = 0.25               ! integral step (the default value is 0.25)
/
-----------------------------------------

5	run the program by ./qsl par (par is the parameter file)

6	use the idl program read_qsl.pro to read the result.
    IDL> .com read_qsl.pro
    IDL> result = read_qsl(par_name); par_name is a string contains the name of the par file, e.g., 'par'

In addition:

This parallel version is based on FORTRAN OPENMP.

If one want to use N threads for the calculation, just change the value of parameter 'nthreads' in parameter file. If this parameter is defined as 0, then the max number of threads in the computer will be used as default.

The field line integral method is Runge-Kutta 4(5).