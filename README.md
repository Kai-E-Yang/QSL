# K-QSL
A package to calculate squashing degree Q in both Cartesian and spherical coordinates

## What to do
The file main.f90 is the main program.

The file mod_io.f90 contains all the routines used for reading and writing.

The file mod_param.f90 contains all the definitions of global parameters.

The file mod_operator.f90 contains all the small parts of some calculations, like rk4, cross product ...

The file mod_solver.f90 contains the main pathway for the calculation.

This program is used for calculating the squashing factor Q.

The old method is proposed by Pariat & Demoulin (2012, A&A), which is not used anymore. A legacy can be found [here](https://github.com/njuguoyang/magnetic_modeling_codes/tree/main/code/QSL).

For the latest version, the method is that proposed in Scott (2017, ApJ) and Tassev & Savcheva (2017, ApJ).

## History
- v1 K.Y. 2016@NJU

- v2 K.Y. Jun/2017@MSU change the compiler from ifort to gfortran.

   K.Y. Mar/2017@MSU change the reading file from the initial formatted to unformatted.
   
- v3 K.Y. Jun/2018@NJU change the most of names of the subroutines to make it look better.

- v4 K.Y. Oct/2018@NJU use Scott's method in ApJ 2017

   K.Y. Nov/2018@NJU add the namelist format for reading the parameters.
   
- v7 K.Y. May/2019@SIFA use SDF lib for io, only used in this version.

- v8 K.Y. Oct/2019@SIFA collect all the files in mod format.

   K.Y. Nov/2019@SIFA correct the version in spherical coordinate system.

- v9 C.X. Apr/2021@YNU add vtk format for output data.

- Yang Guo Dec/2022@NJU fixed some bugs for the position index.

- Chloe Wilkins 2023@newcastle fixed bugs in the test field generation.

## How to use it
1, compile the program by using the makefile
```
>make -f Makefile
```

2, create a 3d magnetic field, e.g., in test1, run bxyz.pro in IDL
```
> IDL
IDL> .r bxyz.pro
```

3, make a soft link of the program at the target directory, like ./test1
```
cd ./test1
ln -s ../code/QSL
```

4, run the code with parameter file "par"
```
./QSL par
```

5, visulize the result, you can use any method you want to read the binary file for the result, I give an example by a python file "read_qsl.py"
```
>python read_qsl.py
```

## Parameter file
The format of the par file is the namelist in fortran, e.g.,

```
&filename_par
	BfieldName = '*****.binary'  ! file name for the Bxyz companent 
	OutFileName = '*****.binary' ! file name for the out put file
        indataformat = 'binary' or 'blk' ! data format of the input file
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
        x_min = -10.0                ! minimal position in x coordinate
        x_max =  10.0                ! maximal position in x coordinate
        y_min = -20.0                ! minimal position in y coordinate
        y_max =  20.0                ! maximal position in y coordinate
        z_min =   0.0                ! minimal position in z coordinate
        z_max =  12.0                ! maximal position in z coordinate
	nlevel = 3                   ! grid level factor to amplifiy grid resolution 
	delta_s = 0.25               ! integral step (the default value is 0.25)
/
```

## Other Note
This parallel version is based on FORTRAN OPENMP.

If one want to use N threads for the calculation, just change the value of parameter 'nthreads' in parameter file. If this parameter is defined as 0, then the max number of threads in the computer will be used as default. The default value in the program is 1.

The field line integral method is Runge-Kutta 4(5).


## Checking Accuracy
Based on the PFSS dipole field:

$$B_r=\frac{R_{\odot}^3}{r^3}(\frac{2Rs^3+r^3}{R_{\odot}^3+2Rs^3}+1)\cos(\theta)$$

$$B_{\theta}=\frac{R_{\odot}^3}{r^3}(\frac{Rs^3-r^3}{R_{\odot}^3+2Rs^3})\sin(\theta)$$

With the source surface assumption of $Rs$= 2.5 $R_{\odot}$:

$$B_r=\frac{4}{129}(\frac{125}{4r^3}+1)\cos(\theta)$$

$$B_{\theta}=\frac{4}{129r}(\frac{125}{8r^2}-r)\sin(\theta)$$

It is very clear that the $Q_{\perp}$ in the closed field region should be constantly 2, since the field line induced map is symmetric and uniform:

$$(R_{cal},\theta,\phi)\mapsto(R_{cal},\pi-\theta,\phi)$$

Let's look at the open field line, $\Phi$ component is zero, the field line equation is: $`\frac{dr}{B_r}=\frac{rd\theta}{B_{\theta}}`$, 
since the field components are radius-angle separatable, the equation can be convert to:

$$\frac{125/(8r^2)-r}{125/(4r^3)+1}dr=\frac{\cos(\theta)}{\sin(\theta)}d\theta$$

Then the LHS and RHS can be integral separately, lead to $`\ln(\sin(\theta_{Rs})/\sin(\theta_{Rcal}))=Const.`$, where $`Const.=\int\frac{125/(8*r^3)-1}{125/(4*r^2)+r}dr`$.
By using NIntegrate in Mathematica, it is 0.266657, when we start from $R_{cal}=1.01R_{\odot}$ and end at $Rs=2.5R_{\odot}$, where $R_{cal}$ is the spherical surface where we calculate the squashing factor $Q_{\perp}$.
Then $`\sin(\theta_{Rs})=e^{Const.}\sin(\theta_{Rcal})`$. We can obtain the separatrix layer locates on $R_{cal}$ surface, $\theta_{SL}=\arcsin(e^{-Const.})\approx 0.8724955346952092$.

By considering the $Q_{\perp}$ defined by Titov 2007 ApJ 660 863, from Eqs. (36)--(41).
The determinant, $\Delta_{\perp}$, of the map is $|\frac{B_{Rcal}}{B_{Rs}}|$.
The matrix $G^*_{\perp}$ is

```math
R_S^2\begin{bmatrix}
\sin^2{\theta_{Rs}} & 0\\
0 & 1
\end{bmatrix}
```

The matrix $G^{\perp}$ is
```math
R_{cal}^{-2}\begin{bmatrix}
\sin^{-2}(\theta_{Rcal}) & 0\\
0 & 1+\frac{B^2_{\theta,Rcal}}{B^2_{r,Rcal}}
\end{bmatrix}
```

The matrix $D$ is
```math
\begin{bmatrix}
1 & 0 \\
0 & e^{Const.}\frac{\cos(\theta_{Rcal})}{\cos(\theta_{Rs})}
\end{bmatrix}
```

Then the square of the matrix norm $`N_{\perp}^2=tr(D^{\top}G^*_{\perp}DG^{\perp})`$ is

$$N_{\perp}^2=\frac{R_S^2}{R_{cal}^2}e^{2Const.}(1+\frac{\cos^2(\theta_{Rcal})}{1-e^{2Const.}\sin^2(\theta_{Rcal})}(1+\frac{B^2_{\theta,Rcal}}{B^2_{r,Rcal}})) $$

Then $Q_{\perp}=N_{\perp}^2/\Delta_{\perp}$. The results are shown in the following figure:
<p align="center"><img src="https://raw.githubusercontent.com/Kai-E-Yang/QSL/master/fig/K-QSL_vs_Theory.png" /></p>

The upper-left panel reveals a close match between K-QSL results and theoretical values. Along the latitude, $Q_{\perp}$ exhibits a narrow layer where its value exceeds 10. While the separatrix layer should ideally display an infinite value, a denser mesh grid is needed to accurately capture this structure; this is illustrated in the bottom panel, which compares K-QSL calculations based on 200,000 points with theoretical calculations using $10^7$ points. The upper-right panel presents findings from varied delta_s values in K-QSL, suggesting that delta_s < 1 could give a better result in practice.

## TODO List
- non-uniform magnetic field (progress ... 100%).
- adaptive mesh refinement for capturing the fine structure (progress ... 90%).
- A new method will come soon (progress ... 50%).
