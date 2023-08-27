## This code is written to replace bxyz.pro
import numpy as np
import f90nml
import os

## Define current path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

## Read in necessary inputs from par file
nml = f90nml.read('par')
output_file_name = nml['filename_par']['bfieldname']
dimx = nml['cal_par']['dimx']
dimy = nml['cal_par']['dimy']
dimz = nml['cal_par']['dimz']
minx = nml['cal_par']['x_min']
maxx = nml['cal_par']['x_max']
miny = nml['cal_par']['y_min']
maxy = nml['cal_par']['y_max']
minz = nml['cal_par']['z_min']
maxz = nml['cal_par']['z_max']

## Establish Cartesian grid points
x = np.linspace(minx, maxx, dimx)
y = np.linspace(miny, maxy, dimy)
z = np.linspace(minz, maxz, dimz)

## Set up coords as a 3D array to feed into B
posix = np.zeros((dimx,dimy,dimz))
posiy = np.zeros((dimx,dimy,dimz))
posiz = np.zeros((dimx,dimy,dimz))

for i in range(0, dimx):
    for j in range(0, dimy):
        for k in range(0, dimz):
            posix[i,j,k] = x[i]
            posiy[i,j,k] = y[j]
            posiz[i,j,k] = z[k]

## Hence define B
bx = np.zeros((dimx,dimy,dimz))
by = np.zeros((dimx,dimy,dimz))
bz = np.zeros((dimx,dimy,dimz))

k = 0.5
bx = k*posix
by = (1-k)*posiy
bz = -1.0*posiz

## Write bx, by, bz to binary file that can be read by Fortran
with open(output_file_name, 'wb') as file:
    file.write(bx.tobytes('F'))
    file.write(by.tobytes('F'))
    file.write(bz.tobytes('F'))
