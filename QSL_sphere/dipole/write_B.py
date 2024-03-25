import numpy as np
import f90nml
import os

## Define current path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

## Read in necessary inputs from par file
nml = f90nml.read('par')
output_file_name = nml['filename_par']['bfieldname']
n_r = nml['cal_par']['dim_ra']
n_t = nml['cal_par']['dim_th']
n_p = nml['cal_par']['dim_ph']
ra_start = nml['cal_par']['ra_start']
th_start = nml['cal_par']['th_start']
ph_start = nml['cal_par']['ph_start']
ra_end = nml['cal_par']['ra_end']
th_end = nml['cal_par']['th_end']
ph_end = nml['cal_par']['ph_end']

# Create radius, theta, and phi arrays
radius = 1.5*np.arange(n_r)/(n_r-1)+1 # range 1 to 2.5
theta = np.pi*np.arange(n_t)/(n_t-1) # range 0 to pi
phi = 2*np.pi*np.arange(n_p)/(n_p-1) # range 0 to 2pi
# phi = np.pi*np.arange(n_p)/(n_p-1) # range 0 to 2pi
# Set up coords as a 3D array to feed into B
posir = np.zeros((n_r,n_t,n_p))
posit = np.zeros((n_r,n_t,n_p))
posip = np.zeros((n_r,n_t,n_p))

for i in range(0, n_r):
    for j in range(0, n_t):
        for k in range(0, n_p):
            posir[i,j,k] = radius[i]
            posit[i,j,k] = theta[j]
            posip[i,j,k] = phi[k]

## Hence define B
br = np.zeros((n_r,n_t,n_p))
bt = np.zeros((n_r,n_t,n_p))
bp = np.zeros((n_r,n_t,n_p))
            
br = 4/129*(125/(4*posir**3)+1)*np.cos(posit)
bt = 4/(129*posir)*(125/(8*posir**2)-posir)*np.sin(posit)

## Write B-field components to binary file
with open(output_file_name, 'wb') as file:
    file.write(radius.tobytes('F'))
    file.write(theta.tobytes('F'))
    file.write(phi.tobytes('F'))
    file.write(br.tobytes('F'))
    file.write(bt.tobytes('F'))
    file.write(bp.tobytes('F'))

