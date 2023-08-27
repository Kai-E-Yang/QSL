import numpy as np
import f90nml
import os

## Define current path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

## Read in necessary inputs from par file
nml = f90nml.read('par')
infile = nml['filename_par']['outfilename']
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

## Define charges and charge vectors
n_q = 5
q1 = 1
q2 = 1
q3 = 1
q4 = 1
q5 = -3
qc = [q1, q2, q3, q4, q5] 

## Define r, theta, phi charge coordinates
dt = 0.1
dp = 0.1
rc = [0.9, 0.9, 0.9, 0.9, 0.9]
tc = (np.array([0.4, 0.5, 0.6, 0.5, 0.5]))*np.pi
pc = (np.array([1.0, 1.1, 1.0, 0.9, 1.0]))*np.pi

## Initialize arrays
er = np.zeros((n_r, n_t, n_p, 3))
ep = np.zeros((n_r, n_t, n_p, 3))
et = np.zeros((n_r, n_t, n_p, 3))
bxyz = np.zeros((n_r, n_t, n_p, 3))
posic = np.zeros((5, 3))
posi = np.zeros((n_r, n_t, n_p, 3))
rmap = np.zeros((n_r, n_t, n_p))
tmap = np.zeros((n_r, n_t, n_p))
pmap = np.zeros((n_r, n_t, n_p))

## Calculate charge positions in Cartesian coordinates
posic[:, 0] = rc * np.sin(tc) * np.cos(pc)
posic[:, 1] = rc * np.sin(tc) * np.sin(pc)
posic[:, 2] = rc * np.cos(tc)

## Create radius, theta, and phi arrays
radius = 2*np.arange(n_r)/(n_r-1)+1  # range 1 to 3
theta = np.pi*np.arange(n_t)/(n_t-1) # range 0 to pi
phi = 2*np.pi*np.arange(n_p)/(n_p-1) # range 0 to 2pi

## Populate rmap, tmap, and pmap arrays
for i in range(n_r):
    rmap[i, :, :] = radius[i]
for j in range(n_t):
    tmap[:, j, :] = theta[j]
for k in range(n_p):
    pmap[:, :, k] = phi[k]

## Calculate x,y,z grid Cartesian coordinates
posi[:, :, :, 0] = rmap * np.sin(tmap) * np.cos(pmap)
posi[:, :, :, 1] = rmap * np.sin(tmap) * np.sin(pmap)
posi[:, :, :, 2] = rmap * np.cos(tmap)

## Calculate Cartesian B-field components
for i in range(n_q):
    dr = np.sqrt((posi[:, :, :, 0] - posic[i, 0]) ** 2 +(posi[:, :, :, 1] - posic[i, 1]) ** 2 +(posi[:, :, :, 2] - posic[i, 2]) ** 2) # dr is distance between points "posi" and the q charge "posic"
    for idir in range(3):
        bxyz[:, :, :, idir] += qc[i] * (posi[:, :, :, idir] - posic[i, idir]) / (dr**3)

## Calculate unit vectors in radial direction
posir = np.sqrt(np.sum(posi**2, axis=3))
for idir in range(3):
    er[:, :, :, idir] = posi[:, :, :, idir] / posir

## Calculate unit vectors in phi direction
ep[:, :, :, 0] = -1 * posi[:, :, :, 1]
ep[:, :, :, 1] = posi[:, :, :, 0]
ep[:, :, :, 2] = 0.0

epsilon = 1e-10 # to avoid division by zero
tmp = np.sqrt(np.sum(ep**2, axis=3))
for idir in range(3):
    ep[:, :, :, idir] = ep[:, :, :, idir] / (tmp+epsilon)

## Calculate unit vectors in theta direction (cross product of unit vectors in r and phi)
et[:, :, :, 0] = ep[:, :, :, 1] * er[:, :, :, 2] - ep[:, :, :, 2] * er[:, :, :, 1]
et[:, :, :, 1] = ep[:, :, :, 2] * er[:, :, :, 0] - ep[:, :, :, 0] * er[:, :, :, 2]
et[:, :, :, 2] = ep[:, :, :, 0] * er[:, :, :, 1] - ep[:, :, :, 1] * er[:, :, :, 0] # correct cross product formula
# et[:, :, :, 2] = ep[:, :, :, 0] * er[:, :, :, 1] - ep[:, :, :, 0] * er[:, :, :, 1] # incorrect cross product formula

## Calculate B-field components in spherical coordinates
br = np.sum(bxyz * er, axis=3)
bt = np.sum(bxyz * et, axis=3)
bp = np.sum(bxyz * ep, axis=3)

## Set B-field components at poles to zero
bt[:, 0, :] = 0.0
bp[:, 0, :] = 0.0
bt[:, -1, :] = 0.0
bp[:, -1, :] = 0.0

## Write B-field components to binary file; 'F'=Fortran contiguous
with open(output_file_name, 'wb') as file:
    file.write(radius.tobytes('F'))
    file.write(theta.tobytes('F'))
    file.write(phi.tobytes('F'))
    file.write(br.tobytes('F'))
    file.write(bt.tobytes('F'))
    file.write(bp.tobytes('F'))

