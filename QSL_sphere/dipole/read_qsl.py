import numpy as np
import f90nml
import matplotlib.pyplot as plt
import os

## Define current path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

nml    = f90nml.read('par')
infile = nml['filename_par']['outfilename']
rdim = nml['cal_par']['dim_ra']
tdim = nml['cal_par']['dim_th']
pdim = nml['cal_par']['dim_ph']
thstart = nml['cal_par']['th_start']
phstart = nml['cal_par']['ph_start']
rastart = nml['cal_par']['ra_start']
thend   = nml['cal_par']['th_end']
phend   = nml['cal_par']['ph_end']
raend   = nml['cal_par']['ra_end']
nlevel = nml['cal_par']['nlevel']

## Define radius, theta, phi arrays for B field calculation
# radius = 1.5*np.arange(rdim)/(rdim-1)+1 # range 1 to 2.5
# theta = np.pi*np.arange(tdim)/(tdim-1) # range 0 to pi
# phi = 2*np.pi*np.arange(pdim)/(pdim-1) # range 0 to 2pi

## Redefine dims
pdim = (pdim-1)*nlevel-1
tdim = (tdim-1)*nlevel-1
rdim = 1

f    = open(infile,'rb')
data = np.fromfile(f, dtype=np.double, count=rdim*tdim*pdim*3)
q    = np.reshape(data,(3,pdim,tdim,rdim))

## Plots of theta against phi
plt.figure(figsize=(5,5))
a1=plt.subplot(131)
plt.imshow(q[0,:,:,0],vmin=-5, vmax=5,origin='lower',cmap='bwr')
plt.title('sign(Br)log(Q) for $r=%s$' %rastart)
plt.xlabel('$\Theta$')
plt.ylabel('$\phi$')

a2=plt.subplot(132)
plt.imshow(q[1,:,:,0],vmin=0, vmax=5,origin='lower',cmap='jet')
plt.title('Field Line Length for $r=%s$' %rastart)
plt.xlabel('$\Theta$')
plt.ylabel('$\phi$')

a3=plt.subplot(133)
plt.imshow(q[2,:,:,0],vmin=0, vmax=5,origin='lower',cmap='hsv')
plt.title('Mask for $r=%s$' %rastart)
plt.xlabel('$\Theta$')
plt.ylabel('$\phi$')

plt.show()
