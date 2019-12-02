import numpy as np
import f90nml
import matplotlib.pyplot as plt

nml    = f90nml.read('par')
infile = nml['filename_par']['outfilename']
xstart = nml['cal_par']['x_start']
ystart = nml['cal_par']['y_start']
zstart = nml['cal_par']['z_start']
xend   = nml['cal_par']['x_end']
yend   = nml['cal_par']['y_end']
zend   = nml['cal_par']['z_end']
nlevel = nml['cal_par']['nlevel']

xdim   = nlevel*(xend-xstart) + 1
ydim   = nlevel*(yend-ystart) + 1
zdim   = nlevel*(zend-zstart) + 1

f    = open(infile,'rb')
data = np.fromfile(f, dtype=np.double, count=xdim*ydim*zdim*3)
q    = np.reshape(data,(3,zdim,ydim,xdim))

a1=plt.subplot(131)
plt.imshow(q[0,0,:,:],vmin=-5, vmax=5,origin='lower')
plt.title('sign(Bz)log(Q)')

a2=plt.subplot(132)
plt.imshow(q[1,0,:,:],vmin=0, vmax=1e2,origin='lower')
plt.title('Field Line Length')

a3=plt.subplot(133)
plt.imshow(q[2,0,:,:],vmin=0, vmax=5,origin='lower')
plt.title('Mask')

plt.show()

