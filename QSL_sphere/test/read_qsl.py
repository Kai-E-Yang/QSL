import numpy as np
import f90nml
import matplotlib.pyplot as plt

nml    = f90nml.read('par')
infile = nml['filename_par']['outfilename']
thstart = nml['cal_par']['th_start']
phstart = nml['cal_par']['ph_start']
rastart = nml['cal_par']['ra_start']
thend   = nml['cal_par']['th_end']
phend   = nml['cal_par']['ph_end']
raend   = nml['cal_par']['ra_end']
nlevel = nml['cal_par']['nlevel']

# rdim   = np.int(nlevel*(raend-rastart) + 1)
# tdim   = np.int(nlevel*(thend-thstart) + 1)
# pdim   = np.int(nlevel*(phend-phstart) + 1)

rdim=1
tdim=179
pdim=360

f    = open(infile,'rb')
data = np.fromfile(f, dtype=np.double, count=rdim*tdim*pdim*3)
q    = np.reshape(data,(3,pdim,tdim,rdim))

a1=plt.subplot(131)
plt.imshow(q[0,:,:,0],vmin=-5, vmax=5,origin='lower',cmap='bwr')
plt.title('sign(Bz)log(Q)')

a2=plt.subplot(132)
plt.imshow(q[1,:,:,0],vmin=0, vmax=5,origin='lower',cmap='jet')
plt.title('Field Line Length')

a3=plt.subplot(133)
plt.imshow(q[2,:,:,0],vmin=0, vmax=5,origin='lower',cmap='hsv')
plt.title('Mask')

plt.show()

