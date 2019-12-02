import numpy as np
import f90nml
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

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

fig, axs = plt.subplots(1, 3, figsize=(15, 5))
fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.05)

im1 = axs[0].imshow(q[0,:,0,:], interpolation='nearest', vmin=-5, vmax=5,origin='lower',cmap='bwr')
fig.colorbar(im1, ax=axs[0])
axs[0].set_title("sign(B$_z$)log$_{10}$(Q)")

im2 = axs[1].imshow(q[1,:,0,:], interpolation='nearest', vmin=0, vmax=80,origin='lower',cmap='jet')
fig.colorbar(im2, ax=axs[1])
axs[1].set_title("Field Line Length")


colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0)]
n_bins = 4
cmap_name = 'my_list'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)


im2 = axs[2].imshow(q[2,:,0,:], interpolation='nearest', vmin=0, vmax=5,origin='lower',cmap=cm)
fig.colorbar(im2, ax=axs[2])
axs[2].set_title("Mask")

# # a1=plt.subplot(131)
# img1=plt.imshow(q[0,:,0,:],vmin=-5, vmax=5,origin='lower',cmap='bwr')
# plt.title('sign(B$_z$)log$_{10}$(Q)')


# # a2=plt.subplot(132)
# plt.imshow(q[1,:,0,:],vmin=0, vmax=80,origin='lower',cmap='jet')
# plt.title('Field Line Length')

# # a3=plt.subplot(133)
# plt.imshow(q[2,:,0,:],vmin=0, vmax=5,origin='lower')
# plt.title('Mask')

plt.show()

