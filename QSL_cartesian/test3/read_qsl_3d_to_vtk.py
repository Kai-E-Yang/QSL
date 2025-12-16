from pyevtk.hl import gridToVTK
import numpy as np
import f90nml

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

x3d = np.linspace(xstart,xend,xdim)
y3d = np.linspace(ystart,yend,ydim)
z3d = np.linspace(zstart,zend,zdim)


x3d_contiguous = np.ascontiguousarray(x3d.astype(np.float32))
y3d_contiguous = np.ascontiguousarray(y3d.astype(np.float32))
z3d_contiguous = np.ascontiguousarray(z3d.astype(np.float32))

qsl = np.ascontiguousarray(q[0].transpose(2, 1, 0).astype(np.float32))
length = np.ascontiguousarray(q[1].transpose(2, 1, 0).astype(np.float32))
marker = np.ascontiguousarray(q[2].transpose(2, 1, 0).astype(np.float32))

vtrname=f'qsl_3d'
gridToVTK(vtrname, x3d_contiguous, y3d_contiguous, z3d_contiguous, 
          pointData={"q": qsl,"length": length, "marker":marker})