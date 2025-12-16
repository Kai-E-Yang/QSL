# pro testfield
# give a test field for testing the null scaning program
# we use the magnetic charges model, 4 charges with their magnitude of unit: 1(0,0,0),-1(0,1,0),-1(cos(7pi/6),sin(7pi/6),0),-1(cos(11pi/6),sin(11pi/6),0)
# domain dimension is (2*2*2)
# per pixel is 0.02
import numpy as np
pixel = 0.02
bx = np.zeros((50,50,50))
by = np.zeros((50,50,50))
bz = np.zeros((50,50,50))
origin=np.array([25.5,25.5,0])
c1 = np.array([0.,0.,-0.05])*0.4#positive
c2 = np.array([0.,1.,-0.05])*0.4#negative
c3 = np.array([np.cos(7*np.pi/6),np.sin(7*np.pi/6),-0.05])*0.4#negative
c4 = np.array([np.cos(11*np.pi/6),np.sin(11*np.pi/6),-0.05])*0.4#negative

for ix in range(50):
	for iy in range(50):
		for iz in range(50):
			r = (1.*np.array([ix,iy,iz])-origin)*pixel
			r1 = r - c1
			r2 = r - c2
			r3 = r - c3
			r4 = r - c4
			b1 = 1.*r1/(sqrt(r1[0]**2 + r1[1]**2 + r1[2]**2)**3)
			b2 = -1.*r2/(sqrt(r2[0]**2 + r2[1]**2 + r2[2]**2)**3)
			b3 = -1.*r3/(sqrt(r3[0]**2 + r3[1]**2 + r3[2]**2)**3)
			b4 = -1.*r4/(sqrt(r4[0]**2 + r4[1]**2 + r4[2]**2)**3)
			bx[i,j,k] = b1[0] + b2[0] + b3[0] + b4[0]
			by[i,j,k] = b1[1] + b2[1] + b3[1] + b4[1]
			bz[i,j,k] = b1[2] + b2[2] + b3[2] + b4[2]

## Write B-field components to binary file
output_file_name=f'bxyz_3d.binary'
with open(output_file_name, 'wb') as file:
    file.write(bx.tobytes('F'))
    file.write(by.tobytes('F'))
    file.write(bz.tobytes('F'))

