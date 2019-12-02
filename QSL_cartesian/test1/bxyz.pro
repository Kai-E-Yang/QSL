nx = 32
ny = 32
nz = 32

posix = dblarr(nx+1,ny+1,nz+1)
posiy = dblarr(nx+1,ny+1,nz+1)
posiz = dblarr(nx+1,ny+1,nz+1)

bx = dblarr(nx+1,ny+1,nz+1)
by = dblarr(nx+1,ny+1,nz+1)
bz = dblarr(nx+1,ny+1,nz+1)


for k=0,32 do begin
for j=0,32 do begin
for i=0,32 do begin
   posix[i,j,k]=0.1*i - 0.1*nx*0.5;+0.05
   posiy[i,j,k]=0.1*j - 0.1*ny*0.5;+0.05
   posiz[i,j,k]=0.1*k - 0.1*nz*0.5;+0.05
endfor
endfor
endfor

k = 0.5

xc = 1.5d0
zc = 1.5d0
bx = k*posix
by = (1-k)*posiy
bz = -1.0*posiz

z0 = posiz[0,0,15]
z1 = posiz[0,0,0]
x0 = reform(posix[*,*,15])
y0 = reform(posiy[*,*,15])

bbn = sqrt( 0.25*(z0/z1)*(x0^2+y0^2) + z1^2 )
bb0 = sqrt( 0.25*(x0^2+y0^2) + z0^2 )

qb = 2*(z0^2/z1^2)/(bb0/bbn)

; vec_2vtk,bx,by,bz,'bxyz.vtk'
openw,1,'bxyz.binary'
writeu,1,bx
writeu,1,by
writeu,1,bz
close,1

zz = -0.5

end