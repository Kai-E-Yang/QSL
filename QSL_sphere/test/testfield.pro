; 3d null point
n_r=41
n_t=180
n_p=360
n_q=5

q1=1
q2=1
q3=1
q4=1
q5=-3

qc=[q1,q2,q3,q4,q5]
dt=0.1
dp=0.1
rc=[0.9,0.9,0.9,0.9,0.9]
tc=[0.5-dt,0.5,0.5+dt,0.5,0.5]*!dpi
pc=[1.0,1.0+dp,1,1-dp,1.0]*!dpi

er=dblarr(n_r,n_t,n_p,3)
ep=dblarr(n_r,n_t,n_p,3)
et=dblarr(n_r,n_t,n_p,3)

bxyz=dblarr(n_r,n_t,n_p,3)


posic=dblarr(5,3)
posi=dblarr(n_r,n_t,n_p,3)
rmap=dblarr(n_r,n_t,n_p)
tmap=dblarr(n_r,n_t,n_p)
pmap=dblarr(n_r,n_t,n_p)

posic[*,0]=rc*sin(tc)*cos(pc)
posic[*,1]=rc*sin(tc)*sin(pc)
posic[*,2]=rc*cos(tc)

radius=2*dindgen(n_r)/(n_r-1)+1
theta=!pi*dindgen(n_t)/(n_t-1)
phi=2*!pi*dindgen(n_p)/(n_p-1)

for i=0,n_r-1 do rmap[i,*,*]=radius[i]
for j=0,n_t-1 do tmap[*,j,*]=theta[j]
for k=0,n_p-1 do pmap[*,*,k]=phi[k]

posi[*,*,*,0]=rmap*sin(tmap)*cos(pmap)
posi[*,*,*,1]=rmap*sin(tmap)*sin(pmap)
posi[*,*,*,2]=rmap*cos(tmap)

for i=0,n_q-1 do begin
	dr=sqrt((posi[*,*,*,0]-posic[i,0])^2+(posi[*,*,*,1]-posic[i,1])^2+(posi[*,*,*,2]-posic[i,2])^2)
	for idir=0,2 do begin
		bxyz[*,*,*,idir]=bxyz[*,*,*,idir]+qc[i]*(posi[*,*,*,idir]-posic[i,idir])/dr^3
	endfor
endfor


posir=sqrt(total(posi^2,4))
for idir=0,2 do er[*,*,*,idir]=posi[*,*,*,idir]/posir

ep[*,*,*,0]=-1*posi[*,*,*,1]
ep[*,*,*,1]=posi[*,*,*,0]
ep[*,*,*,2]=0.0
tmp=sqrt(total(ep^2,4))
for idir=0,2 do ep[*,*,*,idir]=ep[*,*,*,idir]/tmp

et[*,*,*,0]=ep[*,*,*,1]*er[*,*,*,2]-ep[*,*,*,2]*er[*,*,*,1]
et[*,*,*,1]=ep[*,*,*,2]*er[*,*,*,0]-ep[*,*,*,0]*er[*,*,*,2]
et[*,*,*,2]=ep[*,*,*,0]*er[*,*,*,1]-ep[*,*,*,0]*er[*,*,*,1]

br=total(bxyz*er,4)
bt=total(bxyz*et,4)
bp=total(bxyz*ep,4)

bt[*,0,*]=0.0d0
bp[*,0,*]=0.0d0

bt[*,n_t-1,*]=0.0d0
bp[*,n_t-1,*]=0.0d0

openw,1,'bxyz_test.binary'
writeu,1,radius
writeu,1,theta
writeu,1,phi
writeu,1,br
writeu,1,bt
writeu,1,bp
close,1
end
