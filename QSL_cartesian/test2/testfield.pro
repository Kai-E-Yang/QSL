pro testfield
;give a test field for testing the null scaning program
;we use the magnetic charges model, 4 charges with their magnitude of unit: 1(0,0,0),-1(0,1,0),-1(cos(7pi/6),sin(7pi/6),0),-1(cos(11pi/6),sin(11pi/6),0)
;domain dimension is (2*2*2)
;per pixel is 0.02
pixel = 0.02
bx = dblarr(50,50,50)
by = dblarr(50,50,50)
bz = dblarr(50,50,50)
origin=[25.5,25.5,0]
c1 = [0.,0.,-0.05]*0.4;positive
c2 = [0.,1.,-0.05]*0.4;negative
c3 = [cos(7*!dpi/6),sin(7*!dpi/6),-0.05]*0.4;negative
c4 = [cos(11*!dpi/6),sin(11*!dpi/6),-0.05]*0.4;negative
for i=0,49 do begin
	for j=0,49 do begin
		for k=0,49 do begin
			r = (1.*[i,j,k]-origin)*pixel
			r1 = r - c1
			r2 = r - c2
			r3 = r - c3
			r4 = r - c4
			b1 = 1.*r1/(sqrt(r1[0]^2 + r1[1]^2 + r1[2]^2)^3)
			b2 = -1.*r2/(sqrt(r2[0]^2 + r2[1]^2 + r2[2]^2)^3)
			b3 = -1.*r3/(sqrt(r3[0]^2 + r3[1]^2 + r3[2]^2)^3)
			b4 = -1.*r4/(sqrt(r4[0]^2 + r4[1]^2 + r4[2]^2)^3)
			bx[i,j,k] = b1[0] + b2[0] + b3[0] + b4[0]
			by[i,j,k] = b1[1] + b2[1] + b3[1] + b4[1]
			bz[i,j,k] = b1[2] + b2[2] + b3[2] + b4[2]
		endfor
	endfor
endfor

openw,1,'bxyz.binary'
writeu,1,bx
writeu,1,by
writeu,1,bz
close,1
end
