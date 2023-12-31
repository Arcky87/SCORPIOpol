;NS substraction
wdir='e:\sbs1419+538_190216\'
cube=readfits(wdir+'obj_lin.fts',h)
cube_corr=cube

o=0
for e=0, 15 do begin
	for p=0,3 do begin
		print, p, e, o
		ima=cube(*,*,p,e,o)
		Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2')
		;definition substraction region
		wx=15
		V=total(ima(Nx/2-wx:Nx/2+wx,*),1)/(2*wx+1)
		window,2
		y=findgen(Ny)
		!P.multi=[0,1,1]
		plot,V,xst=1
		robomean,V,3,0.5,mean,rms
		Ndeg=2 & tresh=3
	   reg=[0,30,50,Ny-1]
		reg1=(y ge reg(0) and y le reg(1))
		reg2=(y ge reg(2))and(y le reg(3))
		R=where(reg1 or reg2,cc) & print, cc
		oplot,y(R),V(R),psym=6,symsize=0.25, color=1e5
		for i=0,Nx-1 do begin
			V(*)=ima(i,*)
			f=goodpoly(y(R),V(R),Ndeg,tresh,yfit)
			sky=0 & for j=0,Ndeg do sky=sky+f(j)*y^j
			ima(i,*)=V-sky
		endfor
		cube_corr(*,*,p,e,o)=ima
	endfor
endfor

for o=1, 2  do begin
for e=0, 4 do begin
	for p=0,3 do begin
		ima=cube(*,*,p,e,o)
		Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2')
		;definition substraction region
		wx=20
		V=total(ima(Nx/2-wx:Nx/2+wx,*),1)/(2*wx+1)
		window,2
		y=findgen(Ny)
		!P.multi=[0,1,1]
		plot,V,xst=1
		robomean,V,3,0.5,mean,rms
		Ndeg=2 & tresh=3
		reg=[0,20,60,Ny-1]
		reg1=(y ge reg(0) and y le reg(1))
		reg2=(y ge reg(2))and(y le reg(3))
		R=where(reg1 or reg2,cc) & print, cc
		oplot,y(R),V(R),psym=6,symsize=0.25, color=1e5
		for i=0,Nx-1 do begin
			V(*)=ima(i,*)
			f=goodpoly(y(R),V(R),Ndeg,tresh)
			sky=0 & for j=0,Ndeg do sky=sky+f(j)*y^j
			ima(i,*)=V-sky
		endfor
		cube_corr(*,*,p,e,o)=ima
	endfor
endfor
endfor

	writefits,wdir+'obj-sky.fts',cube_corr,h

window,3
robomean,cube_corr(*,*,0,0,0),3,0.5,avg_ima,rms_ima
map=255-bytscl(cube_corr(*,*,0,0,0),avg_ima-rms_ima*3,avg_ima+5*rms_ima)
TV,map

end