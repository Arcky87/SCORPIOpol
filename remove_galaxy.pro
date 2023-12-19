;remove galaxy
wdir='h:\red_data.pol\Sy1\Mkn79_151207\'
goto,cont
;формирование прямого изображения
;N=1024 & Nf=6
;files=+'maps\'+'s135813'+string(findgen(Nf)+3,format='(I2.2)')+'.fts'
;cube=fltarr(N,N,Nf)
;for k=0,Nf-1  do cube(*,*,k)=FLOAT(readfits(wdir+files(k),h)-1000)
;;определение смещения изображений
;xc=557  & yc=620 & wc=25
;w=fltarr(Nf)
;map=fltarr(2*wc+1,2*wc+1)
;for k=0,Nf-1 do begin
;map(*,*)=cube(xc-wc:xc+wc,yc-wc:yc+wc,k)
;yfit = MPFIT2DPEAK(map, A,/titl,/moffat );[, X, Y, /TILT ...] )
;xs=A(4)-wc  & ys=A(5)-wc & FWHM=sqrt(A(2)^2+a(3)^2)/sqrt(2)*0.357
;cube(*,*,k)=shift_image(cube(*,*,k),-xs,-ys)
;yfit = MPFIT2DPEAK(map, A,/titl,/moffat );[, X, Y, /TILT ...] )
;xs=A(4)-wc  & ys=A(5)-wc & FWHM=sqrt(A(2)^2+a(3)^2)/sqrt(2)*0.357
;print,xs,ys,FWHM,A(0),A(1)
;w(k)=FWHM
;endfor
;map=fltarr(N,N)
;;построение среднего изображения
;for kx=0,N-1 do begin
;for ky=0,N-1 do begin
;map(kx,ky)=median(cube(kx,ky,*))
;endfor & endfor
;writefits,wdir+'map.fts',map,h
;Nf=10
;file_flat='s135914'+string(findgen(Nf)+1,format='(I2.2)')+'.fts'
;file_flat=wdir+'flat-i\'+file_flat
;Nf=10
;cube_flat=fltarr(N,N,Nf)
;norm=fltarr(Nf)
;for k=0,Nf-1 do begin
;cube_flat(*,*,k)=FLOAT(readfits(file_flat(k) ,/silent)-1000)
;robomean,cube_flat(N/2-100:N/2+100,N/2-100:N/2+100,k),3,0.5,mean
;norm(k)=mean
;endfor
;flat=fltarr(N,N)
;robomean,norm,3,0.5,avg_norm &  norm=norm/avg_norm
;for kx=0,N-1 do begin
;for ky=0,N-1 do begin
;flat(kx,ky)=median(cube_flat(kx,ky,*)/norm)
;endfor & endfor
;robomean,flat(N/2-100:N/2+100,N/2-100:N/2),3,0.5,avg_flat & flat=flat/avg_flat
;window,0,xsize=N,ysize=N
;tv,255-bytscl(flat,0.9,1.1)
;writefits,wdir+'flat_map.fts',flat
;xc=515 & yc=516 & NS=1720
;map=map/flat-NS
;writefits,wdir+'map.fts',map(xc-N/4:xc+N/4-1,yc-N/4:yc+N/4-1),h
;cont:
;чтение спектра
obj_cube=readfits(wdir+'obj-sky.fts',h)
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2') & Np=sxpar(h,'NAXIS3') & wx=100
Nt=sxpar(h,'NAXIS5') & Nexp=intarr(Nt) & gain=sxpar(h,'GAIN')
for k=0,Nt-1 do Nexp(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))

obj=fltarr(Nx,Ny,Np,1,Nt)
;нормировка спектров

		for t=0,Nt-1 do begin
	for p=0,Np-1 do begin
;определение смещенией спектров вдоль щели

vector=fltarr(Ny,Nexp(t)) & dx=fltarr(Nexp(t)) & norm=fltarr(Nexp(t))

for k=0,Nexp(t)-1 do begin
vector(*,k)=total(obj_cube(Nx/2-wx:Nx/2+wx,*,p,k,t),1)
norm(k)=total(vector(*,k))
endfor
robomean,norm,3,0.5,avg_norm  & norm=norm/avg_norm
print,norm,format='('+string(Nexp(t)+1,format='(I2)')+'F10.3)'
;смещения и нормировка

for k=0,Nexp(t)-1 do obj_cube(*,*,p,k,t)=obj_cube(*,*,p,k,t)/norm(k)
	endfor
		endfor

	print,systime()
	   for t=0,Nt-1 do begin
	for p=0,Np-1 do begin

for kx=0,Nx-1 do begin
for ky=0,Ny-1 do begin
obj(kx,ky,p,0,t)=median(obj_cube(kx,ky,p,0:Nexp(t)-1,t))*Nexp(t)*gain
endfor
endfor

	endfor
	endfor
	print,systime()
writefits,wdir+'avg_obj.fts',obj,h

;mapa=fltarr(Nx,Ny)
;p=0 & t=2
;mapa(*,*)=obj(*,*,p,t)
;window,2,xsize=1000 & ysize=500
;tv,255-bytscl(congrid(mapa,1000,500))
;cont:


map=readfits(wdir+'map.fts')
window,0,xsize=N,ysize=N
!P.multi=[0,1,1]
xc=255 & w=3 & w2=2*w
gal=total(map(xc-w+w2:xc+w+w2,*),1)/(2*w+1)
nuc=total(map(xc-w:xc+w,*),1)/(2*w+1)
plot,nuc-gal,xst=1;,/ylog;,yrange=[1,1e5]
oplot,gal,color=3e5
oplot,[1,1]*275,[-1,1]*2e4,linestyle=2
w=10
print,total(gal(xc-w:xc+w)),total(nuc(xc-w:xc+w))
;исправление атмосферной дисперсии
;cont:    ;**********экстракция спектров*********

obj=readfits(wdir+'avg_obj.fts',h)
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2') & Np=sxpar(h,'NAXIS3') & wx=100
Nt=sxpar(h,'NAXIS5')

y=findgen(Ny)  & x=findgen(Nx)

spectra=fltarr(Nx,Np,Nt)  & wy=5
for kx=0,Nx-1 do spectra(kx,*,*)=TOTAL(obj(kx,Ny/2-wy:Ny/2+wy,*,0,*),2)
window,2,xsize=900,ysize=1100
!p.multi=[0,1,4]
t=1
;for k=0,Np-1 do plot,spectra(*,k,t),xst=1
;вычисление параметров Стокса звезд

Q_zero=(spectra(*,1,1)-spectra(*,0,1))/(spectra(*,1,2)+spectra(*,0,2))
Q_null=LOWESS(findgen(Nx),Q_zero,Nx/4,2,1)
U_zero=(spectra(*,3,1)-spectra(*,2,1))/(spectra(*,3,2)+spectra(*,2,2))
U_null=LOWESS(findgen(Nx),U_zero,Nx/4,2,1)
Q_star=(spectra(*,1,2)-spectra(*,0,2))/(spectra(*,1,2)+spectra(*,0,2));-Q_null
U_star=(spectra(*,3,2)-spectra(*,2,2))/(spectra(*,3,2)+spectra(*,2,2));-U_null

plot,Q_zero,xst=1,charsize=2
plot,U_zero,xst=1,charsize=2
plot,SQRT(Q_zero^2+U_zero^2),xst=1,charsize=2
FI=calc_atan(Q_zero,U_zero)/2
plot,FI,xst=1,charsize=2
;plot,Q_star,xst=1,charsize=2
;plot,U_star,xst=1,charsize=2
;plot,SQRT(Q_star^2+U_star^2),xst=1,charsize=2
;FI=calc_atan(Q_star,U_star)/2
;plot,FI,xst=1,charsize=2
;вычисление параметров Стокса объекта
Q_obj=(spectra(*,1,0)-spectra(*,0,0))/(spectra(*,1,0)+spectra(*,0,0)) ;-Q_null-0.1
U_obj=(spectra(*,3,0)-spectra(*,2,0))/(spectra(*,3,0)+spectra(*,2,0))  ;-U_null
window,1,xsize=900,ysize=1100
!p.multi=[0,1,4]
plot,Q_obj,xst=1,charsize=2
plot,U_obj,xst=1,charsize=2
plot,SQRT(Q_obj^2+U_obj^2),xst=1,charsize=2
FI=calc_atan(Q_obj,U_obj)/2
plot,FI,xst=1,charsize=2,yrange=[150,170]
;p=0 & t=0
;mapa(*,*)=obj(*,*,p,t)
;window,2,xsize=1000 & ysize=500
;tv,255-bytscl(congrid(mapa,1000,500),0,max(mapa)/20)

cont:
;чтение изображения
map=readfits(wdir+'map.fts')
N=512 & y=findgen(N)
;Window,0;,xsize=512,ysize=512
;tv,255-bytscl(map,0,500)
window,2,xsize=600,ysize=900
!P.multi=[0,1,2]
xc=255 & wx=2 & dx=2*wx
slice=total(map(xc-wx:xc+wx,*),1)
slice_gal=(total(map(xc-wx-dx:xc+wx-dx,*),1)+total(map(xc-wx+dx:xc+wx+dx,*),1))/2
slice_nuc=slice-slice_gal
plot,slice,xst=1,xrange=[xc-80,xc+80],yst=1,yrange=[-1,15]*1e4
oplot,slice_gal,color=3e5
oplot,slice_nuc,color=1e5

profile=MPFITPEAK(y,slice_nuc, A,/MOFFAT) & print,A(0:2)
oplot,profile,linestyle=1
oplot,slice_nuc-profile
yc=256  & wy=10
oplot,[1,1]*(yc-wy),[-1,1]*1e6,linestyle=2
oplot,[1,1]*(yc+wy),[-1,1]*1e6,linestyle=2
print,total(slice_gal(yc-wy:yc+wy)),total(slice_nuc(yc-wy:yc+wy))
;чтение среднего спектра
 obj=readfits(wdir+'avg_obj.fts',h)
a=size(avg_obj) & Nx=a(1) & Ny=a(2) & Npol=a(3) & Nt=a(4)
tot_obj=total(avg_obj(*,*,*,0),3)

kx=Nx/2
slice_obj=tot_obj(kx,*)-min(tot_obj(kx,*))
;slice_obj=congrid(slice_obj(yc-Ny/4:yc+Ny/4-1),Ny)
plot,slice_obj,xst=1,yrange=[-0.6,10]*1e3,yst=1
galaxy=congrid(slice_gal(yc-Ny/4:yc+Ny/4-1),Ny)
cross=congrid(slice(yc-Ny/4:yc+Ny/4-1),Ny) &
;формирование PSF
y=findgen(Ny)
print,Ny
PSF=gaussian(y,[1,Ny/2,4]) & PSF=PSF/total(PSF)

RES_slice=ABS(shift(FFT(FFT(PSF,1)*FFT(cross,1),-1),Ny/2))/13
RES_galaxy=ABS(shift(FFT(FFT(PSF,1)*FFT(galaxy,1),-1),Ny/2))/13
oplot, RES_slice,color=3e5
oplot,RES_galaxy,color=3e5
RES_nuc=slice_obj-RES_galaxy(*,0)


fit=MPFITPEAK(y,RES_nuc, A,/moffat) & print,A(2)
oplot, fit,color=3e5
oplot, RES_nuc,linestyle=1
yc=Ny/2-2
oplot,[1,1]*(yc-2*wy),[-1,1]*1e6,linestyle=2
oplot,[1,1]*(yc+2*wy),[-1,1]*1e6,linestyle=2

print,total(RES_galaxy(yc-2*wy:yc+2*wy,0)),total(RES_nuc(yc-2*wy:yc+2*wy))
window,0,xsize=1500
!p.multi=[0,1,1]
print,total(RES_galaxy(yc-2*wy:yc+2*wy,0)),total(RES_galaxy(yc-6*wy:yc-2*wy,0))+total(RES_galaxy(yc+2*wy:yc+6*wy,0))
plot,total(total(obj(*,Ny/2-2*wy:Ny/2+2*wy,*,0),3),2),xst=1
oplot,(total(total(obj(*,Ny/2-6*wy:Ny/2-2*wy,*,0),3),2)+total(total(obj(*,Ny/2+2*wy:Ny/2+6*wy,*,0),3),2))*2.8

end