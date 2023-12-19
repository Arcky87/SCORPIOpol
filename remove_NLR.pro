function  angle_calculation,Q,U
teta_0=0
Teta=ATAN(U/Q)*180/!PI
if Q gt 0 and U ge 0 then TETA_0=0
if Q gt 0 and U lt 0 then TETA_0=360
if Q lt 0 then TETA_0=180
if Q eq 0 and U gt 0 then TETA_0=90
if Q eq 0 and U lt 0 then TETA_0=270

return, (TETA+TETA_0)/2
END

;analyz broad line
function remove_NLR,spectra,REP=rep,WIDTH=width
if not(keyword_set(rep)) then rep=4
if not(keyword_set(width))then width=0.1
N=N_elements(spectra)
for k=0,rep-1 do begin
cont=LOWESS(findgen(N),spectra,N*width,3,3)
robomean,spectra-cont,3,0.5,avg_cont,rms_cont
R=where(spectra-cont gt -rms_cont, ind)
if ind gt 0 then spectra(R)=cont(R)
endfor
return,spectra
end
;*********************************
set_plot,'WIN'
dir='h:\red_data.pol\Sy1\Mkn744_160307\\
;goto,cont
;чтение плоского поля
flat=readfits(dir+'avg_flat.fts')
obj=readfits(dir+'obj_lin.fts',h,/silent)

;
Nz=sxpar(h,'NAXIS5')
PA=fltarr(Nz) & for j=0,Nz-1 do PA(j)=sxpar(h,'PA'+string(j+1,format='(I1)'))
name=strarr(Nz) & for j=0,Nz-1 do name(j)=sxpar(h,'NAME'+string(j+1,format='(I1)'))
Nexp=intarr(Nz) & for j=0,Nz-1 do Nexp(j)=sxpar(h,'NUMEXP'+string(j+1,format='(I1)'))
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2') & Npol=4
y=findgen(Ny)  & d=50  & V=fltarr(Ny) & x=findgen(Nx)
; исправление плоского поля
bias=0
	for k=0,Nz-1 do begin
for j=0,Nexp(k)-1 do obj(*,*,*,j,k)=(obj(*,*,*,j,k)+bias)/flat
	endfor
;изменение нумерации каналов поляризации
num=[3,2,1,0]
tmp=obj
obj=tmp(*,*,num,*,*)

wave=findgen(Nx)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
z=sxpar(h,'Z') & print,z & z=z
Ha=6562.8
R=where(wave ge 6000 and wave le 7100,N)

wave=wave(R)
;обработка нулевого стандарта
;вычитание неба
		for e=0,Nexp(1)-1 do begin
	for p=0,Npol-1 do begin
for kx=0,Nx-1 do begin
V(*)=obj(kx,*,p,e,1)
f=goodpoly([y(0:d),y(Ny-d:Ny-1)],[V(0:d),V(Ny-d:Ny-1)],2,2)
obj(kx,*,p,e,1)=obj(kx,*,p,e,1)-f(0)-f(1)*y-f(2)*Y^2
endfor
	endfor
		endfor
yc=162  & wy=20
star=TOTAL(obj(*,*,*,*,1),4)
star=total(star(*,yc-wy:yc+wy,*),2)
for j=0,Npol-1 do star(*,j)=LOWESS(x,star(*,j),Nx/2,2,2)
window,0
!P.multi=[0,1,4]
for k=0,Npol-1 do plot,star(*,k),xst=1
Window,3
!P.multi=[0,1,2]
Q_o=(star(*,0)-star(*,1))/(star(*,0)+star(*,1))
U_o=(star(*,2)-star(*,3))/(star(*,2)+star(*,3))
plot,Q_o,xst=1,yrange=[0,-0.2],yst=1
;Q_o=LOWESS(x,Q_o,Nx/2,2,2)
oplot,Q_o,thick=3,color=3e5
plot,U_o,xst=1,yrange=[0,0.2],yst=1
;U_o=LOWESS(x,U_o,Nx/2,2,2)
oplot,U_o,thick=3,color=3e5

obj=obj(R,*,*,*,*)
Q_o=Q_o(R) & U_o=U_o(R)
sxaddpar,h,'CRVAL1',wave(0)
writefits,dir+'spectra_2D.fts',obj,h
Nx=N



window,0,xsize=Nx,ysize=Ny
p=0
& t=0
tv,255- bytscl(obj(*,*,p,0,t),0,1000)
;window,1
xc=210 & wx=80

norm=total(obj(xc-wx:xc+wx,*,*,*,t),1) & norm=total(norm,1) & norm=total(norm,1)
norm=norm/total(norm)*Nexp(t)
;plot,norm,psym=6,yrange=[0.9,1.1],yst=1
; вычистка частиц
ww=1
spectra_Ha=fltarr(Nx,Ny,Npol)
		for kp=0,Npol-1 do begin
	;for kx=ww,Nx-1-ww do begin
	for kx=0,Nx-1 do begin
for ky=0,Ny-1 do begin
;robomean,obj(kx-ww:kx+ww,ky,kp,*,t),3,0.5,avg_val
spectra_Ha(kx,ky,kp)=median(obj(kx,ky,kp,*,t))
spectra_Ha(kx,ky,kp)=avg_val
endfor
	endfor
		endfor
writefits,dir+'avg_spectra_2D.fts',spectra_Ha,h
cont:

obj=readfits(dir+'avg_spectra_2D.fts',h,/silent)
Nz=3
PA=fltarr(Nz) & for j=0,Nz-1 do PA(j)=sxpar(h,'PA'+string(j+1,format='(I1)'))
name=strarr(Nz) & for j=0,Nz-1 do name(j)=sxpar(h,'NAME'+string(j+1,format='(I1)'))
Nexp=intarr(Nz) & for j=0,Nz-1 do Nexp(j)=sxpar(h,'NUMEXP'+string(j+1,format='(I1)'))
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2') & Npol=4
;Ne=10
wave=findgen(Nx)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
z=sxpar(h,'Z') & print,z & z=z
Ha=6562.8
x=findgen(Nx)
yc=162   & wy=10
Flux=total(obj(*,yc-wy:yc+wy,*),2)
profile=fltarr(Nx,4)
window,0
!P.multi=[0,1,4]

for p=0,3 do begin

;удаление узких линий
BLR=remove_NLR(FLUX(*,p),rep=20,width=0.05)
V=Flux(*,p)-BLR
V_max=max(V)
V=V/V_max
Ndeg=3
lines=[6085,6300,6364,6548,6563,6583,6716,6731]*(1+z)
width=[   5,   5,   5,   5,   5,   5,   5,   5]
lines=[6548,6563,6583,6716,6731]*(1+z)
width=[   5,   5,   5,   5,   5]
res=MULTIGAUS ( wave,V,lines,FWHM=width,yfit=fit);,/plot)
print,res.flux
plot,wave,FLUX(*,p),xst=1
profile(*,p)=FLUX(*,p)-fit*V_max
;profile(*,p)=BLR
oplot,wave,profile(*,p),color=3e5
;profile(*,p)=LOWESS(x,profile(*,p),Nx/20,3,3)
;oplot,wave,profile(*,p)
endfor
Q=(profile(*,0)-profile(*,1))/(profile(*,0)+profile(*,1))
U=(profile(*,2)-profile(*,3))/(profile(*,3)+profile(*,3))
window,1,ysize=1000
!P.multi=[0,1,5]
!P.charsize=1.5
N=N_elements(wave) & angle=fltarr(N)
Flux=total(Flux,2)
plot,wave,Flux,xst=1

;Q=Q-Q_o
;U=U-U_o;
Q=Q-LOWESS(x,Q,Nx/2,3,3)+0.002;Q_o
U=U-LOWESS(x,U,Nx/2,3,3)-0.008;U_o;
plot,wave,Q,xst=1
plot,wave,U,xst=1
Pol=SQRT(U^2+Q^2)
;tmp_R=shift(Pol,-300)
;Pol(670:N-1)=tmp_R(670:N-1)
;tmp_L=shift(Pol,300)
;Pol(250:450)=tmp_L(250:450)

plot,Pol,xst=1


;вычисление угла
for k=0,N-1 do begin
angle(k)= angle_calculation(Q(k),U(k))

endfor
angle=PA(0)-angle+180;-20
;R=where(angle gt 360,ind) & if ind gt 0 then angle(R)=angle(R)-360
;R=where(angle gt 150,ind) & if ind gt 0 then angle(R)=angle(R)-180
;tmp_R=shift(angle,-300)
;angle(670:N-1)=tmp_R(670:N-1)
;angle(570:670)=angle(570:670)-10
;tmp_L=shift(angle,300)
;angle(250:450)=tmp_L(250:450)
plot, angle,xst=1
;формирование выходного файла объекта
RA=where(ABS(wave-Ha*(1+z)) lt 500, NA)
wave=wave(RA)
flux=flux(RA)

Pol=Pol(RA)*100
angle=angle(RA)
Q=Pol*cos(angle*!PI/90)
U=Pol*sin(angle*!PI/90)
window,2,xsize=400 ,ysize=1100
!P.multi=[0,1,5]
plot,wave,flux,xst=1 	&	oplot,[1,1]*(1+z)*Ha,[-1,1]*1e6,linestyle=2
plot,wave,Q,xst=1		&	oplot,[1,1]*(1+z)*Ha,[-1,1]*1e6,linestyle=2
plot,wave,U,xst=1		&	oplot,[1,1]*(1+z)*Ha,[-1,1]*1e6,linestyle=2
plot,wave,Pol,xst=1		& oplot,[1,1]*(1+z)*Ha,[-1,1]*1e6,linestyle=2
;PA=90-PA
plot,wave,angle,xst=1		& oplot,[1,1]*(1+z)*Ha,[-1,1]*1e6,linestyle=2
;запись таблицы результата

;obj='NGC4151'

;openw,1,'h:\red_data.pol\Sy1\result\'+obj+'.txt'
;for k=0,Na-1 do printf,1,wave(k),flux(k),Q(k),U(k),Pol(k),angle(k),format='(I4,E12.3,3F7.2,F8.1)'
;close,1
end
