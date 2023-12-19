;***************************

dir='D:\SCORPIO\red_data.pol\NGC4151_160307\'
dir='D:\Sy1\NGC4151_160307\'
goto,cont
obj=readfits(dir+'obj_lin.fts',h)

a=size(obj)
Nx=a(1) & Ny=a(2) & Npol=a(3) & Nexp=a(4) & Nobj=a(5)

;correction flat-field
flat=readfits(dir+'avg_flat.fts')
for e=0,Nexp-1 do begin
	for p=0,Npol-1 do begin
		for o=0,Nobj-1 do begin
			obj(*,*,p,e,o)=obj(*,*,p,e,o)/flat(*,*,p)
		endfor
	endfor
endfor

;sky substraction
y=findgen(Ny)

d=40
ima=fltarr(Nx,Ny)
for o=0,Nobj-1 do begin
	for p=0,Npol-1 do begin
		for e=0,Nexp-1 do begin
			ima=obj(*,*,p,e,o)
			if ima(0,0) ne 0 then begin
				for kx=0,Nx-1 do begin
					V=ima(kx,*)
					;plot,y,V,xst=1,yrange=[0,max(V/10)],yst=1
					yy=[y(0:d),y(Ny-d/2:Ny-1)]
					VV=[V(0:d),V(Ny-d/2:Ny-1)]
					;if VV(0) ne 0 then robomean,VV,1,0.5,sky
					f=goodpoly(yy,VV,1,1)
					sky=f(0)+f(1)*y
					ima(kx,*)=V-sky
				endfor
			endif
		obj(*,*,p,e,o)=ima(*,*)
 		endfor
    endfor
 endfor

sxaddhist,'flat-field corrected',h
sxaddhist,'sky substracted',h
writefits,dir+'obj-sky.fts',obj,h


obj=readfits(dir+'obj-sky.fts',h)
a=size(obj) & Nx=a(1) & Ny=a(2) & Npol=a(3) & Nexp=a(4) & Nobj=a(5)
wave_0=sxpar(h,'CRVAL1') & d_wave=sxpar(h,'CDELT1')
wave=findgen(Nx)*d_wave+wave_0

YC=155 & WY=75
;create spectra
spectra_nuc=fltarr(Nx,Npol,Nexp,Nobj)
spectra_gal=spectra_nuc
for e=0,Nexp-1 do begin
	for p=0,Npol-1 do begin
		for o=0,Nobj-1 do begin
			spectra_nuc(*,p,e,o)=TOTAL(obj(*,yc-wy:yc+wy,p,e,o),2)
			spectra_gal(*,p,e,o)=TOTAL(obj(*,0:yc-wy,p,e,o),2)+TOTAL( obj(*,yc+wy:Ny-1,p,e,o),2)
		endfor
	endfor
endfor

cgdisplay, wid=0, xsize=1000, ysize=300
avg_gal=fltarr(Nx,4)
for kx=11,Nx-11 do begin
	for p=0,3 do begin
		robomean,(spectra_gal(kx-10:kx+10,p,*,0)),3,0.5,av,rms
		avg_gal(kx,p)=av
	endfor
endfor
cgplot,(avg_gal(*,0)-avg_gal(*,1))/(avg_gal(*,0)+avg_gal(*,1)), color='gray', yrange=[-1,1]
cgoplot,spectra_nuc(*,0,0,0)/1e5



;определение деполяризации по спектру объекта
time=fltarr(7,4)
openr,1,dir+'obj_start.txt'
readf,1,time
close,1
TR=1
x0=700 & dx=600
avg_R=fltarr(Npol,Nexp) & rms_r=avg_R
for k=0,Nexp-1 do begin
for j=0,Npol-1 do begin
;R=spectra_gal(x0:x0+dx,0,j,k)/spectra_gal(x0:x0+dx,1,j,k)
R=spectra_nuc(x0:x0+dx,0,j,k)/spectra_nuc(x0:x0+dx,1,j,k)
robomean,R,TR,0.5,mean,rms
avg_R(j,k)=mean & rms_R(j,k)=rms
endfor & endfor
;window,1,ysize=800
;!P.multi=[0,1,4]
;print,time(*,0)
;for j=0,Npol-1 do begin
;plot,time(1:5,j),avg_R(j,*) , psym=6,xrange=[5.5,8],xst=1,yrange=[1.08,1.11],yst=1
;oploterr,time(1:5,j),avg_R(j,*),rms_R(j,*),psym=6
;robomean,avg_R(j,*),3,0.5,mean,rms & print,mean,rms
;ENDFor


Q_0=0.8
U_0=-0.1
F_nuc=fltarr(Nx,Nexp) & F_gal=F_nuc
Q_nuc=F_nuc & Q_gal=Q_nuc
U_nuc=F_nuc & U_gal=U_nuc
S=fltarr(Nx,Nray,Npol)


tmp=readfits(dir+'spectra.fit')
spectra_nuc=tmp(*,*,*,0)
Ns=6 & Nx=2231 & Nexp=10
S=fltarr(Nx,4)
stoks_nuc=fltarr(Nx,Ns,Nexp)



cont:

stoks=readfits(dir+'stoks.fit')
a=size(stoks) & Nx=a(1) & Nexp=a(3)
stoks_nuc=fltarr(Nx,6,Nexp)
for k=0,Nexp-1 do begin
	stoks_nuc(*,0,k)=stoks(*,0,k) ;I
	stoks_nuc(*,1,k)=stoks(*,1,k) ;Q
	stoks_nuc(*,2,k)=stoks(*,2,k) ;U
	stoks_nuc(*,3,k)=SQRT(stoks_nuc(*,1,k)^2+stoks_nuc(*,2,k)^2)
	stoks_nuc(*,4,k)=stoks_nuc(*,0,k)*stoks_nuc(*,3,k)/100
	stoks_nuc(*,5,k)=calc_atan(stoks_nuc(*,1,k),stoks_nuc(*,2,k))/2
endfor


;window,0, ysize=800
set_plot,'PS'
device,file=dir+'polarization.ps',xsize=18,ysize=26,xoffset=2,yoffset=1,/portrait
!P.multi=[0,1,1]
!P.charsize=1
z=0.003319
Ha=6563*(1+z)
Hb=4862*(1+z)
cube=fltarr(Nx,Nexp)
Amp=[1e6,1,1,1,1e4,1]
sw=[0,1,1,0,0,1]
ss=[0,1,1,1,0,1]
ww=[2,5,5,5,2,5]
xcsz=fltarr(Ns)+1e-5 & xcsz(Ns-1)=1
ER=[0,1,1,1,1,1]
yrng=fltarr(2,Ns)
yrng(*,0)=[0,2]*0.9
yrng(*,1)=[-1,1]*0.9
yrng(*,2)=[-1,1]*0.9
yrng(*,3)=[0,2]*0.9
yrng(*,4)=[0,2]*0.9
yrng(*,5)=[-1,1]*15*0.9
xrng=[5000,7200]
ys=0.999 & dy=0.16


for k=0,Ns-1 do begin
cube(*,*)=stoks_nuc(*,k,*)
res=robust_estimate_cube(cube,W=ww(k),TRESH=2)
a=size(res); & print,a
robomean,res(a(1)/2-500/ww(k):a(1)/2+500/ww(k),0)/amp(k),3,0.5, mean
yrng(*,k)=yrng(*,k)+mean*sw(k)

plot,wave(res(*,2)),res(*,0)/amp(k) ,xrange=xrng,xst=1,psym=10,xcharsize=xcsz(k),$
	yrange=yrng(*,k),yst=1,position=[0.1,ys-(k+1)*dy,0.995,ys-k*dy],/norm,noerase=er(k),$
	xtitle='Wavelength, '+string(197B)
if SS(k) ne 0  then oploterr,wave(res(*,2)),res(*,0)/amp(k),res(*,1)/amp(k),psym=10,HATLENGTH=0, ERRCOLOR = 100
oplot,[1,1]*Hb,[-1,1]*1e8,linestyle=2
oplot,[1,1]*Ha,[-1,1]*1e8,linestyle=2
;mean_FI=LOWESS(res(*,2),res(*,0),Nx/ww(5)/2,2,2)
;;f=goodpoly(res(*,2),res(*,0),2,2,mean_FI)
;oplot,wave(res(*,2)),mean_FI,thick=2
endfor
avg_FI=res(*,0) & rms_FI=res(*,1)
mean_FI=LOWESS(res(*,2),res(*,0),Nx/ww(5)/2,2,2)

;oplot,wave(res(*,2)),mean_FI,thick=2
ytitl=['F!B!7k!3!N, 10!U6!N ADU','Q!B!7k!3!N, %','U!B!7k!3!N, %',$
	   'P!B!7k!3!N, %','F!B!7k!3!N!9. !3P!B!7k!3!N, 10!U4!N ADU',$
	   '!9P!3!B!7k!3!N, deg']
for k=0,5 do xyouts,0.03,ys-dy/2-dy*k,ytitl(k),/norm,align=0.5,orientation=90
device,/close
set_plot,'Win'

end
;***********************************************
;зависимость LOG(v/c))_VS_log(tan(FI))
;***********************************************
res=robust_estimate_cube(stoks_nuc(*,0,*),W=2,TRESH=2)
flux=res(*,0) &  lambda=wave(res(*,2))
;Ha=Hb
window,0
!P.multi=[0,1,2]
vel=(wave(res(*,2))-Ha)/Ha*3e5
plot,vel,flux,xst=1,xrange=[-1,1]*20000

res=robust_estimate_cube(stoks_nuc(*,5,*),W=ww(5),TRESH=1)
;avg_FI=res(*,0) & rms_FI=res(*,1)
vel=(wave(res(*,2))-Ha)/Ha*3e5
FI_0=80
;FI_0= 78
R=where(ABS(vel) lt 25000, ind)
if ind gt 1 then begin
vel=vel(R)
avg_FI=Avg_FI(R)
rms_FI=rms_FI(R)
endif
plot,vel,avg_FI-FI_0,xst=1,psym=10,xrange=[-1,1]*18000,xtickinterval=5000,xtickformat='(I6)'
oploterr,vel,avg_FI,rms_FI,psym=10
oplot,-[1,1]*2000,[-1,1]*30,linestyle=2

; линия На

;set_plot,'PS'
;device,file=dir+ 'Vel_VS_tanFI.ps'
smz=1.5

crsz=1.5
Vmin=0
sign=-1 & TRESH=3
Rsc=215
;vel=lambda
dV=1000
dN=3
RL=where(Vel lt -Vmin,ind_RL)  & RR=where(Vel gt Vmin,ind_RR)
print,ind_RL,ind_RL
avg_FI=(avg_FI-FI_0)*sign
avg_FI=shift(avg_FI,dN)
T=TAN(avg_FI*!PI/180)


;rms_T=rms_FI/180.*!PI/(cos(rms_FI/180.*!PI))^2
VR=VEL(RR)/3e5   & TR=-T(RR) & ;rms_TR=rms_T(RR)
R=where(TR gt 0.0) & TR=TR(R) & VR=VR(R)  & ;rms_TR=rms_TR(R)
VL=-VEL(RL)/3e5   & TL=T(RL) & ;rms_TL=rms_T(RR)
R=where(TL gt 0,ind)
if ind gt 0 then begin
TL=TL(R)  & VL=VL(R) & ;rms_TL=rms_TL(R)
endif
A=-2.5 & B=-0.5

f=goodpoly([ALOG10(TL),ALOG10(TR)],[ALOG10(VL),ALOG10(VR)],1,2,Fit)
print,'f: ',f(0),f(1)
d_logV=0
LogT=[ALOG10(TL),ALOG10(TR)] & LogV=[ALOG10(VL)+d_logV,ALOG10(VR)-d_logV]
LogV_Kepler=A+B*LogT

robomean,LogV-LogV_Kepler,tresh,0.5,avg_A,rms_A
rms_A=rms_A/2
avg_A=A+avg_A
;avg_A=-2.41  & rms_A=0.11
print,avg_A,rms_A
VDELT=Vel(1)-Vel(0)
print,Vdelt

;avg_A=avg_A-0.3
if Rsc ne 0 then begin
Log_M_BH=ALOG10(1.78)+10+2*avg_A+ALOG10(Rsc)
;Log_M_BH=ALOG10(0.58)+10+2*avg_A +ALOG10(Rsc)
err_Log_M_BH=2*rms_A
print,Log_M_BH,err_Log_M_BH
endif
kep=1
LogV_Kepler=avg_A+B*LogT
Window,1,ysize=900
!P.multi=[0,1,2]
titl='Broad H!7a!3 Fairall 9'
plot,Vel-VDELT*dN,-avg_FI,xst=1,xrange=[-1,1]*14999 ,charsize=crsz,yrange=[-1,1]*29,yst=1,$
	title=titl,xtickformat='(I6)',psym=10,xtitle='Velocity, km/s',ytitle='!7D!9P!3'

;oplot,[-1,1]*2e4,[1,1]*45,linestyle=1
;oplot,[-1,1]*2e4,-[1,1]*45,linestyle=1
;xyouts,0,80,'broad H!7a!3',align=0.5,charsize=1
oplot,[0,0] ,[-1,1]*100,linestyle=2
oplot,[-1,1]*2e4,[0,0],linestyle=2
print,'sign ',sign
;if key_plot eq  0 then col=4e4 else col=150
Vmin=2000
RPOS=WHERE(Vel ge Vmin*1.5)  & RNEG=WHERE(Vel le -Vmin)
M_BH=Log_M_BH
col=180
FI_Kepler=(ATAN(10.0^(-2*avg_A)*Vel(RPOS)^2/9E10/2)*180/!PI-90)/2
oplot,Vel(RPOS)-VDELT*dN ,-FI_Kepler*Kep  ,color=col,thick=20
FI_Kepler=(90-ATAN(10.0^(-2*avg_A)*Vel(RNEG)^2/9E10)*180/!PI)/2

oplot,Vel(RNEG)-VDELT*dN ,-FI_Kepler*Kep   ,color=col,thick=20
oplot,Vel-VDELT*dN,-avg_FI,psym=10 & oplot,[-1,1]*2e4,[0,0],linestyle=2
oploterr,vel-VDELT*dN,-avg_FI,rms_FI,xst=1,psym=10
oplot,Vel(RR)-VDELT*dN,-avg_FI(RR),psym=10,thick=2
oplot,Vel(RL)-VDELT*dN,-avg_FI(RL),psym=10,thick=2


dd_V=-0.00
d_V=0.00
smz=1
TETA=findgen(32)*!PI/16 & teta=[teta,teta(0)]
USERSYM,cos(teta) ,sin(teta)
R=where(ABS(ALOG10(VR)-avg_A-B*ALOG10(TR)) lt tresh*rms_A)  ;-0.2 was added!!!

TR_out=ALOG10(TR(R)) & VR_out=ALOG10(VR(R))-d_V+dd_V & print,N_elements(TR_out)

;openw,1,dir+'TR.txt'

;for i=0,21 do printf,1,TR_out(i),VR_out(i)
;close,1
plot,TR_out,VR_out,psym=8,yrange=[-3,0],charsize=crsz,$
	title='Log(V/c) = a + 0.5 Log(tan!9P!3)',symsize=smz,$
	xtickinterval=1,xtitle='Log(tan!9P!3)',ytitle='Log(V/c)',xrange=[-3,0.5],xst=1
oplot,[-2.75],[-2.47],psym=8,symsize=smz
xyouts,[-2.85],[-2.52],'( ) - Velocity > 0',charsize=crsz
R=where(ABS(ALOG10(VL)-avg_A-B*ALOG10(TL)) lt tresh*rms_A)
TL_out=ALOG10(TL(R)) & VL_out=ALOG10(VL(R))+d_V+dd_V & print,N_elements(TL_out)
;openw,1,dir+'TL.txt'
;for i=0,38 do printf,1,TL_out(i),VL_out(i)
;close,1

USERSYM,cos(teta),sin(teta),/fill
oplot,TL_out,VL_out,psym=8,symsize=smz
oplot,[-2.75],[-2.67],psym=8,symsize=smz
xyouts,[-2.85],[-2.72],'( ) - Velocity < 0',charsize=crsz
index=sort(LogT)
oplot,LogT(index),LogV_Kepler(index)+dd_V,linestyle=2
;avg_A=-2.14 & rms_A=0.10  & Log_M_BH=8.35 & err_Log_M_BH=0.26

xyouts,-1.3,-0.3,'a = '+string(avg_A,format='(F5.2)')+'!9+!3'+string(rms_A,format='(F4.2)'),$
	charsize=crsz ,align=0.5

if Rsc ne 0 then xyouts,-1.3,-0.7,'Log(M!DBH!N/M!D!9n!N!3)='$
	+string(Log_M_BH,format='(F5.2)')+'!9+!3'+string(err_Log_M_BH,format='(F5.2)'),$
	charsize=crsz ,align=0.5
;arrow,0,-1.5,0,-2,color=0,/data
;xyouts,0,-1.3,'!9P!3!Bmax!N',/data,align=0.5,charsize=crsz
;device,/close
;set_plot,'PS'
END;
;******************************************************************************