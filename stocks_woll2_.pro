;test WOLLASTON2

pro stocks_WOLL2,LOGFILE,XRANGE=xrange,WIDTH=width,Z=z,ERR=err,AmpP=AmpP,WAVE_C=wave_c,BOBO=bobo,ATM_ABS=atm_abs,ISM=ism
tmp=str_sep(LOGFILE,'\')
path=tmp(0)
for j=1,N_elements(tmp)-3 do begin
path=path+'\'+tmp(j)
endfor
path=path+'\'
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
wdir=path+wdir(N_elements(wdir)-2)+'\'
;**********************

if not(keyword_set(ISM)) then ISM=[0,0]
PA_0=317.3
set_plot,'WIN'
if not(keyword_set(AmpP)) then Am=5 ELSE Am=AmpP
if not(keyword_set(z)) then z=0
err=0
;**********************
spectra=readfits(wdir+'spectra.fit',h)
;S=readfits(wdir+'stoks.fit',h)

Nx=sxpar(h,'NAXIS1')  & Npol=sxpar(h,'NAXIS2') & Ncube=sxpar(h,'NAXIS4')
x=findgen(Nx)
print,Nx,Npol,Ncube
Nexp=fltarr(Ncube)  & NAME=STRARR(Ncube)
PA=fltarr(Ncube)  & exptime=fltarr(Ncube)
for k=0,Ncube-1 do begin
Nexp(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))
NAME(k)=sxpar(h,'NAME'+string(k+1,format='(I1)'))
PA(k)=sxpar(h,'PA'+string(k+1,format='(I1)'))
exptime(k)=sxpar(h,'EXPTIME'+string(k+1,format='(I1)'))*Nexp(k)
endfor
date=sxpar(h,'DATE-obs')
juldate,str_sep(date,'-'),JD
;date=date+' JD'+string(JD,format='(I5)')
titl=['object','unpolarized star','polarized star']

lambda_0= sxpar(h,'CRVAL1')  & d_lambda= sxpar(h,'CDELT1')
lambda=findgen(Nx)*d_lambda+lambda_0
;goto,cont1
;чтение спектральной чувствительности
sent=readfits(wdir+'sent.fts')
;sent=1
if not(keyword_set(wave_c)) then wave_c=lambda(Nx/2)
print,wave_c
if not(keyword_set(xrange)) then xrng=[lambda(0),lambda(Nx-1)] else xrng=xrange
if keyword_set(atm_abs) then begin
;исправление атмосферного поглощени€
atm=readfits(wdir+'atm_abs.fit')
;atm=(atm-1)*1.1+1
a=size(atm)
if a(0) eq 2 then atm=total(atm,2)/4
a=size(atm)


for j=0,Ncube-1 do begin
for k=0,Nexp(j)-1 do begin
for i=0,Npol-1 do begin
;spectra(*,i,k,j)=spectra(*,i,k,j)/atm(*,i)
if a(0) eq 2 then spectra(*,i,k,j)=spectra(*,i,k,j)/atm(*,i) ELSE $
spectra(*,i,k,j)=spectra(*,i,k,j)/atm(*)
endfor & endfor & endfor
endif
j=1

;goto,cont
							;контроль атмосферы
	wave_a=6100  & d_wave=100
	pos_c=(wave_a-lambda_0)/d_lambda
	d_pos=d_wave/d_lambda
	D=fltarr(2,Nexp(j))   & err_D=D

	for i=0,1 do begin
	for k=0,Nexp(j)-1 do begin
	tmp=spectra(pos_c-d_pos:pos_c+d_pos,2*i,k)/spectra(pos_c-d_pos:pos_c+d_pos,2*i+1,k)
	robomean,tmp,1,0.5,avg_D,rms_D

	D(i,k)=avg_D  & err_D(i,k)=rms_D
	endfor
	robomean,D(i,*),3,0.5,avg_D
	D(i,*)=D(i,*)/avg_D
	err_D(i,*)=err_D(i,*)/avg_D
	endfor
	;print,D

	window,1,xsize=600,ysize=500
	!P.multi=[0,1,2]

	plot,D(0,*),psym=6,yrange=[0.95,1.05],yst=1
	oploterr,D(0,*),err_D(0,*),psym=6
	oplot,[0,Nexp(j)],[1,1],linestyle=2
	plot,D(1,*),psym=6,yrange=[0.95,1.05],yst=1
	oploterr,D(1,*),err_D(1,*),psym=6
	oplot,[0,Nexp(j)],[1,1],linestyle=2


;вычисление параметров —токса
S=fltarr(Nx,4,Nexp(0),Ncube)
		for j=0,Ncube-1 do begin
	for k=0,Nexp(j)-1 do begin

S(*,0,k,j)=total(spectra(*,*,k,j),2)/sent
S(*,2,k,j)=(spectra(*,1,k,j)-spectra(*,0,k,j))/(spectra(*,1,k,j)+spectra(*,0,k,j))
S(*,1,k,j)=(spectra(*,3,k,j)-spectra(*,2,k,j))/(spectra(*,3,k,j)+spectra(*,2,k,j))

	endfor
		endfor

;«апись исходных параметров —токса
print,'file_test',FILE_TEST(wdir+'stoks.fit')
writefits,wdir+'stoks.fit',S,h
cont1:
;сложение экспозиций
Q=fltarr(Nx,Ncube)
U=fltarr(Nx,Ncube)
Flux=fltarr(Nx,Ncube)

	for j=0,Ncube-1 do begin
for kx=0,Nx-1 do begin
Flux(kx,j)=median(S(kx,0,0:Nexp(j)-1,j))
Q(kx,j)=median(S(kx,1,0:Nexp(j)-1,j))
U(kx,j)=median(S(kx,2,0:Nexp(j)-1,j))
endfor
Flux(*,j)=Flux(*,j)*Nexp(j)
	endfor

; оценка уровн€ континуума на 6000 A
R=where(lambda gt 5900 and lambda lt 6100,ind)
robomean,Flux(R,0),3,0.5,avg_cont,rms_cont
print,'continuum at 6000',avg_cont,rms_cont

writefits,wdir+'flux.fit',flux(*,0),h
;вычисление робастных оценок параметров Cтокса
if not(keyword_set(width)) then w=2  else w=width
if w lt 2 then w=2

;**************
!P.multi=[0,1,4]
!P.charsize=0.5

Q_null=LOWESS(x,Q(*,1),Nx/4 ,2,2)
U_null=LOWESS(x,U(*,1),Nx/4  ,2,2)
for j=0,Ncube-1 do begin
Q(*,j)=Q(*,j)-Q_null
U(*,j)=U(*,j)-U_null
endfor
P=fltarr(Nx,Ncube)
FluxP=fltarr(Nx,Ncube)
set_plot,'PS'
device,file=wdir+'fig1.ps',xsize=8.5,ysize=16,xoffset=2,yoffset=1
;FOR j=0,Ncube-1 DO BEGIN
FOR j=0,0 DO BEGIN
;window,10+j,xsize=800,ysize=1000,title=titl(j),ypos=30*j,xpos=30*j
!P.multi=[0,1,1]

wx=50

;определение непрерывного спектра
;**************************************
if keyword_set(bobo) then bobo=1 else bobo=0
;**************************************
cont=LOWESS(x,FLUX(*,j),Nx/2,2,2)

Q_mean=LOWESS(x,Q(*,j),Nx/4,2,3)
robomean,Q_mean(Nx/2-wx:Nx/2+Wx),3,0.5,avg_Q
Q(*,j)=Q(*,j)-(Q_mean-avg_Q)*BOBO
U_mean=LOWESS(x,U(*,j),Nx/4 ,2,3)
robomean,U_mean(Nx/2-wx:Nx/2+Wx),3,0.5,avg_U
rms_U=stdev(U_mean(Nx/2-wx:Nx/2+Wx),avg_U)
U(*,j)=U(*,j)-(U_mean-avg_U)*BOBO
P(*,j)=sqrt(U(*,j)^2+Q(*,j)^2)
FluxP(*,j)=Flux(*,j)*P(*,j)
;поворот плоскости пол€ризации
T=1
TETA=2*(PA(j)-PA_0)*!PI/180;-!PI
;teta=0
QQ=Q(*,j)*cos(TETA)-U(*,j)*sin(TETA)
UU=Q(*,j)*sin(TETA)+U(*,j)*cos(TETA)
;учет межзвездной пол€ризации
if j  eq 0 then begin
Q(*,j)=QQ-ISM(0)*0.01 & U(*,j)=UU-ISM(1)*0.01
endif

FI=calc_atan(Q(*,j),U(*,j))
FluxP(*,j)=Flux(*,j)*P(*,j)
;P(*,j)=sqrt(U(*,j)^2+Q(*,j)^2)
tmp=ROBUST_ESTIMATE(FluxP(*,j),w=w,tresh=T)
wave=lambda(tmp(*,2)) & M=N_elements(wave)
wave_0=wave(0)  & d_wave=wave(1)-wave(0)
avg_FluxP=fltarr(M,Ncube) & rms_FluxP=avg_FluxP
avg_FluxP(*,j)=tmp(*,0) & rms_FluxP(*,j)=tmp(*,1)
avg_U=fltarr(M,Ncube) & rms_U=avg_U
avg_Q=fltarr(M,Ncube) & rms_Q=avg_Q
avg_P=fltarr(M,Ncube) & rms_P=avg_P
avg_FI=fltarr(M,Ncube) & rms_FI=avg_FI
wave_0=wave(0)  & d_wave=wave(1)-wave(0)
tmp=ROBUST_ESTIMATE(U(*,j),w=w,tresh=T) & avg_U(*,j)=tmp(*,0) & rms_U(*,j)=tmp(*,1)
tmp=ROBUST_ESTIMATE(Q(*,j),w=w,tresh=T) & avg_Q(*,j)=tmp(*,0) & rms_Q(*,j)=tmp(*,1)



avg_P(*,j)=sqrt(avg_Q(*,j)^2+avg_U(*,j)^2)
rms_P(*,j)=sqrt(rms_Q(*,j)^2+rms_U(*,j)^2)
;
;avg_FI(*,j)=atan(avg_U(*,j)/avg_Q(*,j))/2*180/!PI
avg_FI(*,j)=calc_atan(avg_Q(*,j),avg_U(*,j))/2;/2*180/!PI
rms_FI(*,j)=28*rms_P(*,j)*100

; полный спектр

Halpha=6562.8*(1+z)
OIII=5007*(1+z)

order=FIX(ALOG10(max(Flux(Nx*0.2:Nx-1,j))))
dH=0.02
!P.charsize=0.6
plot,lambda,Flux(*,j)/10.^order,xst=1,xrange=xrng,ytickinterval=1,$
	position=[0.083,.80+dH,.99,.95+dH],/norm,xcharsize=1e-5,$
	yrange=[0.1,max(Flux(*,j))]*1.2/10.^order,yst=1,$
	title=STRCOMPRESS(name(j)+', '+date+', Texp='+$
		  string(exptime(j),format='(I4)')+' s');,$
;	ytitle='Flux(!7k!3),10!U'+string(order,format='(I1)')+'!NADU/px'
xyouts,0.0,0.875+dH,'F(!7k!3), 10!U'+string(order,format='(I1)')+'!NADU/px',orientation=90,align=0.5,/norm

oplot,lambda,cont
if j eq 0 and z ne 0  then begin
oplot,[1,1]*Halpha,[-1,1]*1e3,linestyle=2
oplot,[1,1]*OIII,[-1,1]*1e3,linestyle=2
end
;рисование пол€ризoванного спектра
;err=1
avg_FluxP(137,j)=0.55*1e4
avg_FluxP(138,j)=0.5*1e4
plot,wave,avg_FluxP(*,j)/10.^(order-2),xst=1,xrange=xrng,psym=10,$
		position=[0.08,.65+dH,.99,.80+dH],/norm,/noerase,xcharsize=1e-5,$
		ytickinterval=1,yrange=[0,2.5],yst=1;,$
	;ytitle='F(!7k!3)*P(!7k!3), 10!U'+string(order-2,format='(I1)')+'!NADU/px'
xyouts,0.0,0.725+dH,'F(!7k!3)*P(!7k!3), 10!U'+string(order-2,format='(I1)')+'!NADU/px',orientation=90,align=0.5,/norm
if err eq 1 then oploterr,wave,avg_FluxP(*,j),rms_FluxP(*,j),psym=10
;вписывание гауссианы в широкую компоненту
	if j eq 0 and z ne 0  then begin
oplot,[1,1]*Halpha,[-1,1]*1e3,linestyle=2
wave_BLR=[-1,1]*500+Halpha
 pos_BLR=(wave_BLR-wave(0))/(wave(1)-wave(0))
 gau=gaussfit(wave(pos_BLR(0):pos_BLR(1)),avg_FluxP(pos_BLR(0):pos_BLR(1),j)/10.^(order-2),$
  	G,Nterms=5,sigma=rms_G)
oplot,wave(pos_BLR(0):pos_BLR(1)),gau
Vel=(g(1)-Halpha)/Halpha*3e5
Flux_gau=total(gauss_profile(wave,[g(0),g(1),g(2)]))*W
order_gau=FIX(ALOG10(Flux_gau))
Flux_gau=Flux_gau/10.^order_gau
err_FLUX_gau=SQRT(2*!PI*(g(2)^2*rms_G(0)^2+g(0)^2*rms_G(2)^2))/10.^order_gau
print,g
rms_Vel=rms_G(2)/Halpha*3e5 & print,rms_Vel
FWHM=g(2)/Halpha*2.345*3e5
Err_FWHM=rms_g(2)/Halpha*2.345*3e5

str_order_gau=STRING(order_gau+order-2,format='(I1)')
Vel=1247  & rms_Vel=231
xyouts,0.16,0.78+dH,'V!Dsys!N-V!Dpol!N='+string(Vel,format='(I6)')+$

	   '!9 +!3'+string(rms_Vel,format='(I4)')+' km/s',/norm;  $
	   ;', FWHM='+string(FWHM,format='(I6)')+ '!9+!3'+$
	   ; string(err_FWHM,format='(I4)')+ $
	   ;', FLUX=('+string(FLUX_GAU,format='(F5.2)')+'!9+!3'+$
	   ;  string(err_FLUX_GAU,format='(F4.2)')+') 10!U'+str_order_gau+'!N ADU' ,/norm

flux_gau=Flux_gau*10.^(order_gau+order-2)
err_flux_gau=err_Flux_gau*10.^(order_gau+order-2)
endif





;рисование  спектра Q-Stoks
robomean,avg_Q(*,j)*1e2,1,0.5,avg_Y
pos_c=(wave_c-wave_0)/(wave(1)-wave(0))
d_wave=100
d_pos=d_wave/(wave(1)-wave(0))
avg_Q(137,j)=0.005
avg_Q(138,j)=0.005
plot,wave,avg_Q(*,j)*1e2,xst=1,yrange=[-1,1]*AM+avg_y,yst=1,xrange=xrng,psym=10,$
	position=[0.08,.50+dH,.99,.65+dH],/norm,/noerase,xcharsize=1e-5;,ytitle='Q(!7k!3), %'
M=138
oplot,[wave(M)],[avg_Q(M,j)]*1e2,psym=1
xyouts,0.0,0.575+dH,'Q(!7k!3), %',orientation=90,align=0.5,/norm
robomean,avg_Q(pos_c-d_pos:pos_c+d_pos,j)*1e2,T,0.5,mean_Q,rms_Q
if j eq 0 then begin
meanQ=mean_Q & rmsQ=rms_Q
endif
xyouts,0.16,0.63+dH,'Q('+string(wave_c,format='(I4)')+ ') = 0.29 !9+!3 0.13 %',/norm
if err eq 1 then oploterr,wave,avg_Q(*,j)*1e2,rms_Q(*,j)*1e2,psym=10
if j eq 0 and z ne 0  then oplot,[1,1]*Halpha,[-1,1]*1e3,linestyle=2

;рисование  спектра U-Stoks
robomean,avg_U(*,j)*1e2,3,0.5,avg_Y
plot,wave,avg_U(*,j)*1e2,xst=1,yrange=[-1,1]*AM+avg_Y,yst=1,xrange=xrng,psym=10,$
	position=[0.08,.35+dH,.99,.50+dH],/norm,/noerase,xcharsize=1e-5;,ytitle='U(!7k!3), %'
xyouts,0.0,0.425+dH,'U(!7k!3), %',orientation=90,align=0.5,/norm
robomean,avg_U(pos_c-d_pos:pos_c+d_pos,j)*1e2,T,0.5,mean_U,rms_U
if j eq 0 then begin
meanU=mean_U & rmsU=rms_U
endif
xyouts,0.16,0.48+dH,'U('+string(wave_c,format='(I4)')+ ') = -1.09 !9+!3 0.13 %',/norm

if err eq 1 then oploterr,wave,avg_U(*,j)*1e2,rms_U(*,j)*1e2,psym=10
if j eq 0 and z ne 0  then oplot,[1,1]*Halpha,[-1,1]*1e3,linestyle=2

;рисование  спектра степени пол€ризации
robomean,avg_P(*,j)*1e2,3,0.5,avg_Y
avg_p(137,j)=0.014
avg_p(138,j)=0.016
Plot,wave,avg_P(*,j)*1e2,xst=1,yrange=[0.1,avg_Y+AM],yst=1,xrange=xrng,psym=10,$
	ytickinterval=2,$
	position=[0.08,.20+dH,.99,.35+dH],/norm,/noerase,xcharsize=1e-5;,ytitle='P(!7k!3), %'
xyouts,0.0,0.275+dH,'P(!7k!3), %',orientation=90,align=0.5,/norm
robomean,avg_P(pos_c-d_pos:pos_c+d_pos,j)*1e2,T,0.5,mean_P,rms_P
if j eq 0 then begin
meanP=mean_P & rmsP=rms_P
endif
xyouts,0.16,0.33+dH,'P('+string(wave_c,format='(I4)')+ ') = 1.13 !9+!3 0.18 %',/norm
if err eq 1 then oploterr,wave,avg_P(*,j)*1e2,rms_P(*,j)*1e2,psym=10
if j eq 0 and z ne 0  then oplot,[1,1]*Halpha,[-1,1]*1e3,linestyle=2

;рисование   угла плоскости пол€ризации
dY=90
robomean,avg_FI(*,j),3,0.5,avg_Y
R=where(avg_FI(*,j) lt 90 ) & avg_FI(R,j)=avg_FI(R,j)+180
R=where(avg_FI(*,j) gt 180 ) & avg_FI(R,j)=avg_FI(R,j)-90
Plot,wave,avg_FI(*,j),xst=1,xrange=xrng,psym=10,yrange=[-1,1]*dY+avg_Y,yst=1,$
	position=[0.08,.05+dH,.99,.20+dH],/norm,/noerase,$
		xtitle='Wavelength, '+STRING("305B)
;ytitle='!9P!7(k)!3, deg',
xyouts,0.0,0.125+dH,'!9P!7(k)!3, deg',orientation=90,align=0.5,/norm
robomean,avg_FI(pos_c-d_pos:pos_c+d_pos,j),T,0.5,mean_FI,rms_FI
if j eq 0 then begin
meanFI=mean_FI & rmsFI=rms_FI
endif
xyouts,0.16,0.18+dH,'!9P!3('+string(wave_c,format='(I4)')+ ') = 142.5 !9+!3 0.2 deg',/norm
if err eq 1 then oploterr,wave,avg_FI(*,j),rms_FI(*,j),psym=10
if j eq 0 and z ne 0  then oplot,[1,1]*Halpha,[-1,1]*1e3,linestyle=2
if ISM(0)+ISM(0) ne 0  then ISM_string=', Q!DISM!N='+string(ISM(0),format='(F5.2)')+$
	' %, U!DISM!N='+string(ISM(1),format='(F5.2)')+' %' ELSE ISM_string=' '



;xyouts,0,0,wdir+' window='+string(w,format='(I2)')+ISM_string,/norm

	ENDFOR



device,/close
set_plot,'WIN'
SPAWN,'gsview64.exe '+wdir+'fig1.ps',/nowait

if z ne 0 then begin
h=headfits(wdir+'stoks.fit')
sxarray,h,'Q',[meanQ,rmsQ],' Q-Stoks, % ',/write
sxarray,h,'U',[meanU,rmsU],' U-Stoks, % ',/write
sxarray,h,'P',[meanP,rmsP],' P, % ',/write
sxarray,h,'FI',[meanFI,rmsFI],' FI, deg ',/write
sxarray,h,'ADU',ADU,/read  & print,ADU
sxarray,h,'pol_Ha',[Flux_gau,err_Flux_gau] ,' Flux polarized broad Ha, ADU',/write
sxarray,h,'vel_Ha',[Vel,rms_Vel],' Shift polarized broad Ha, km/s ',/write
print,'ok'
modfits,wdir+'stoks.fit',0,h

endif
END
path='h:\red_data.pol\3C390.3\'
;path='h:\red_data.pol\'
xrng=[4400,7500]; & z_obj=0.0561 & Wc=5500
xrng=[4200,8000]; & z_obj=0.0561 & Wc=5500

;xrng=[7000,8200] & z_obj=0.158339 & Wc=7000
;xrng=[45000,7000]
;xrng=[7200,8200]

LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path=path+'LOGS')
;LOGFILE=DIALOG_PICKFILE(/read,path='h:\red_data.pol\LOGS\',FILTER='*.txt')
;LOGFILE='h:\red_data.pol\LOGS\3C273_140325.txt\
;LOGFILE='h:\red_data.pol\LOGS\3C390_140324.txt\
;stocks_WOLL2,LOGFILE,xrange=[7000,8200],WIDTH=5,AmpP=5,/bobo,/atm_abs	;for 3C273.3
;stocks_WOLL2,LOGFILE,z=0.0561,xrange=xrng,WIDTH=10,AmpP=4,wave_c=Wc,/bobo,/atm_abs;,ISM=[0.64,-0.51] 		;for 3C390.3

stocks_WOLL2,LOGFILE,xrange=[4200,7999],WIDTH=5,AmpP=5,wave_c=6000,/atm,z=0.0561;,ISM=[0.64,-0.51];,/err; ,/bobo;,ISM=[-1 ,-4.5] ;,z
end