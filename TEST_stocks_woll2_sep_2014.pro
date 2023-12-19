;test WOLLASTON2
function corr_cont,spectra,NDEG,NOISE
N=N_elements(spectra)
x=findgen(N)
cont=LOWESS(x,spectra,N/2,NDEG,NOISE)
f=goodpoly(x(0:N/2),spectra(0:N/2),1,3,fit)
fit=f(0)+x*f(1)
corr_spectra=spectra-cont+fit
return,corr_spectra
end
function median_estimate,cube
a=size(cube)
vector=fltarr(a(1))
for j=0,a(1)-1 do vector(j)=median(cube(j,*))

return,vector
end

pro stocks_WOLL2,LOGFILE,XRANGE=xrange,WIDTH=width,Z=z,ERR=err,AmpP=AmpP,WAVE_C=wave_c,$
	ATM_ABS=atm_abs,ISM=ism,dPA=dpa,CORR=corr,PS=ps,CONT=cont
tmp=str_sep(LOGFILE,'\')
path=tmp(0)
for j=1,N_elements(tmp)-3 do begin
path=path+'\'+tmp(j)
endfor
path=path+'\'
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
wdir=path+wdir(N_elements(wdir)-2)+'\'
;**********************
T=1
if not(keyword_set(ISM)) then ISM=[0,0]
if not(keyword_set(dPA)) then PA_0=317.3 ELSE PA_0=317.3+dPA
set_plot,'WIN'
if not(keyword_set(AmpP)) then Am=5 ELSE Am=AmpP
if not(keyword_set(z)) then z=0

if keyword_set(err) then err=1 else err=0
;**********************
spectra=readfits(wdir+'spectra.fit',h)


Nx=sxpar(h,'NAXIS1')  & Npol=sxpar(h,'NAXIS2') & Ncube=sxpar(h,'NAXIS4')
x=findgen(Nx)

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
date=date+' JD'+string(JD,format='(I5)')
titl=['object','unpolarized star','polarized star']

lambda_0= sxpar(h,'CRVAL1')  & d_lambda= sxpar(h,'CDELT1')
lambda=findgen(Nx)*d_lambda+lambda_0
;goto,cont1
;чтение спектральной чувствительности
if file_test(wdir+'sent.fts') then   sent=readfits(wdir+'sent.fts') ELSE sent=1

if not(keyword_set(wave_c)) then wave_c=lambda(Nx/2)

if not(keyword_set(xrange)) then xrng=[lambda(0),lambda(Nx-1)] else xrng=xrange
if keyword_set(atm_abs) then begin
;исправление атмосферного поглощения
atm=readfits(wdir+'atm_abs.fit')
;atm=(atm-1)*1.05+1
a=size(atm)
if a(0) eq 2 then atm=total(atm,2)/4
a=size(atm)


for j=0,Ncube-1 do begin
for k=0,Nexp(j)-1 do begin
for i=0,Npol-1 do begin
;atm(*,i)=shift_S(atm(*,i), 0.1)
;spectra(*,i,k,j)=spectra(*,i,k,j)/atm(*,i)
if a(0) eq 2 then spectra(*,i,k,j)=spectra(*,i,k,j)/atm(*,i) ELSE $
spectra(*,i,k,j)=spectra(*,i,k,j)/atm(*)
endfor & endfor & endfor
endif
j=0

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
;выделение спектрального интервала
R=where(lambda lt xrng(1) and lambda gt xrng(0))
;R=findgen(Nx)
lambda=lambda(R)
;вычисление параметров Стокса
Nx=N_elements(lambda)

SP=fltarr(Nx,3,Nexp(0),Ncube)
		for j=0,Ncube-1 do begin
	for k=0,Nexp(j)-1 do begin

SP(*,0,k,j)=total(spectra(R,*,k,j),2)/sent(R)
SP(*,2,k,j)=(spectra(R,1,k,j)-spectra(R,0,k,j))/(spectra(R,1,k,j)+spectra(R,0,k,j))
SP(*,1,k,j)=(spectra(R,3,k,j)-spectra(R,2,k,j))/(spectra(R,3,k,j)+spectra(R,2,k,j))
SP(*,2,k,j)=(spectra(R,1,k,j)-spectra(R,0,k,j)/D(1,k))/(spectra(R,1,k,j)+spectra(R,0,k,j)/D(1,k))
SP(*,1,k,j)=(spectra(R,3,k,j)-spectra(R,2,k,j)/D(0,k))/(spectra(R,3,k,j)+spectra(R,2,k,j)/D(0,k))
	endfor
		endfor

;Запись исходных параметров Стокса
print,'file_test',FILE_TEST(wdir+'stoks.fit')
writefits,wdir+'stoks.fit',SP,h
cont1:
;робастное сложение экспозиций
if not(keyword_set(width)) then w=2  else w=width
if w lt 2 then w=2
;медианное значение потоков
Flux=fltarr(Nx,Ncube)
for j=0,Ncube-1 do begin
for kx=0,Nx-1 do FLUX(kx,j)=median(SP(kx,0,0:Nexp(j)-1,j))
endfor
;новая шкала длин волн
cube=fltarr(Nx,Nexp(0))  & cube(*,*)=SP(*,0,0:Nexp(0)-1)
tmp=ROBUST_ESTIMATE_CUBE(cube,w=w,tresh=2); & avg_F(*,0)=tmp(*,0) & rms_F(*,0)=tmp(*,1) & wave=tmp(*,2)
NN=N_elements(tmp(*,2))  & wave=fltarr(NN) & wave(*)=lambda(tmp(*,2))
avg_S=fltarr(NN,3,Ncube)  & rms_S=avg_S
avg_P=fltarr(NN,Ncube)  & rms_P=avg_P
avg_FI=fltarr(NN,Ncube)
;
Ha=6562.8*(1+z)
T=1
;вычисление робастных оценок параметров Стокса и их ошибок
M=N_elements(wave)
for j=0,Ncube-1 do begin
cube=fltarr(Nx,Nexp(j))

for k=0,2 do begin
cube(*,*)=SP(*,k,0:Nexp(j)-1,j)

tmp=ROBUST_ESTIMATE_CUBE(cube,w=w,tresh=T)
;vector=median_estimate(cube)
;tmp=ROBUST_ESTIMATE(vector,w=w,tresh=T)
avg_S(*,k,j)=tmp(*,0) & rms_S(*,k,j)=tmp(*,1)
endfor & endfor

;поправка за стандарт нулевой поляризации
Q_null=LOWESS(findgen(NN),avg_S(*,1,1),NN/4,1,1) & Q_null=LOWESS(findgen(NN),Q_null,NN/4 ,1,1)
U_null=LOWESS(findgen(NN),avg_S(*,2,1),NN/4,1,1) & U_null=LOWESS(findgen(NN),U_null,NN/4 ,1,1)

;window,2
;!P.multi=[0,1,2]
;plot,Q_null
;plot,U_null
;goto,fin
for j=0,Ncube-1 do begin
avg_S(*,1,j)=avg_S(*,1,j)-Q_null
avg_S(*,2,j)=avg_S(*,2,j)-U_null

if keyword_set(corr) then begin
avg_S(*,1,j)=corr_cont(avg_S(*,1,j),3,1)
avg_S(*,2,j)=corr_cont(avg_S(*,2,j),3,1)
endif
;исправление наклона континуума в векторах Стокса объекта
if keyword_set(cont) and j eq 0 then begin
for k=1,2 do begin
f=goodpoly(wave,avg_S(*,k,j),1,1,fit)
avg_S(*,k,0)=avg_S(*,k,j)-fit+fit(NN/2)
endfor
endif
;приведение  плоскости поляризации к небесной сфере

TETA=2*(PA(j)-PA_0)*!PI/180.
;teta=0
print,j,teta,PA(j),PA_0,sin(TETA),cos(TETA)
tmp_Q=avg_S(*,1,j)*cos(TETA)-avg_S(*,2,j)*sin(TETA)
tmp_U=avg_S(*,1,j)*sin(TETA)+avg_S(*,2,j)*cos(TETA)
avg_S(*,1,j)=tmp_Q  & avg_S(*,2,j)=tmp_U


;учет межзвездной поляризации
if j  eq 0 then begin
avg_S(*,1,j)=avg_S(*,1,j)-ISM(0)*0.01 & avg_S(*,2,j)=avg_S(*,2,j)-ISM(1)*0.01
endif
avg_FI(*,j)=calc_atan(avg_S(*,1,j),avg_S(*,2,j))/2
avg_P(*,j)=sqrt(avg_S(*,1,j)^2+avg_S(*,2,j)^2)
rms_P(*,j)=sqrt(rms_S(*,1,j)^2+rms_S(*,2,j)^2)
endfor







;goto,fin
str_cube=['object','unpolarized star','polarized star']
if keyword_set(ps) then key_plot=1 else key_plot=0
set_plot,'WIN'
if key_plot eq 1 then begin
set_plot,'PS'
device,file=wdir+'res.ps',xsize=16,ysize=28,xoffset=2,yoffset=1
!P.charsize=0.9
endif
FOR jj=0,Ncube-1 DO BEGIN
if key_plot eq 0 then j=Ncube-jj-1 else j=jj
if key_plot eq 0 then begin
WINDOW,10+jj,xsize=600,ysize=1000,xpos=jj*40,ypos=jj*40,title=str_cube(j)
!P.charsize=1
endif
!P.multi=[0,1,1]
;рисование суммарного спектра
order=FIX(ALOG10(max(Flux(*,j))))
yrng=[0,max(Flux(*,j))]/10^order
plot,lambda,Flux(*,j)/10^order,xst=1,yrange=yrng,ytickinterval=1,yminor=5,$
		position=[0.1,0.80,0.995,0.95],/norm,xcharsize=1e-5,$
		title=STRCOMPRESS(name(j)+', '+date+', Texp='+$
		  string(exptime(j),format='(I4)')+' s, PA!Dslit!N= '+string(PA(j),format='(F6.1)')+' deg'),$
	ytitle='FLux, 10!U'+string(order,format='(I1)')+'!N, ADU/px'
if j eq 0 then oplot,[1,1]*Ha,[-1,1]*1e6,linestyle=2
;рисование поляризованного суммарного спектра
order=FIX(ALOG10(max(avg_S(*,0,j)*avg_P(*,j))))
yrng=[0,max(avg_S(*,0,j)*avg_P(*,j))]/10^order

plot,wave,avg_S(*,0,j)*avg_P(*,j)/10^order, psym=10,xst=1,thick=1+err,$
	yrange=yrng,ytickinterval=1,yminor=5,$
	position=[0.1,0.65,0.995,0.80],/norm,/noerase,xcharsize=1e-5,$
	ytitle=STRCOMPRESS('Flux*P, 10!U'+string(order)+'!N, ADU')
if j eq 0 then oplot,[1,1]*Ha,[-1,1]*1e6,linestyle=2
 dW=300
R=where(ABS(wave-Wave_c) lt dW)
;рисование параметров Стокса
titl_S=['Q-Stoks, %','U-Stoks, %']
title_S=['Q('+string(Wave_c,format='(I4)')+')=','U('+string(Wave_c,format='(I4)')+')=',$
	'P('+string(wave_c,format='(I4)')+')=','!9P!3('+string(Wave_c,format='(I4)')+')=']
for k=1,2 do begin


robomean,avg_S(R,k,j)*100,3,0.5,mean,rms
yrng=[mean-1*Am,mean+1*Am]
plot,wave,avg_S(*,k,j)*100,psym=10,yrange=yrng,yst=1, thick=err+1,xst=1,$
	ytickinterval=1,yminor=5,$
	position=[0.1,0.50-0.15*(k-1),0.995,0.65-0.15*(k-1)],/norm,/noerase,xcharsize=1e-5,	ytitle=titl_S(k-1)
	ytitle=titl_S(k-1)

xyouts,(max(wave)+min(wave))/2,yrng(0)+(yrng(1)-yrng(0))*0.85,/data,align=0.5,$
	title_S(k-1)+string(mean,format='(F6.2)')+'!9 +!3'+string(rms,format='(F5.2)')
if j eq 0 then oplot,[1,1]*Ha,[-1,1]*1e6,linestyle=2
if j eq 0 then oplot,[0,1]*1e6,[0,0],linestyle=2
if err eq 1 then oploterr,wave,avg_S(*,k,j)*100,rms_S(*,k,j)*100,psym=10
endfor
;рисование степени поляризации
robomean,avg_P(R,j)*100,2,0.5,mean
yrng=[0,mean+2*Am]
plot,wave,avg_P(*,j)*100,psym=10 ,yst=1, thick=err+1,xst=1,yrange=yrng,ytickinterval=1,yminor=5,$
	position=[0.1,0.20,0.995,0.35],/norm,/noerase,xcharsize=1e-5,ytitle='P, %'
xyouts,(max(wave)+min(wave))/2,yrng(0)+(yrng(1)-yrng(0))*0.85,/data,align=0.5,$
	title_S(2)+string(mean,format='(F6.2)')+'!9 +!3'+string(rms,format='(F5.2)')
if err eq 1 then oploterr,wave,avg_P(*,j)*100,rms_P(*,j)*100,psym=10
if j eq 0 then oplot,[1,1]*Ha,[-1,1]*1e6,linestyle=2
;рисование угла плоскости поляризации
;if j eq 0 then avg_FI(*,j)=calc_atan(avg_S(*,1,j),avg_S(*,2,j))/2
dFI=28*Am
robomean,avg_FI(R,j),3,0.5,mean,rms
if j eq 0 then begin
;R=where(avg_FI(*,j) lt 50,ind) & if ind gt 0 then avg_FI(R,j)=avg_FI(R,j)+180
endif
if j eq 0 then begin
RR=where(wave gt 6578 and wave lt 6598); correction NGC3227
;RR=where(wave gt 7610 and wave lt 7640); correction 3C273; & avg_FI(RR,j)=avg_FI(RR,j)-10  ; коррекци IRAS
;avg_FI(RR,j)=-(avg_FI(RR,j)-70)+70; correction 3C273
;avg_FI(RR,j)=-(avg_FI(RR,j)-155)+155; correction NGC3227
avg_FI(RR,j)=-(avg_FI(RR,j)-80)+80; correction NGC4051
endif
yrng=[mean-dFI,mean+dFI]
;yrng=[0,mean+dFI]
plot,wave,avg_FI(*,j),psym=10, thick=1+err,xst=1,yrange=yrng,yst=1,$
	position=[0.1,0.05,0.995,0.20],/norm,/noerase,$
	ytitle='!9P!3, deg',xtitle='Wavelength, A'
if j eq 0 then oplot,wave(RR),avg_FI(RR,j),psym=2
xyouts,(max(wave)+min(wave))/2,yrng(0)+(yrng(1)-yrng(0))*0.85,/data,align=0.5,$
	title_S(3)+string(mean,format='(F6.1)')+'!9 +!3'+string(rms,format='(F4.1)')
if err eq 1 then oploterr,wave,avg_FI(*,j),rms_P(*,j)*28*100,psym=10
if j eq 0 then begin
oplot,[1,1]*Ha,[-1,1]*1e6,linestyle=2
if ISM(0)+ISM(0) ne 0  then ISM_string=', Q!DISM!N='+string(ISM(0),format='(F5.2)')+$
	' %, U!DISM!N='+string(ISM(1),format='(F5.2)')+' %' ELSE ISM_string=' '
xyouts,0,0,wdir+' window='+string(w,format='(I2)')+ISM_string,/norm
endif
ENDFOR

if key_plot eq 1 then device,/close ELSE set_plot,'WIN'

;goto,fin


set_plot,'WIN'
if key_plot eq 0 then Window,3,xsize=600,ysize=900
!P.multi=[0,1,2]
if key_plot eq 1 then begin
set_plot,'PS'
device,file=wdir+'Angle_VS_Vel.ps'


endif
;*****************************************
Rsc=452 & titl=' Akn 120, Rsc=452 ligth days'  ;

;Rsc=963 & titl=' 3C273, Rsc=963 ligth days'  ;
Rsc=38 & titl=' NGC4051, Rsc=38 ligth days'
;Rsc=1130 & titl='IRAS13349+2438 , Rsc=1130 ligth days'
;Rsc=142 & titl='Mkn335 , Rsc=142 ligth days'
;Rsc=60& titl=' NGC5548, Rsc=60 ligth days'
;Rsc=44 & titl=' NGC4151, Rsc=44 ligth days'
;Rsc=180 & titl=' Mkn817, Rsc=180 ligth days'
 ;
;*****************************************

Ha=6562.8*(1+z)
d_wave=100
R=where(ABS(wave-Ha) lt d_wave,Nw)
;print,wave(R)
;print,Nw
Q=avg_S(R,1,0)  & U=avg_S(R,2,0)  & P=sqrt(Q^2+U^2)
rms_FI=sqrt(rms_S(R,1,0)^2+rms_S(R,2,0)^2)*1e2*28
wave=wave(R)
dV=500
Vel=(wave-Ha)/Ha*3e5-dV
FI_cont=65

avg_FI=(avg_FI(R,0)-FI_cont)
crsz=1.2
;tmp=avg_FI(76:82)-20 & tmp=-tmp+20 & avg_FI(76:82)=tmp  ; for NGC4051



plot,Vel+dV,avg_FI,xst=1 ,charsize=crsz,yrange=[-60,60]*1.5,$
	title=titl,xtickformat='(I6)',psym=10
;oplot,Vel(76:82),avg_FI(76:82),psym=2
oploterr,Vel+dV,avg_FI,rms_FI;,psym=10
oplot,[0,0]+dV,[-1,1]*100,linestyle=2
oplot,[-1,1]*2e4,[0,0],linestyle=2




;goto,fin
sign= -1
;зависимость LOG(v/c))_VS_log(tan(FI))
RL=where(Vel lt 0)  & RR=where(Vel gt 0)
avg_FI=avg_FI*sign
T=TAN(avg_FI*!PI/180)/2
rms_T=rms_FI/180.*!PI/(cos(rms_FI/180.*!PI))^2
VR=VEL(RR)/3e5  & TR=-T(RR) & rms_TR=rms_T(RR)
R=where(TR gt 0) & TR=TR(R) & VR=VR(R)*2 & rms_TR=rms_TR(R)
VL=-VEL(RL)/3e5  & TL=T(RL) & rms_TL=rms_T(RR)
R=where(TL gt 0) & TL=TL(R)  & VL=VL(R) & rms_TL=rms_TL(R)
;отбрасывание нереальных точек
;R=where(VR GT 0.006) & VR=VR(R)*1.5 & TR=TR(R)
;R=where(VL GT 0.006) & VL=VL(R)*1.5 & TL=TL(R)

;сортировка
;KL=sort(TL) & TL=TL(KL) & VL=VL(KL) & rms_TL=rms_TL(KL)
;KR=sort(TR) & TR=TR(KR) & VR=VR(KR) & rms_TR=rms_TR(KR)

;N=min([N_elements(KL),N_elements(KR)])
;TR=TR(0:N-1) & TL=TL(0:N-1)  & rms_TL=rms_TL(0:N-1)
;VR=VR(0:N-1) & VL=VL(0:N-1)  & rms_TR=rms_TR(0:N-1)



A=-2.5 & B=-0.5

f=goodpoly([ALOG10(TL),ALOG10(TR)],[ALOG10(VL),ALOG10(VR)],1,1,Fit)
print,f

LogT=[ALOG10(TL),ALOG10(TR)] & LogV=[ALOG10(VL),ALOG10(VR)]
LogV_Kepler=A+B*LogT
robomean,LogV-LogV_Kepler,1.8  ,0.5,avg_A,rms_A

avg_A=A+avg_A
print,avg_A,rms_A
Log_M_BH=ALOG10(1.78)+10+2*avg_A+ALOG10(Rsc)
err_Log_M_BH=2*rms_A
print,Log_M_BH,err_Log_M_BH
LogV_Kepler=avg_A+B*LogT
TRESH=3
R=where(ABS(ALOG10(VR)-avg_A-B*ALOG10(TR)) lt TRESH*rms_A)
plot,ALOG10(TR(R)),ALOG10(VR(R)),psym=6,yrange=[-3,0],charsize=crsz,$
	title='Log(V/c) = a + 0.5 Log(!9P!3)',$
	xtickinterval=1,xtitle='Log(tan!9P!3)',ytitle='Log(V/c)',xrange=[-3,0],xst=1
R=where(ABS(ALOG10(VL)-avg_A-B*ALOG10(TL)) lt TRESH*rms_A)
oplot,ALOG10(TL(R)),ALOG10(VL(R)),psym=2
oplot,LogT,LogV_Kepler
xyouts,-1.5,-0.3,'a = '+string(avg_A,format='(F5.2)')+'!9+!3'+string(rms_A,format='(F4.2)'),$
	charsize=crsz*1.5,align=0.5
xyouts,-1.5,-0.7,'Log(M!DBH!N/M!D!9n!N!3)='$
	+string(Log_M_BH,format='(F5.2)')+'!9+!3'+string(err_Log_M_BH,format='(F5.2)'),$
	charsize=crsz*1.5,align=0.5
if key_plot eq 1 then device,/close


set_plot,'WIN'
Window,1,xsize=600,ysize=1100
!P.multi=[0,1,6]

flux=INTERPOL(flux(*,0),lambda,wave)
plot,wave,flux,xst=1,charsize=2
plot,wave,Q,xst=1,charsize=2
plot,wave,U,xst=1,charsize=2
plot,wave,flux*P,xst=1,charsize=2
plot,wave,P,xst=1,charsize=2
plot,wave,avg_FI+FI_cont,xst=1,charsize=2

;вывод результата
;openw,1, wdir+'tanFI_VS_Vel.txt'
;for k=0,N-1 do printf,1,TL(k),VL(k),TR(k),VR(k)
;close,1

N=N_elements(wave)
openw,1, wdir+'result.txt'
for k=0,N-1 do printf,1,wave(k),vel(k),flux(k),flux(k)*P(k),Q(k)*1e2,U(k)*1e2,P(k)*1e2,avg_FI(k)+FI_cont,$
	format='(2I7,2E12.3,3F7.2,F8.1)'
close,1
fin:

END
;ath='h:\red_data.pol\3C390.3\'
path='h:\red_data.pol\AGN\'
xrng=[4600,7200]; & z_obj=0.0561 & Wc=5500


LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path=path+'LOGS')

;stocks_WOLL2,LOGFILE,xrange=[6500,8500],WIDTH=2,AmpP=5,/atm_abs,/bobo,wave_c=7200	;for 3C273.3
xrng=[7200,8000]

;stocks_WOLL2,LOGFILE,xrange=xrng,WIDTH=2,AmpP=5,wave_c=7200,z=0.1583,/atm,dPA=45,ISM=[0.51,0.1],/ps;,/err;for 3C273;,/corr;
;xrng=[6300,6800]
;stocks_WOLL2,LOGFILE,xrange=xrng,WIDTH=2,AmpP=2,wave_c=6000,z=0.0315,/atm   ,dPA= 0,/corr,ISM=[-0.5,-1],/ps;,/err;for Mkn817
;stocks_WOLL2,LOGFILE,xrange=xrng,WIDTH=2,AmpP=2,wave_c=7200,z=0.0258   ,dPA= 0,/corr,ISM=[0 ,-1 ],/ps;,/err;for Mkn335
;stocks_WOLL2,LOGFILE,xrange=xrng,WIDTH=4,AmpP=3,wave_c=6200,z=0.0172,/atm   ,dPA=00,ISM=[-0.2 ,0.0 ],/ps;,/err;for NGC5548
;stocks_WOLL2,LOGFILE,xrange=xrng,WIDTH=4,AmpP=3,wave_c=6200,z=0.0033   ,dPA=00,ISM=[ 0.0 ,0.0 ],/ps;,/err;for NGC4151
;stocks_WOLL2,LOGFILE,xrange=xrng,WIDTH=2,AmpP=5,wave_c=6200,z=0.0039,/atm   ,dPA=00,ISM=[-0.1 , 0.1 ],/ps;,/err;for NGC3227
;stocks_WOLL2,LOGFILE,xrange=xrng,WIDTH=2,AmpP=2,wave_c=5200,z=0.0323   ,dPA=0,/corr,ISM=[-0.8,0.2 ],/ps;,/err;for AKN120
stocks_WOLL2,LOGFILE,xrange=xrng,WIDTH=2,AmpP=4,wave_c=6500;,z=0.0024,/atm   ,dPA=0,ISM=[0 ,-1.4 ],/ps;,/err;for NGC4051
;xrng=[7000,7500]
;stocks_WOLL2,LOGFILE,xrange=xrng,WIDTH=1,AmpP=3,wave_c=7200,/atm,z=0.1076,/cont,dPA=0,ISM=[-1,-3],/ps;,/err;for IRAS 13349+2438
end