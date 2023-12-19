;определение уровня континнума в единицах [OIII]
;goto,mark1
!P.charsize=1
;a='8.7±1		3846±105'
;a_tmp=str_sep(STRCOMPRESS(a),' ')
;R_BLR=FLOAT(str_sep(a_tmp(0),'±'))
;  Vel=FLOAT(str_sep(a_tmp(1),'±'))
;print,R_BLR,Vel
;
;x=[14.10,1.2,2253,85]
;Mbh=ALOG10(R_BLR(0)*Vel(0)^2*0.1952)
;err=(ALOG10((R_BLR(0)+R_BLR(1))*(Vel(0)+Vel(1))^2*0.1952)-ALOG10((R_BLR(0)-R_BLR(1))*(Vel(0)-Vel(1))^2*0.1952))/2
;;print,Mbh,'±',err,format='(
;print,string(Mbh,format='(F4.2)')+'±'+string(err,format='(F4.2)')
;END
mul=2.35482
path='/home/elias/SCORPIO/sppol_pipeline_v2023.8/'
LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path=path+'LOGS')
bobo=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'/')
print,bobo(2)
bobo=bobo(2)

wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'/')

wdir=path+wdir(N_elements(wdir)-2)+'/s'


z=sxpar(read_table(LOGFILE),'redshift')
;z=0.089338
print,'z',z
mul=2.35482
;spectra_in=readfits(wdir+'spectra.fit',h)/readfits(wdir+'sent.fts')
sent=readfits(wdir+'sent.fts')
print,size(sent)
;spectra_in=readfits(wdir+'spectra.fit',h)
tmp=readfits(wdir+'avg_spectra.fit',h)
Nx=sxpar(h,'Naxis1')
spectra_in=fltarr(Nx)
spectra_in(*)=total(tmp(*,*,0),2)
spectra_in=spectra_in/sent
; Nexp=sxpar(h,'NAXIS3')
;spectra_tmp=fltarr(Nx)
;for k=0,Nx-1 do spectra_tmp(k)=median(spectra_in(k,0,*,0))
;
;spectra_in=spectra_tmp

wave_in=findgen(sxpar(h,'NAXIS1'))*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')

; оценка уровня континуума на 6000 A
;R=where(wave_in gt 5400 and wave_in lt 5500,ind)
;robomean,spectra_in(R),3,0.5,avg_cont,rms_cont

d_wave=sxpar(h,'CDELT1')
Hb=4861.0 & d_W=700 &SHI=0
R=WHERE(wave_in lt (Hb+d_W)*(1+z) and wave_in gt (Hb-d_W+SHI)*(1+z))
wave=wave_in(R) & spectra=spectra_in(R)
;определение поправки шкалы длин волн по [0III]5007
Roxy=where(ABS(wave-5007*(1+z)) lt 20)
tmp=max(spectra(Roxy),Nmax) & max_wave=wave(Roxy(Nmax))
Roxy=where(ABS(wave-max_wave) lt 10)
f=goodpoly(wave(Roxy),spectra(Roxy),2,3,Yfit)
dw=-f(1)/f(2)/2-5007*(1+z)
wave=wave-dw
N=N_elements(wave)
D=100
R1=where(wave lt (Hb-d_W+D+shi)*(1+z))
;print,wave
R2=where(wave gt (Hb+d_W-D)*(1+z))
;проведение континуума в области Hbeta
f=goodpoly([wave(R1),wave(R2)],[spectra(R1),spectra(R2)],1,3)
cont1=(f(0)+f(1)*wave)
cont_wave=[4200,4250 ,4450,4770,5090,5450,5530];,5500];,4450,4770 ,
cont_wave=cont_wave*(1+z);,4995
;cont_wave=[4580,4600 ,5150]*(1+z)
Npos=N_elements(cont_wave) & cont_pos=fltarr(Npos)
 for j=0,Npos-1 do begin
R=where(ABS(wave-cont_wave(j)) lt 10)
robomean,spectra(R),3,0.5,avg_cont
cont_pos(j)=avg_cont
cont_pos(j)=min(spectra(R))
endfor
Ndeg=3
f=goodpoly(cont_wave,cont_pos,Ndeg,2)
cont2=0
for C=0,Ndeg do cont2=cont2+f(C)*wave^C
continuum=cont2
;continuum=1
;spectra=spectra-(f(0)+f(1)*wave)

window,0,xsize=650,ysize=900
!P.multi=[0,1,2]
!P.charsize=1
x=findgen(N)
plot,wave,spectra,xst=1,yst=1,yrange=[-0.5,max(spectra)],title=sxpar(h,'NAME1')
oplot,[1,1]*5007*(1+z),[0,1e6],linestyle=1

;oplot,wave(Roxy),Yfit,thick=2
oplot,cont_wave,cont_pos,psym=6
oplot,wave,continuum,linestyle=2
oplot,wave,spectra-continuum,psym=6,symsize=0.1
oplot,[wave(0),wave(N-1)],[0,0],linestyle=2
;oplot,wave,LOWESS(x,cont,N/4,2)
;гаусс-анализ компонент
; line=[4959,4924,4924,4861,4950]
;doubl=[ 48, 48,    94,   0, 0]
;    w=[ 20, 20,   40,  20];, 20]
;Akn120, Mkn1501
 line=[4959,4861,4861];,5018]
doubl=[ 48 ,   0,  0 ];, 0]
    w=[  10,  10, 40];, 40]
;PG0844+349
;line=[4959,4861,4861];,4950]
;doubl=[ 48,   0,   0];, 0]
;    w=[ 10,  10, 40 ];, 20]


line=line*(1+z) & doubl=doubl*(1+z)
c=3.0e5
Nline=N_elements(line)

;narrow component
spectra=spectra-continuum
res_n=MULTIGAUS ( wave, spectra , line, FWHM=w,SIGMA=SIGMA_N,DOUBL_LEN=doubl)
;                          [,/SILENT] [ABSORP=ABSORP] [,/LIM_FWHM=], [,/DOUBL_LEN]
;                          [,/DOUBL_RATIO=] [,/FIXRATIO=],YFIT=YFIT,sigma=sigma)
param=fltarr(4)
par=fltarr(4,Nline)
amp=[1,1/(1+z),1,1]
component=fltarr(N,7)
dz=0
for j=0,Nline-1 do begin
param(0)=res_n(j).max
param(1)=res_n(j).center
param(2)=res_n(j).FWHM/mul
param(3)=res_n(j).ratio

if j eq 4 then begin
;param(2)=param(2)/5
;param(1)=param(1)+25
;param(0)=param(0)/2
endif
;print,param/(1+z)
par(*,j)=param*amp
print,par(*,j)

;определение поправки по длине волны по [OIII]4959

component(*,j)=gauss_profile(wave,param)

param(1)=param(1)+doubl(j)
component(*,j)=component(*,j)+gauss_profile(wave,param)*param(3)
oplot,wave,component(*,j),color=3e5
endfor


ymax=max(spectra+continuum)*0.9
for j=0,Nline-1 do xyouts,wave(0) , ymax-(ymax/5/4)*j,string(par(*,j),format='(4F10.2 )')
print,dz
;dz=(dz(2) -line(0)/(1+z))/line(0)*(1+z)
;dz=(total(dz)/2-line(0)/(1+z))/line(0)*(1+z)
;print,'dz=',dz & dz=0
narrow=fltarr(N,2)
narrow=component(*,0:1)
broad=spectra-total(narrow,2);-component(*,4)
broad=spectra-component(*,0)-component(*,1); -component(*,3)
broad=shift(broad,0)

smooth_broad=LOWESS(findgen(N),broad,5,2,2)
;smooth_broad=broad
; вычисление полуширины Hbeta
wg=100
R=where(broad gt 0.5*max(broad),FWHM)
vel=(wave-Hb*(1+z))/Hb/(1+z)*3.0e5
dz=total(vel(R))/FWHM
d_wave=ABS(wave(1)-wave(0))
FWHM=FLOAT(FWHM)/mul*d_wave/4861.0/(1+z)*3e5
Ng=N_elements(smooth_broad)
;вписывание гауссианы в Hbeta
Hbeta=fltarr(Ng)

Hbeta(Ng/2-wg:Ng/2+wg)=gaussfit(wave(Ng/2-wg:Ng/2+wg),smooth_broad(Ng/2-wg:Ng/2+wg),G,sigma=S,Nterms=3)
print, G(1),G(2)
G(2)=G(2)/4861/(1+z)*3E5*mul
S(2)=S(2)/4861/(1+z)*3E5*mul
;определение полуширины гауссианы

R=where(smooth_broad(Ng/2-wg:Ng/2+wg) gt 0.5*max(smooth_broad(Ng/2-wg:Ng/2+wg)),FWHM)
FWHM=FLOAT(FWHM)*d_wave/4861.0/(1+z)*3e5

;g(1)=dz
;print,G(1)/c ,s(1)/c
;print,G(2) ,S(2)
wave=wave/(1+z)
SHI=90 & AMPL=0.0
broad=broad-SHIFT(Hbeta,-shi)*AMPL
plot,wave,broad,xst=1,yrange=[-0.2,1.1]*max(broad),yst=1,$
	title='FWHM='+string(FWHM,format='(I4)') +' km/s,    !7r!3V(Gauss)*2.355='+$
                        string(G(2),format='(I4)')+'!9+!3'+ string(2*S(2),format='(I3)')+' km/s';+$
      ;  ', Zg='+string(G(1)/c ,format='(F8.5)')+'!9+!3'+string(S(1)/c,format='(F7.5)')
oplot,wave,smooth_broad
oplot,wave,Hbeta,color=3e5,thick=3
;oplot,wave,SHIFT(Hbeta,-shi)*ampl,color=3e5,thick=3
;вычисление потока FeII
R_blu=where(wave gt 4450 and wave lt 4731)
R_red=where(wave gt 5100 and wave lt 5450)
R_blu=total(broad(R_blu))/total(Hbeta) ; & R_blu=0.266
R_red=total(broad(R_red))/total(Hbeta)
xyouts,4600,-0.12*max(broad),'FeII/H!7b!3='+string(R_blu,format='(F5.3)'),align=0.5,charsize=1.5
xyouts,5247,-0.12*max(broad),'FeII/H!7b!3='+string(R_red,format='(F5.3)'),align=0.5,charsize=1.5
;вывод результатов
openw,1,'h:\red_data.pol\Sy1\result\'+bobo+'.txt'
for k=0,Ng-1 do printf,1,wave(k)*(1+z),spectra(k)
close,1
END
















;for j=0,3 do broad=broad-component(*,j)
;broad=broad
;plot,wave,broad,xst=1
;;проведение континнума прямой линией
;R1=where(wave lt 4900) & R2=where(wave gt 5250)
;f=goodpoly([wave(R1),wave(R2)],[broad(R1),broad(R2)],1,3)
;cont=f(0)+f(1)*wave
;oplot,wave,cont
;line=[5080,5100, 5250]
;w=[100,150,100]
;res_b=MULTIGAUS ( wave, broad-cont, line, FWHM=w,SIGMA=SIGMA_b)
;
;component=fltarr(N,3)
;for j=0,2 do begin
;param(0)=res_b(j).max
;param(1)=res_b(j).center
;param(2)=res_b(j).FWHM/mul
;component(*,j)=gauss_profile(wave,param)
;oplot,wave,component(*,j)+cont
;endfor
;broad=fltarr(N,2)
;broad=component
;
;END
;window,0,xsize=600,ysize=550
;!P.multi=[0,1,1]
;plot,wave,spectra,xst=1,yrange=[-0.1,1.1]*max(spectra),yst=1,psym=10
;oplot,wave,cont
;sub=spectra-cont
;for j=0,2 do begin
;oplot,wave,narrow(*,j)+cont
;sub=sub-narrow(*,j)
;endfor
;for j=0,2 do begin
;oplot,wave,broad(*,j)+cont
;sub=sub-broad(*,j)
;endfor
;END
;
;
;oplot,wave,sub
;oplot,[wave(0),wave(N-1)],[0,0],linestyle=2
;oxy=res_n(2).flux+res_n(3).flux
;avg_broad=(total(res_b.flux)+res_n(4).flux)/oxy
;rms_broad=SQRT(total(sigma_b.flux^2)+sigma_n(4).flux^2)/oxy/4/3
;
;print,' BROAD: flux',avg_broad,rms_broad
;;oxy=total(narrow(*,3)+narrow(*,2))*d_wave
;oxy=oxy/100.*d_wave
;kk=1.7002
;print,LOGFILE,kk*avg_cont/oxy,kk*rms_cont/oxy
;end
;output=[1.96,2.62,2.03,1.56,2.21,1.72,1.65,1.33]
; input=[0.58,1.54,1.19,0.92,0.76,1.01,0.97,0.78]
; window,2
; plot,input,output,psym=6
; f=goodpoly(input,output,1,1,Yfit)
; oplot,input,Yfit
; print,f
;
;end
