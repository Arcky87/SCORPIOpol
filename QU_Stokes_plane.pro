;QU_Stokes_plane
wdir='h:\red_data.pol\AGN\PG0844+349_141121\'
stoks=readfits(wdir+'stoks.fit',h)
Nx=sxpar(h,'NAXIS1') & Npol=4 & Ntarget=3
Nexp=findgen(Ntarget)
wave=findgen(Nx)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
x=findgen(Nx)
wave_c=6200  & d_wave=500 & R=where(ABS(wave-wave_c) LT d_wave)
for k=0,Ntarget-1 do Nexp(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))
print,Nexp
spectra=fltarr(Nx,Nexp(0))  & Q_stoks=spectra  & U_stoks=Q_stoks
avg_spectra=fltarr(Nx)
spectra(*,*)=stoks(*,0,0:Nexp(0)-1,0)
Q_stoks(*,*)=stoks(*,1,0:Nexp(0)-1,0)
U_stoks(*,*)=stoks(*,2,0:Nexp(0)-1,0)
window,0,xsize=600,ysize=800
!P.multi=[0,1,2]
plot,[min(wave),max(wave)],[100,8E4],xst=1,yst=1,/nodata,/ylog
for k=0,Nexp(0)-1 do begin
mean_spectra=LOWESS(x,spectra(*,k),50,3,3)
norm=total(mean_spectra(R))/7e6
spectra(*,k)=spectra(*,k)/norm
oplot,wave,spectra(*,k)
endfor

rms_spectra=fltarr(Nx)
;удаление космических частиц
for k=0,Nx-1 do begin
avg_spectra(k)=median(spectra(k,*))

robomean,spectra(k,*)-avg_spectra,3,0.5,mean,rms
rms_spectra(k)=rms
RR=where(abs(spectra(k,*)-avg_spectra-mean) gt rms,ind) & if ind gt 1 then spectra(k,RR)=avg_spectra(k)
print,ind
endfor
plot,[min(wave),max(wave)],[100,8E4],xst=1,yst=1,/nodata,/ylog
for k=0,Nexp(0)-1 do oplot,wave,spectra(*,k)
;исправление атмосферного поглощения
atm_abs=readfits(wdir+'atm_abs.fit')
atm_abs=total(atm_abs,2)/4
co=0.94
atm_abs=((atm_abs-1)*co+1)
atm_abs=shift_s(atm_abs,1.5)
wave_c=7200  & d_wave=500 & R=where(ABS(wave-wave_c) LT d_wave)
window,2,xsize=1200,ysize=500
!P.multi=[0,1,1]
R=findgen(Nx)
plot,wave(R),avg_spectra(R),xst=1,yrange=[-5E3,7e4],yst=1
avg_spectra=avg_spectra/atm_abs
oplot,wave(R),avg_spectra(R),color=3e5
oplot,wave(R),rms_spectra(R)
window,3
plot,wave(R),Q_stoks(R,0),xst=1,yrange=[-1,1]*0.05,yst=1
END