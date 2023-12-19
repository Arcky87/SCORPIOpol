
function create_continuum,wave,spectra,wave_cont,NDEG=Ndeg
d_wave=wave(1)-wave(0)
Nc=N_elements(wave_cont)
value_cont=fltarr(Nc)
;проведение континуума
dc=5 ;интервал усреднения
if not(keyword_set(Ndeg)) then Ndeg=2

for j=0,Nc-1 do begin
Rc=where(ABS(wave-wave_cont(j)) lt d_wave*dc,ind)
if ind gt 0 then begin
robomean,spectra(Rc),3,0.5,mean
value_cont(j)=mean
endif
endfor
Rc=where(value_cont gt 0)

Fc=goodpoly(wave_cont(Rc),value_cont(Rc),Ndeg,3,fit)
continuum=0 & for j=0,Ndeg do continuum=continuum+wave^j*Fc(j)
return,continuum
end

;oplot,wave_cont ,y_cont, psym=6

;CREATE AND CORRECTION TOTAL SPECTRA
line_neon=[5015.68,6678.15]
 line_sky=[5577.34,6300.30]
 ion=['[NeIII]','H!7f!3' ,'H!7e!3','H!7d!3','H!7c!3','[OIII]','HeII','H!7b!3','[OIII]','[0III]',$
 	 '[FeVII]', 'HeI','[FeVII]','[OI]','[OI]','[NII]','H!7a!3','[NII]','[SII]','          [SII]']
line=[  3868.7,3889.05,3970.1, 4104.74, 4340.47, 4363.21,  4686, 4861.3,  4958.9,  5006.9,$
		  5721,5875.7,     6087,6300.2,6363.9, 6548.1,  6562.8, 6583.4, 6716.4, 6730.3]
     N_line=N_elements(line)
ind_print=1
date='151210' & bin=2 & titl='NGC7469, 10 dec 2015, resolution 5 '+STRING(197B)+', seeing 0.9", total exposure 14400 s'
;date='151013'
;date='141123'
date='101102' & bin=1 & titl='NGC7469, 2 nov 2010, resolution 8 '+STRING(197B)+', seeing 1.5", total exposure 8000 s'
;z=0.016317
;dir='h:\red_data.pol\AGN\NGC7469_'+date+'\'
date='131105'
;date='140324'

dir='h:\red_data.pol\AGN\E1841+643_'+date+'\'

;goto,cont
;чтение FITS-заголовка
h=headfits(dir+'obj-sky.fts')
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2') & Npol=sxpar(h,'NAXIS3')
Nexp=sxpar(h,'NAXIS4') & Ntar=sxpar(h,'NAXIS5')
d_wave=sxpar(h,'CDELT1')
wave=findgen(Nx)*d_wave+sxpar(h,'CRVAL1')
obj=readfits(dir+'obj-sky.fts')
; sky=readfits(dir+'obj_lin.fts')
;neon=readfits(dir+'neon_lin.fts')
;atm_abs=readfits(dir+'atm_abs.fit')
;atm_abs=shift_s(1+(atm_abs-1)*1.1,-1)
;atm_abs=1
;print,Nx,Ny,Npol,Nexp
;END

;центрирование и нормировка спектров
;y=dindgen(Ny)
;yc=Ny/2 & wc=10 & dwc=3
;w=100
;norm=fltarr(Npol,Nexp) & shi_y=fltarr(Nexp)
;for j=0,Npol-1 do begin
;for k=0,Nexp-1 do begin
; yslice=total(obj(Nx/2-w:Nx/2+w,*,j,k,0),1)
;max_value=max(yslice(yc-wc:yc+wc),Nmax)
;yrng=findgen(2*dwc+1)+yc-wc+Nmax-dwc
;f=goodpoly(y(yrng),yslice(yrng),2,3,fit)
;dy=-f(1)/f(2)/2-yc & shi_y(k)=dy
;obj(*,*,j,k,0)=shift_image(obj(*,*,j,k,0),0,-dy)
;yslice=median(yslice,5)
;norm(j,k)=total(yslice(wc:Ny-wc))
;endfor
;robomean,norm(j,*),3,0.5,avg_norm
;norm(j,*)=norm(j,*)/avg_norm
;if ind_print eq 1 then begin
;print,'Npol=',j
;print,'shift along slit',shi_y,format='(A17,'+string(Nexp,format='(I2)')+'F10.3)'
;print,'transmission',norm(j,*),format='(A15,'+string(Nexp,format='(I2)')+'F10.3)'
;endif
;endfor
;spectra_pol=fltarr(Nx,Ny,Npol)
;;исправление вариции прозрачности и чистка космических частиц
;for k=0,Npol-1 do begin
;for x=0,Nx-1 do begin
;for y=0,Ny-1 do begin
;spectra_pol(x,y,k)=median(obj(x,y,k,*)/norm(k,*))*Nexp
;endfor & endfor
;endfor
;spectra=total(spectra_pol,3)*Npol
;print,size(spectra)
;END
;исправление спектральной чувствительности и атмосферного поглащения
;atm_abs=readfits(dir+'atm_abs.fit')
;atm_abs=shift_s(1+(atm_abs-1)*1.1,-1)
;atm_abs=1
;sent=readfits(dir+'sent.fts') & N=N_elements(sent)
;sent=1
;for k=0,Ny-1 do begin
;spectra(R,k)=spectra(R,k)/sent(R)/atm_abs(R)
;endfor
;spectra=spectra(R,*)
;sxaddpar,h,'CRVAL1',wave(R(0))
;writefits,dir+'total_spectra.fts',spectra,h
window,2
!P.multi=[0,1,1]
wy=4
plot,[wave(0),wave(Nx-1)],[0,1]*1.45E4,xst=1,yst=1,/nodata
for k=0,Nexp-1 do oplot,wave,total(obj(*,Ny/2-wy:Ny/2+wy,3,k,0),2)
END
