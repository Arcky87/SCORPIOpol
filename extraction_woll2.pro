pro extraction_WOLL2,dir,PLOT=plot,YC=yc,WY=wy,ATM_ABS=atm_abs,AMPL=ampl,SHIFT_X=shift_x
!P.charsize=1
print,Yc,Wy
; extraction spectra WOLL2
if not(keyword_set(shift_x)) then shift_x=0
if not(keyword_set(wy)) then wy=20
if not(keyword_set(ampl)) then ampl=1
cube=readfits(dir+'obj-sky.fts',h)
a=size(cube)
if not(keyword_set(yc)) then yc=a(2)/2

Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2')
Npol=sxpar(h,'NAXIS3') & Ntarget=sxpar(h,'NAXIS5')
Nexp=intarr(Ntarget) & wx=50  & Name=strarr(Ntarget) & PA=fltarr(Ntarget)

for k=0,Ntarget-1 do begin
Nexp(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))
Name(k)=sxpar(h,'NAME'+string(k+1,format='(I1)'))
PA(k)=sxpar(h,'PA'+string(k+1,format='(I1)'))
endfor

spectra=fltarr(Nx,Npol,Nexp(0),Ntarget)

;исправление атмосферного поглощения
atm_index=0
if keyword_set(atm_abs) then begin
atm=readfits(dir+'atm_abs.fit')
atm=(atm-1)*ampl+1
for j=0,Npol-1 do atm(*,j)=shift_s(atm(*,j),shift_x)
atm_index=1
print,'correction atmosphere'
endif
mean_spectra=fltarr(Nx,Npol,Ntarget)


for j=0,Ntarget-1 do begin
for i=0,Nexp(j)-1 do begin

!P.multi=[0,1,1]
for k=0,Npol-1 do begin
spectra(*,k,i,j)=total(cube(*,yc-wy/2:yc+wy/2,k,i,j),2)
if atm_index eq 1 then spectra(*,k,i,j)=spectra(*,k,i,j)/atm(*,k)
endfor & endfor & endfor

;проверка совместимости спектров
norm=fltarr(Npol,Nexp(0),Ntarget)
for k=0,Ntarget-1 do begin
	for j=0,Npol-1 do begin
		for i=0,Nexp(k)-1 do begin
			robomean,spectra(Nx/2-wx:Nx/2+wx,j,i,k),3,0.5,avg_spectra,rms_spectra
			norm(j,i,k)=avg_spectra
			spectra(*,j,i,k)=spectra(*,j,i,k)/norm(j,i,k)
		endfor
	spectra(*,j,*,k)=spectra(*,j,*,k)*max(norm(j,*,k))
	;for x=0,Nx-1 do mean_spectra(x,j,k)=MEDIAN(spectra(x,j,0:Nexp(k)-1,k))

for x=0,Nx-1 do begin
	if total(spectra(x,j,0:Nexp(k)-1,k)) eq 0 then begin
	mean_spectra(x,j,k)=MEDIAN(spectra(x,j,0:Nexp(k)-1,k))
	endif else begin
	robomean,[spectra(x,j,0:Nexp(k)-1,k)],1,0.5,av
	mean_spectra(x,j,k)=av
	endelse
endfor

	endfor
endfor
;window,2
;!P.multi=[0,1,4]
;for k=0,Npol-1 do plot,norm(k,*,0),PSYM=6




title_target=['object','unpolarized star','polarized star']
angle=['0','90','45','135']+' deg'
window,3,xsize=1800,ysize=600
!P.multi=[0,1,1]
for k=0,Ntarget-1 do begin
;y_max=max(median(mean_spectra(*,*,k),3))
robomean,mean_spectra(*,*,k),3,0.5,avg_val
if k gt 0 then y_max=2*avg_val ELSE y_max=5*avg_val
plot,[0,1],[0,1],/nodata,position=[0.005+0.33*k,0.01,0.005+0.33*(k+1),0.97],/norm, noerase=k,$
   TICKLEN=0,xcharsize=1e-5,ycharsize=1e-5,$
   title=title_target(k)+'  '+name(k)+' '+'Nexp='+string(Nexp(k),format='(I2)')
for j=0,Npol-1 do begin
plot,[0,Nx-1],[0,y_max],xst=1,yst=1,$
	position=[0.005+0.33*k,0.01+0.24*j,0.005+0.33*(k+1),0.01+0.24*(j+1)],$
	/norm,/noerase,/nodata,TICKLEN=0,charsize=1e-5
for i=0,Nexp(k)-1 do oplot,spectra(*,j,i,k),color=3e5
oplot,mean_spectra(*,j,k)
xyouts,Nx*0.85,y_max*0.85,angle(j)
endfor
endfor
;вычисление параметров Стокса
mode=0
S=fltarr(Nx,4,Nexp(0),Ntarget)
		for j=0,Ntarget-1 do begin
	for k=0,Nexp(j)-1 do begin
S(*,0,k,j)=total(spectra(*,*,k,j),2)
S(*,2,k,j)=(spectra(*,1,k,j)-spectra(*,0,k,j))/(spectra(*,1,k,j)+spectra(*,0,k,j))
S(*,1,k,j)=(spectra(*,3,k,j)-spectra(*,2,k,j))/(spectra(*,3,k,j)+spectra(*,2,k,j))
if mode eq 1 then begin
Q_tmp=S(*,2,k,j)
U_tmp=-S(*,1,k,j)
S(*,1,k,j)=Q_tmp
S(*,2,k,j)=U_tmp
endif

	endfor
		endfor
;проверка типа эталона



;Запись исходных параметров Стокса
print,'file_test',FILE_TEST(dir+'stoks.fit')
writefits,dir+'stoks.fit',S,h

writefits,dir+'spectra.fit',spectra,h
writefits,dir+'avg_spectra.fit',mean_spectra,h
END

wdir='h:\red_data.pol\E1841+643_140324\'

wdir='h:\red_data.pol\Sy1\NGC3516_170131\'
;wdir='h:\red_data.pol\Sy1\Mkn79_151207\'
wdir='e:\sbs1419+538_190216\'
extraction_WOLL2,wdir,/plot,yc=41,Wy=40,/atm_abs,ampl=1 ,shift_x=0

end
