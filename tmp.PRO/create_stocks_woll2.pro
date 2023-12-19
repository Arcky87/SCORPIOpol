;test WOLLASTON2

pro create_stocks_WOLL2,wdir,XRANGE=xrange,WIDTH=width,ERR=err,AmpP=AmpP,WAVE_C=wave_c,BOBO=bobo,ATM_ABS=atm_abs,ISM=ism

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
date=date+' JD'+string(JD,format='(I5)')
titl=['object','unpolarized star','polarized star']

lambda_0= sxpar(h,'CRVAL1')  & d_lambda= sxpar(h,'CDELT1')
lambda=findgen(Nx)*d_lambda+lambda_0
;goto,cont1
;чтение спектральной чувствительности
;sent=readfits(wdir+'sent.fts')
sent=1
if not(keyword_set(wave_c)) then wave_c=lambda(Nx/2)
print,wave_c
if not(keyword_set(xrange)) then xrng=[lambda(0),lambda(Nx-1)] else xrng=xrange
if keyword_set(atm_abs) then begin
;исправление атмосферного поглощения
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




;контроль атмосферы
j=0
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
set_plot,'PS'
device,file=wdir+'depolarization.ps'
	!P.multi=[0,1,2]
		plot,D(0,*),psym=6,yrange=[0.95,1.05],yst=1,$
			title='wollaston 0-90 deg',$
			xtitle='Number of exposure',$
			ytitle='depolarization'
		oploterr,D(0,*),err_D(0,*),psym=6
		oplot,[0,Nexp(j)],[1,1],linestyle=2
			plot,D(1,*),psym=6,yrange=[0.95,1.05],yst=1,$
			title='wollaston 45-135 deg',$
			xtitle='Number of exposure',$
			ytitle='depolarization'
			oploterr,D(1,*),err_D(1,*),psym=6
			oplot,[0,Nexp(j)],[1,1],linestyle=2
		subtitle=wdir
	device,/close
	set_plot,'WIN'






;вычисление параметров Стокса
S=fltarr(Nx,4,Nexp(0),Ncube)
		for j=0,Ncube-1 do begin
	for k=0,Nexp(j)-1 do begin
sent=1
S(*,0,k,j)=total(spectra(*,*,k,j),2)/sent
;S(*,2,k,j)=(spectra(*,1,k,j)-spectra(*,0,k,j))/(spectra(*,1,k,j)+spectra(*,0,k,j))
;S(*,1,k,j)=(spectra(*,3,k,j)-spectra(*,2,k,j))/(spectra(*,3,k,j)+spectra(*,2,k,j))

S(*,2,k,j)=(spectra(*,1,k,j)-spectra(*,0,k,j))/(spectra(*,1,k,j)+spectra(*,0,k,j))
S(*,1,k,j)=(spectra(*,3,k,j)-spectra(*,2,k,j))/(spectra(*,3,k,j)+spectra(*,2,k,j))

	endfor
		endfor

;Запись исходных параметров Стокса
print,'file_test',FILE_TEST(wdir+'stoks.fit')
writefits,wdir+'stoks.fit',S,h






END

path='h:\red_data.pol\'
xrng=[4400,7500]; & z_obj=0.0561 & Wc=5500
xrng=[4200,8000]; & z_obj=0.0561 & Wc=5500

;path='h:\red_data.pol\AGN\'
LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path=path+'LOGS')
;LOGFILE=DIALOG_PICKFILE(/read,path='h:\red_data.pol\LOGS\',FILTER='*.txt')
;LOGFILE='h:\red_data.pol\LOGS\3C273_140325.txt\
;LOGFILE='h:\red_data.pol\LOGS\3C390_140324.txt\
;stocks_WOLL2,LOGFILE,xrange=[6500,8500],WIDTH=2,AmpP=5,/atm_abs,/bobo,wave_c=7200	;for 3C273.3
;stocks_WOLL2,LOGFILE,z=0.0561,xrange=xrng,WIDTH=10,AmpP=4,wave_c=Wc,/bobo,/atm_abs;,ISM=[0.64,-0.51] 		;for 3C390.3
tmp=str_sep(LOGFILE,'\')
path=tmp(0)
for j=1,N_elements(tmp)-3 do begin
path=path+'\'+tmp(j)
endfor
path=path+'\'
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
wdir=path+wdir(N_elements(wdir)-2)+'\'
print,wdir
create_stocks_WOLL2,wdir,WIDTH=2,AmpP=5,wave_c=5500;,/atm;,ISM=[0 ,-4.5],/bobo;ISM=[0.64,-0.51]; ,/bobo; ;,z
end