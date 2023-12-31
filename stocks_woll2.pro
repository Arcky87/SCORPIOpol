;test WOLLASTON2

pro stocks_WOLL2,LOGFILE,LEFT=left,RIGTH=rigth
tmp=str_sep(LOGFILE,'\')
path=tmp(0)
for j=1,N_elements(tmp)-3 do begin
path=path+'\'+tmp(j)
endfor
path=path+'\'
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
wdir=path+wdir(N_elements(wdir)-2)+'\'
;**********************

spectra=readfits(wdir+'spectra.fit',h)

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
;����������� ������������ ����������
;atm=readfits(wdir+'atm_abs.fit')
;;atm=(atm-1)*1.1+1
;a=size(atm)
;if a(0) eq 2 then atm=total(atm,2)/4
;a=size(atm)
;
;
;for j=0,Ncube-1 do begin
;for k=0,Nexp(j)-1 do begin
;for i=0,Npol-1 do begin
;;spectra(*,i,k,j)=spectra(*,i,k,j)/atm(*,i)
;
;if a(0) eq 2 then spectra(*,i,k,j)=spectra(*,i,k,j)/atm(*,i) ELSE $
;spectra(*,i,k,j)=spectra(*,i,k,j)/atm(*)
;endfor & endfor & endfor
;endif
j=1

;goto,cont
;							;�������� ���������
;	wave_a=6100  & d_wave=100
;	pos_c=(wave_a-lambda_0)/d_lambda
;	d_pos=d_wave/d_lambda
;	D=fltarr(2,Nexp(j))   & err_D=D
;
;	for i=0,1 do begin
;	for k=0,Nexp(j)-1 do begin
;	tmp=spectra(pos_c-d_pos:pos_c+d_pos,2*i,k)/spectra(pos_c-d_pos:pos_c+d_pos,2*i+1,k)
;	robomean,tmp,1,0.5,avg_D,rms_D
;
;	D(i,k)=avg_D  & err_D(i,k)=rms_D
;	endfor
;	robomean,D(i,*),3,0.5,avg_D
;	D(i,*)=D(i,*)/avg_D
;	err_D(i,*)=err_D(i,*)/avg_D
;	endfor
;	;print,D
;
;	window,1,xsize=600,ysize=500
;	!P.multi=[0,1,2]
;
;	plot,D(0,*),psym=6,yrange=[0.95,1.05],yst=1
;	oploterr,D(0,*),err_D(0,*),psym=6
;	oplot,[0,Nexp(j)],[1,1],linestyle=2
;	plot,D(1,*),psym=6,yrange=[0.95,1.05],yst=1
;	oploterr,D(1,*),err_D(1,*),psym=6
;	oplot,[0,Nexp(j)],[1,1],linestyle=2
;
;������������ ���������������� ��������������� �������
;window,1
;!P.multi=[0,1,2]
;x=findgen(Nx)
;beg=200
;cut=1800
;;cut=4200
;& deg=4
;R01=spectra(*,0,1)/spectra(*,1,1)
;R23=spectra(*,2,1)/spectra(*,3,1)
;plot,R01,xst=1,charsize=1,color=1e5,yrange=[0,3]
;fit=LOWESS(x(beg:cut),R01(beg:cut),Nx/2,3,1)
;R01=INTERPOL( fit,x(beg:cut),x)
;oplot,x,R01
;plot,R23,xst=1,charsize=1,color=1e5,yrange=[0,3]
;fit=LOWESS(x(beg:cut),R23(beg:cut),Nx/2,3,1)
;R23=INTERPOL( fit,x(beg:cut),x)
;oplot,x,R23

;���������� ���������� ������
S=fltarr(Nx,4,Nexp(0),Ncube)
		for j=0,Ncube-1 do begin
	for k=0,Nexp(j)-1 do begin

S(*,0,k,j)=total(spectra(*,*,k,j),2)
S(*,2,k,j)=(spectra(*,1,k,j)-spectra(*,0,k,j))/(spectra(*,1,k,j)+spectra(*,0,k,j))
S(*,1,k,j)=(spectra(*,3,k,j)-spectra(*,2,k,j))/(spectra(*,3,k,j)+spectra(*,2,k,j))

	endfor
		endfor

;������ �������� ���������� ������
print,'file_test',FILE_TEST(wdir+'stoks.fit')
writefits,wdir+'stoks.fit',S,h


END
path='h:\red_data.pol\AGN\'
;path='h:\red_data.pol\TINATIN\'
;path='k:\red_data.pol\AGN\'
xrng=[4400,7500]; & z_obj=0.0561 & Wc=5500
xrng=[4200,8000]; & z_obj=0.0561 & Wc=5500

;xrng=[7000,8200] & z_obj=0.158339 & Wc=7000
;xrng=[45000,7000]
;xrng=[7200,8200]
;path='h:\red_data.pol\AGN\'
LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path=path+'LOGS')


stocks_WOLL2,LOGFILE
end