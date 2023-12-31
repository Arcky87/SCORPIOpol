;check_stoks_WOLL2
;�������� �������� ����� ����������� �������� ������
;goto,fin
pro stoks_WOLL2,wdir
;;�������� ����� ���� ����
;neon=readfits(wdir+'neon_lin.fts',h)
;
;Nx=sxpar(h,'NAXIS1')
;Ny=sxpar(h,'NAXIS2')
;yc=Ny/2  & wy=20
;vector=total(neon(*,yc-wy:yc+wy,*),2)
;x_shift=fltarr(4)
;M=10 & xcross=findgen(2*M+1)-M
;for k=0,3 do begin
;ycross=cross_norm(vector(*,k),vector(*,0),M)
;gau=gaussfit(xcross,ycross,G)
;x_shift(k)=G(1)
;endfor
;print,x_shift,format='(4F7.2)'

;vector(*,k)=shift_s(vector(*,k),-x_shift(k))


;������������ �������� �������
;

;s=readfits(wdir+'avg_spectra.fit',h)

spectra=readfits(wdir+'spectra.fit',hs)
Nx=sxpar(hs,'NAXIS1') & Npol=sxpar(hs,'NAXIS2') & Ntarget=sxpar(hs,'NAXIS4')
Nexp=intarr(Ntarget) & wx=50  & Name=strarr(Ntarget) & PA=fltarr(Ntarget)
for k=0,Ntarget-1 do begin
Nexp(k)=sxpar(hs,'NUMEXP'+string(k+1,format='(I1)'))
Name(k)=sxpar(hs,'NAME'+string(k+1,format='(I1)'))
PA(k)=sxpar(hs,'PA'+string(k+1,format='(I1)'))
endfor
x=findgen(Nx)
wave=findgen(Nx)*sxpar(hs,'CDELT1') +sxpar(hs,'CRVAL1')
;star_name=name(2)
;PA_slit=PA(2)
mean_spectra=fltarr(Nx,Npol,Ntarget)
;�������� ������������� ��������
norm=fltarr(Npol,max(Nexp),Ntarget)
for k=0,Ntarget-1 do begin
for j=0,Npol-1 do begin
for i=0,Nexp(k)-1 do begin
robomean,spectra(Nx/2-wx:Nx/2+wx,j,i,k),3,0.5,avg_spectra,rms_spectra
norm(j,i,k)=avg_spectra
spectra(*,j,i,k)=spectra(*,j,i,k)/norm(j,i,k)
;spectra(*,j,i,k)=shift_s(spectra(*,j,i,k),-x_shift(j))
endfor
spectra(*,j,*,k)=spectra(*,j,*,k)*max(norm(j,*,k))
for x=0,Nx-1 do mean_spectra(x,j,k)=MEDIAN(spectra(x,j,0:Nexp(k)-1,k))
endfor
endfor

title_target=['object','unpolarized star','polarized star']
angle=['0','90','45','135']+' deg'
set_plot,'PS'
device,file=wdir+'spectra.ps',/landscape
;window,4,xsize=1200,ysize=600
!P.multi=[0,1,1]
for k=0,Ntarget-1 do begin
y_max=max(mean_spectra(*,*,k))
plot,[0,1],[0,1],/nodata,position=[0.005+0.33*k,0.01,0.005+0.33*(k+1),0.97],/norm, noerase=k,$
   TICKLEN=0,xcharsize=1e-5,ycharsize=1e-5,$
   title=title_target(k)+'  '+name(k)+' '+'Nexp='+string(Nexp(k),format='(I2)')+' PA='+STRING(PA(k),format='(F7.1)'),charsize=0.6
for j=0,Npol-1 do begin
plot,[0,Nx-1],[0,y_max],xst=1,yst=1,$
	position=[0.005+0.33*k,0.01+0.24*j,0.005+0.33*(k+1),0.01+0.24*(j+1)],$
	/norm,/noerase,/nodata,TICKLEN=0,charsize=1e-5
for i=0,Nexp(k)-1 do oplot,spectra(*,j,i,k),color=120
oplot,mean_spectra(*,j,k)
xyouts,Nx*0.85,y_max*0.85,angle(j),charsize=0.6
endfor
endfor
xyouts,0,-0.01,wdir,/norm,charsize=0.6
device,/close
set_plot,'WIN'
;����������� �������������
avg_D=fltarr(2,Nexp(0),Ntarget)  & rms_D=avg_D
wave_c=5500 & d_wave=250
R=where(wave gt wave_c-d_wave and wave lt wave_c+d_wave,ind)
p=0 & t=0

		for t=0,Ntarget-1 do begin
	for j=0,1 do begin
for  k=0,Nexp(t)-1 do begin
robomean,spectra(R,2*j+1,k,t)/spectra(R,2*j,k,t),1,0.5,avg_val,rms_val
avg_D(j,k,t)=avg_val & rms_D(j,k,t)=rms_val
endfor
		robomean,avg_D(j,0:Nexp(t)-1,t),3,0.5,mean
		avg_D(j,0:Nexp(t)-1,t)=avg_D(j,0:Nexp(t)-1,t)/mean
		rms_D(j,0:Nexp(t)-1,t)=rms_D(j,0:Nexp(t)-1,t)/mean
	endfor
		endfor
set_plot,'PS'
device,file=wdir+'depolarization.ps',/landscape

!P.multi=[0,1,1]
CS=1
ytlt=['FLUX(90)/FLUX(0)','FLUX(135)/FLUX(45)']
for t=0,Ntarget-1 do begin
	for i=0,1 do begin

if t+i gt 0 then N=1  ELSE N=0
if t gt 0 then CSY=1e-5 ELSE  CSY=CS
if i gt 0 then CSX=1e-5 ELSE  CSX=CS
if i gt 0 then tlt=Name(t) ELSE tlt=''
plot,avg_D(i,0:Nexp(t)-1,t),psym=6,charsize=CS,ycharsize=CSY,xcharsize=CSX,title=tlt,$
xrange=[-0.5,Nexp(t)-0.5],xst=1,yrange=[0.95,1.05],yst=1,ytitle=ytlt(i),$
position=[0.09+t*0.3,0.07+0.45*i,0.09+(1+t)*0.3,0.07+(i+1)*0.45],/norm,noerase=N
oploterr,avg_D(i,0:Nexp(t)-1,t),rms_D(i,0:Nexp(t)-1,t)
oplot,[-1,Nexp(t)],[1,1],linestyle=2
endfor
	endfor
xyouts,0,0,wdir,/norm,charsize=1.5
device,/close
set_plot,'WIN'

;name=strarr(3) & PA=fltarr(3)
;for j=0,2 do begin
;name(j)=sxpar(h,'NAME'+string(j+1,format='(I1)'))


;endfor

;window,1,xsize=600,ysize=427
;!P.multi=[0,1,2]
;x=findgen(Nx)
;
;;cut=4200
;;& deg=5
;R01=s(*,0,1)/s(*,1,1)
;R23=s(*,2,1)/s(*,3,1)
;plot,R01,xst=1,charsize=1,color=1e5,yrange=[0,3]
;fit=LOWESS(x(beg:cut),R01(beg:cut),Nx/2 ,3,1)
;R01=INTERPOL( fit,x(beg:cut),x)
;oplot,x,R01
;plot,R23,xst=1,charsize=1,color=1e5,yrange=[0,3]
;fit=LOWESS(x(beg:cut),R23(beg:cut),Nx/2 ,3,1)
;R23=INTERPOL( fit,x(beg:cut),x)
;oplot,x,R23
;
;;
;stoks=fltarr(Nx,4,3)
cube_stoks=fltarr(Nx,4,Nexp(0),3)
;Nx=2291
;s=s(200:Nx-200,*,*)
;a=size(s) &
;Nx=a(1) &

Q=fltarr(Nx,3) & U=Q
;***********************
mode=slope_flat(wdir)
;**********************
;avg_D=fltarr(2,Nexp(0),Ntarget)+1
for k=0,Ntarget-1 do begin
;stoks(*,0,k)=s(*,0,k)+R01*s(*,1,k)+s(*,2,k)+R23*s(*,3,k)
;stoks(*,1,k)=(s(*,0,k)-R01*s(*,1,k))/(s(*,0,k)+R01*s(*,1,k))
;stoks(*,2,k)=(s(*,2,k)-R23*s(*,3,k))/(s(*,2,k)+R23*s(*,3,k))
for j=0,Nexp(k)-1 do begin
cube_stoks(*,0,j,k)= spectra(*,0,j,k)+spectra(*,1,j,k)+spectra(*,2,j,k)+spectra(*,3,j,k)
;cube_stoks(*,1,j,k)=(spectra(*,0,j,k)-avg_D(0,j,k)*spectra(*,1,j,k))/$
;					(spectra(*,0,j,k)+avg_D(0,j,k)*spectra(*,1,j,k))
cube_stoks(*,1,j,k)=(avg_D(0,j,k)*spectra(*,0,j,k)-spectra(*,1,j,k))/$
					(avg_D(0,j,k)*spectra(*,0,j,k)+spectra(*,1,j,k))

;cube_stoks(*,2,j,k)=(spectra(*,2,j,k)-avg_D(1,j,k)*spectra(*,3,j,k))/$
;				(spectra(*,2,j,k)+avg_D(1,j,k)*spectra(*,3,j,k))
cube_stoks(*,2,j,k)=(avg_D(1,j,k)*spectra(*,2,j,k)-spectra(*,3,j,k))/$
       			(avg_D(1,j,k)*spectra(*,2,j,k)+spectra(*,3,j,k))
endfor
endfor
;Q_null=Smooth(stoks(*,1,1),wx,/edge_truncate)
;U_null=SMooth(stoks(*,2,1),wx,/edge_truncate)
;Q_null=0 & U_null=0

for k=0,Ntarget-1 do begin
;stoks(*,1,k)=stoks(*,1,k)-Q_null
;stoks(*,2,k)=stoks(*,2,k)-U_null
for j=0,Nexp(k)-1 do begin
cube_stoks(*,1,j,k)=cube_stoks(*,1,j,k);-Q_null
cube_stoks(*,2,j,k)=cube_stoks(*,2,j,k);-U_null

if mode eq 1 then begin
Q_tmp=cube_stoks(*,2,j,k)
U_tmp=-cube_stoks(*,1,j,k)
cube_stoks(*,1,j,k)=Q_tmp
cube_stoks(*,2,j,k)=U_tmp
endif
endfor
endfor
;if mode eq 1 then begin
;Q_tmp=stoks(*,1,*) & U_tmp=stoks(*,2,*)
;stoks(*,1,*)=U_tmp & stoks(*,2,*)=-Q_tmp
;endif
;writefits,wdir+'avg_stoks.fit',stoks,hs
writefits,wdir+'stoks.fit',cube_stoks,hs




end
log_dir='d:\SCORPIO\red_data.pol\'
;log_dir='h:\red_data.pol\AGN\'
;LOGFILE=DIALOG_PICKFILE(/read,path=log_dir+'LOGS\',FILTER='*.txt')
LOGFILE=log_dir+'LOGS\'+'NGC4151_160307.txt'
	name_out=str_sep(FILE_BASENAME(LOGFILE),'.')
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
	wdir=log_dir+wdir(N_elements(wdir)-2)+'\'
stoks_WOLL2,wdir
	end