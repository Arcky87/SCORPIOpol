pro extraction_WOLL2,dir,PLOT=plot,YC=yc,WY=wy,ATM_ABS=atm_abs
; extraction spectra WOLL2
;if not(keyword_set(wy)) then wy=10

cube=readfits(dir+'obj-sky.fts',h)
a=size(cube)
if not(keyword_set(yc)) then yc=a(2)/2
Nx=a(1) & Ny=a(2) & Npol=a(3) & Nmax=a(4)  & Ncube=a(5)
y=findgen(Ny)
Nexp=fltarr(Ncube)
for k=0,Ncube-1 do Nexp(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))
spectra=fltarr(Nx,Npol,Nmax,Ncube)
titl=['object','unpolrized star','polarized star']
angle=['0 deg','90 deg','45 deg','135 deg']
;����������� ������������ ����������
if keyword_set(atm_abs) then begin
atm=readfits(dir+'atm_abs.fit')
print,size(atm)

 for j=0,Ncube-1 do begin
 for k=0,Nexp(j)-1 do begin
for i=0,Npol-1 do begin
spectra(*,i,k,j)=spectra(*,i,k,j)/atm(*,i)

endfor & endfor & endfor
endif



p=0
tmp_1=fltarr(Nx,Npol)  & tmp_2=tmp_1
if keyword_set(plot) then p=1
for j=0,Ncube-1 do begin
for i=0,Nexp(j)-1 do begin
if p eq 1 then WINDOW,3,xsize=1000,ysize=400,title=dir+'  '+titl(j)+'   exposure'+string(i+1)
!P.multi=[0,1,1]
for k=0,Npol-1 do begin
spectra(*,k,i,j)=total(cube(*,yc-wy:yc+wy,k,i,j),2)
if keyword_set(order) then begin
tmp_1(*,k)=total(cube(*,yc-wy:yc+wy,k,i,j),2)
tmp_2(*,k)=total(cube(*,yc-2*wy:yc+2*wy,k,i,j),2)
endif
endfor
if keyword_set(order) then begin
tmp_1(*,1)=tmp_1(*,1)-(tmp_2(*,0)-tmp_1(*,0))
tmp_1(*,2)=tmp_1(*,2)-(tmp_2(*,3)-tmp_1(*,3))
spectra(*,*,i,j)=tmp_1(*,*)
endif
if p eq 1 then begin
for k=0,Npol-1 do begin
plot,spectra(*,k,i,j),xst=1,position=[0.05,0.25*k,1,0.25*(k+1)],noerase=k,$
	xcharsize=1e-5,yticklen=0.002,yst=1,yminor=2
xyouts,0.88,0.2+0.25*k ,angle(k),/norm
endfor
	endif

if p eq 1 then wait,0.5
endfor
endfor
;print,h
writefits,dir+'spectra.fit',spectra,h
END
;log_dir='h:\red_data.pol\TINATIN\'
;log_dir='h:\red_data.pol\AGN\'
;LOGFILE=DIALOG_PICKFILE(/read,path=log_dir+'LOGS\',FILTER='*.txt')
;wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
;wdir=log_dir+wdir(N_elements(wdir)-2)+'\'
;print,wdir
wdir='h:\red_data.pol\NGC7469_151210\'
extraction_WOLL2,wdir,/plot,Wy=20;,/atm_abs


end
;stocks_WOLL2,LOGFILE
;create_stocks_WOLL2,wdir,WIDTH=2,AmpP=5,wave_c=5500


stoks=readfits(wdir+'stoks.fit',h)
nX=SXPAR(H,'naxis1')  & nPOL=4
;atm_abs=readfits(wdir+'atm_abs.fit')
wave=findgen(Nx)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')

avg_stoks=fltarr(Nx,Npol)
for k=0,Npol-1 do begin
for x=0,Nx-1 do begin
avg_stoks(x,k)=median(stoks(x,k,*,1))
endfor
avg_stoks(*,k)=avg_stoks(*,k)
endfor

co=0.9
;atm_abs=(atm_abs-1)*co+1 &
;atm_abs=shift(atm_abs,1)
window,2,xsize=500,ysize=800
z=0.0172  & Ha=6562.8  & dw=500
R=where(wave gt (Ha-dw)*(1+z) and wave lt (Ha+dw)*(1+z))
!P.charsize=2
!P.multi=[0,1,5]
plot,wave(R),avg_stoks(R,0),xst=1
oplot,[1,1]*Ha*(1+z),[-1,1]*1e7,linestyle=2
 ;dU=0.1
; dQ=-0.1
plot,wave(R),avg_stoks(R,1)-dQ,xst=1,yrange=[-1,1]*0.1,yst=1
plot,wave(R),avg_stoks(R,2)-dU,xst=1,yrange=[-1,1]*0.1,yst=1
plot,wave(R),SQRT((avg_stoks(R,1)-dQ)^2+(avg_stoks(R,2)-dU)^2),xst=1,yrange=[0,0.1],yst=1
plot,wave(R),calc_atan(avg_stoks(R,1)-dQ,avg_stoks(R,2)-dU)/2 ,xst=1 ,yst=1
oplot,[1,1]*Ha*(1+z),[-1,1]*1e5,linestyle=2
END
