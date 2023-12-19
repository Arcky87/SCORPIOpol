;test dispersion WOLL-2
function etalon_spectra,spectra,GAMMA=gamma
;формирование еталонного спектра неона, гамма-коррекция и вычитание фона
if not(keyword_set(gamma)) then G=0 else G=gamma
a=size(spectra) & Nx=a(1)
x=findgen(Nx)
R=where(spectra lt 1, ind)
if ind gt 1 then spectra(R)=1
spectra=spectra^G
fi_peak,x ,spectra,0,ipix,xpk,ypk,bkpk,ipk
fon=INTERPOL(bkpk,xpk,x)
fon=LOWESS(x,fon,Nx/4,2,2)
spectra=spectra-fon
R=where(spectra lt 0)
if ind gt 1 then spectra(R)=0
return,spectra
end
;**************************************************
LOGFILE=DIALOG_PICKFILE(/read,path='h:\red_data.pol\LOGS\',FILTER='*.txt')

dir=def_wdir(LOGFILE)

;goto,cont



cube_neon=readfits(dir+'neon.fts',h)

Nx=sxpar(h,'NAXIS1')
Ny=sxpar(h,'NAXIS2')
Npol=4
Ndeg=3
;построение дисперсионной кривой
;чтение таблицы линий
Ntab=numlines('h:\WOLLASTON-2.lib\'+'etalon.txt')
tab=fltarr(5,Ntab)
openr,1,'h:\WOLLASTON-2.lib\'+'etalon.txt'
readf,1,tab
close,1

d_wave=tab(0,0) & min_wave=tab(1,0) & max_wave=tab(2,0)
tab=tab(*,1:Ntab-1) & Ntab=Ntab-1

ident_table=fltarr(2,Ntab)
ident_table(0,*)=tab(0,*)


Disp=dblarr(Ny,Ndeg+1,Npol)
FOR k=0,Npol-1 DO BEGIN
ident_table(1,*)=tab(k+1,*)

DISP(*,*,k)=dispersion_WOLL2(cube_neon(*,*,k),ident_table,N_DEG=Ndeg  ,plot=k+1);,,TRESH=tresh,PLOT=plot)
ENDFOR
writefits,dir+'disp.fts',disp

END
cont:
disp=readfits(dir+'disp.fts')
lambda_0=3800 & d_lambda=2 & Nlin=2150 ; VPHG940@600 + GS-11


;линеаризация  калибровочных спектров

lambda_0=3800 & d_lambda=2 & Nlin=2150 ; VPHG940@600 + GS-11
type=['eta','flat','neon']

for k=0,2 do begin
wait,1
cube_ini=readfits(dir+type(k)+'.fts',h)
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2')
Npol=4
Window,2,xsize=Nlin/2,ysize=(Ny+1)*4,title=dir
cube_lin=fltarr(Nlin,Ny,Npol)
for j=0,Npol-1 do begin
ima=linerisation_WOLL2(cube_ini(*,*,j),DISP(*,*,j),PARAM=[lambda_0,d_lambda,Nlin])
robomean,congrid(ima,Nlin,Ny),3,0.5,avg_ima,rms_ima
tv,255-bytscl(congrid(ima,Nlin/2,Ny),avg_ima-rms_ima*3,avg_ima+rms_ima*3),0,(Ny+1)*j
cube_lin(*,*,j)=ima
endfor
sxaddpar,h,'CRVAL1',lambda_0
sxaddpar,h,'CDELT1',d_lambda
writefits,dir+type(k)+'_lin.fts',cube_lin,h

ENDFOR



;линеаризация  спектра объекта
cube_obj=readfits(dir+'obj.fts',h)

Npol=4
Nexp=sxpar(h,'NAXIS4')
print,Nexp
Ncub=3
cube_lin=fltarr(Nlin,Ny,Npol,Nexp,Ncub)


for i=0,Ncub-1 do begin

for j=0,Nexp-1 do begin
if total(cube_obj(*,*,*,j,i)) ne 0 then begin
;Window,2,xsize=Nlin/2,ysize=(Ny+1)*4,title=dir+'   cube '+string(i,format='(I2)')+'  exp '+string(j,format='(I3)')
for k=0,Npol-1 do begin
ima=linerisation_WOLL2(cube_obj(*,*,k,j,i),DISP(*,*,k),PARAM=[lambda_0,d_lambda,Nlin])
;robomean,ima(*,Ny/2-10:Ny/2+10),3,0.5,avg_ima,rms_ima
;tv,255-bytscl(congrid(ima,Nlin/2,Ny),avg_ima-rms_ima*3,avg_ima+rms_ima*10),0,(Ny+1)*k
cube_lin(*,*,k,j,i)=ima
endfor
print,'linearization', i,j
endif
;wait,1
endfor
endfor
sxaddpar,h,'CRVAL1',lambda_0
sxaddpar,h,'CDELT1',d_lambda
writefits,dir+'obj_lin.fts',cube_lin,h
;исправление плоского поля
flat=readfits(dir+'flat_lin.fts',h)
a=size(flat) & Nx=a(1)  & Ny=a(2)  & Npol=a(3)
;goto,cont
map=fltarr(Nx,Ny)
norm=fltarr(Nx,Npol)

cross=fltarr(Ny,Npol)
w=5
;формирование нормировки плоского поля

for k=0,Npol-1 do begin
map(*,*)=flat(*,*,k)
for ky=0,Ny-1 do map(*,ky)=smooth(map(*,ky),50,/edge_truncate)
norm(*,k)=total(map(*,Ny/2-w:Ny/2),2)/(2*w+1)
cross(*,k)=total(map(Nx/2-w:Nx/2,*),1)/(2*w+1)

endfor


x=findgen(Nx)
for j=0,Npol-1 do norm(*,j)=LOWESS(x,norm(*,j),Nx/32,2,2)
norm=total(norm,2)/4
for k=0,Npol-1 do begin
for ky=0,Ny-1 do begin
flat(*,ky,k)=flat(*,ky,k)/norm
flat(*,ky,k)=LOWESS(x,flat(*,ky,k),Nx/16,2,2)
print,'create flat',k,ky
wait,0.02
end
endfor

writefits,dir+'flat_norm.fts',flat

flat=readfits(dir+'flat_norm.fts')
w=3
norm=fltarr(Npol)
for k=0,Npol-1  do begin
norm=total(flat(Nx/2-w:Nx/2+w,Ny/2-w:Ny/2+w,k))/(2*w+1)/(2*w+1)
flat(*,*,k)=flat(*,*,k)/norm
endfor

writefits,dir+'avg_flat.fts',flat,h
subsky_WOLL2,LOGFILE
spectra_WOLL2,LOGFILE,/plot,Wy=12;,/order
stocks_WOLL2,LOGFILE,WIDTH=2,xrange=xrng,AmpP=5
end