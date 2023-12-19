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

wdir=def_wdir(LOGFILE)


neon=readfits(wdir+'neon.fts',h)
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2')
Npol=4 & Ndeg=3
;построение дисперсионной кривой
;чтение таблицы линий
Ntab=numlines(wdir+'etalon.txt')
tab=fltarr(5,Ntab)
openr,1,wdir+'etalon.txt'
readf,1,tab
close,1
d_wave=tab(0,0) & min_wave=tab(1,0) & max_wave=tab(2,0)
tab=tab(*,1:Ntab-1) & Ntab=Ntab-1
ident_table=fltarr(2,Ntab)
ident_table(0,*)=tab(0,*)
Disp=dblarr(Ny,Ndeg+1,Npol)
FOR k=0,Npol-1 DO BEGIN
ident_table(1,*)=tab(k+1,*)
DISP(*,*,k)=dispersion_WOLL2(neon(*,*,k),ident_table,N_DEG=Ndeg ,plot=k+1)
ENDFOR
;END
;lambda_0=3600 & d_lambda=2 & Nlin=2400
lambda_0=3800 & d_lambda=2 & Nlin=2150 ; VPHG940@600 + GS-11
;lambda_0=5750 & d_lambda=2 & Nlin=1950
j=0
window,2,xsize=1200,ysize=81*4
for j=0,3 do begin
ima=linerisation_WOLL2(neon(*,*,j),DISP(*,*,j),PARAM=[min_wave,d_wave,FIX((max_wave-min_wave)/d_wave)])
tv,255-bytscl(congrid(ima,1200,80),0,1000),0,81*j
print,size(ima)
endfor


;линеаризация  калибровочных спектров
;lambda_0=3600 & d_lambda=2 & Nlin=2400 ; VPHG940@600
lambda_0=5750 & d_lambda=2 & Nlin=1950 ;VPHG1026@735
;lambda_0=3600 & d_lambda=1.5 & Nlin=2400 ; VPHG1200@540
;lambda_0=4200 & d_lambda=2 & Nlin=1900 ; VPHG940@600 + GS-11
;goto,cont
;lambda_0=3800 & d_lambda=2 & Nlin=2150 ; VPHG940@600 + GS-11
type=['eta','flat','neon']

for k=0,2 do begin
wait,1
cube_ini=readfits(wdir+type(k)+'.fts',h)
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
writefits,wdir+type(k)+'_lin.fts',cube_lin,h

ENDFOR


cont:
;линеаризация  спектра объекта
cube_obj=readfits(wdir+'obj.fts',h)
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
endif
;wait,1
endfor
endfor
sxaddpar,h,'CRVAL1',lambda_0
sxaddpar,h,'CDELT1',d_lambda
writefits,wdir+'obj_lin.fts',cube_lin,h
;output
W_DIR=sxpar(read_table(LOGFILE),'w_dir')
LOADFILE,dir=w_dir,ysz=500,xsz=1130

ViewPol_2
end