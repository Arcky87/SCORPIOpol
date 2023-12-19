
;test WOLL-2
dir='h:\red_data.pol\Arp102b_131103\'
cubes=['s114201','s114204','s114203']
Nx=2326 & Ny=532 & dy=20
goto,cont
file_neon='s1142'+string([128,129,315,316,412,413],format='(I4.4)')+'.fts'
file_flat='s1142'+string([131,130,317,318,414,415],format='(I4.4)')+'.fts'
file_eta ='s1142'+string([132,133,316],format='(I4.4)')+'.fts'
Nz=3
cube=fltarr(Nx,Ny,Nz)
for j=0,Nz-1 do begin
tmp=readfits(dir+file_eta(j),h)
bias=def_bias(tmp)
cube(*,*,j)=tmp-bias
print,bias
endfor
ima=fltarr(Nx,Ny)
for x=0,Nx-1 do begin
for y=0,Ny-1 do begin
ima(x,y)=median(cube(x,y,*))
endfor
endfor
ima=shift(ima,0,dy)
writefits,dir+'eta_i.fts',ima,h

eta=readfits(dir+'eta_i.fts',h)
tra=create_traectory(eta,NP=12,WX=20,NDEG=2,X_beg=200,/plot)
writefits,dir+'traectory.fit',tra
;cont:
neon=readfits(dir+'neon_i.fts',h)
tra= readfits(dir+'traectory.fit',h)

H=131

Npol=4
cube_neon=fltarr(Nx,H,Npol)
for j=0,Npol-1 do cube_neon(*,*,j)=neon(*,j*H:(j+1)*H-1)
Ntra=12
cube_tra=fltarr(Nx,Ntra/4,4)
for j=0,Npol-1 do cube_tra(*,*,j)=tra(*,j*3:(j+1)*3-1)-H*j
tmp=fltarr(Nx,H)
Wtitle=['0','90','45','135']
for i=0,Npol-1 do begin
window,i+3,xsize=Nx/2,ysize=H,xpos=0,ypos=(H+35)*(Npol-1)-(H+35)*i,title=Wtitle(i)
tmp(*,*)=cube_neon(*,*,i)
tv,255-bytscl(congrid(tmp,Nx/2,H),0,200)
plot,[0,Nx-1],[0,H-1],xst=1,yst=1,/nodata,/noerase,position=[0,0,1,1],/norm
for j=0,2 do oplot,cube_tra(*,j,i),color=1e6
endfor
;cont:
; формирование куба объекта
hist='INITIAL FILES'
cub=['object','unpolarized star','polarized star']
Nobj=21
file_obj='s114201'+string(findgen(Nobj)+7,format='(I2.2)')+'.fts'
Ns_0=11
file_star_0='s114204'+string(findgen(Ns_0)+1,format='(I2.2)')+'.fts'
Ns=11
file_star='s114203'+string(findgen(Ns_0)+4,format='(I2.2)')+'.fts'
print,file_star
Nexp=[Nobj,Ns_0,Ns]
Ncube=3

header=create_header(Nx,Ny,Nobj,Ncube,1,0,1)

files=strarr(Nobj,Ncube)
files(0:Nobj-1,0)=file_obj
files(0:Ns_0-1,1)=file_star_0
files(0:Ns-1,2)=file_star
cube_obj=fltarr(Nx,Ny,Nobj,Ncube)
for k=0,Ncube-1 do begin
for E=0,Nexp(k)-1 do begin
tmp=FLOAT(readfits(dir+files(E,k),h))
PA=sxpar(h,'PARANGLE')-sxpar(h,'ROTANGLE')+132.5
sxaddpar,header,'name'+string(k+1,format='(I1)'),sxpar(h,'OBJECT')
sxaddpar,header,'RA'+string(k+1,format='(I1)'),sxpar(h,'RA')
sxaddpar,header,'DEC'+string(k+1,format='(I1)'),sxpar(h,'DEC')
sxaddpar,header,'A'+string(k+1,format='(I1)'),sxpar(h,'A')
sxaddpar,header,'Z'+string(k+1,format='(I1)'),sxpar(h,'Z')
sxaddpar,header,'PA'+string(k+1,format='(I1)'),PA
sxaddpar,header,'CUBE'+string(k+1,format='(I1)'),cubes(k)
if e eq 0 then begin
start=sxpar(h,'TIME-OBS')

sxaddpar,header,'START'+string(k+1,format='(I1)'),START
sxaddpar,header,'EXPTIME'+string(k+1,format='(I1)'),sxpar(h,'EXPTIME')*Nexp(k)
endif
	if k eq 0 then begin
sxaddpar,header,'DATE-OBS',sxpar(h,'DATE')
sxaddpar,header,'PROG-ID',sxpar(h,'PROG-ID')
sxaddpar,header,'AUTHOR',sxpar(h,'AUTHOR')
sxaddpar,header,'OBSERVER',sxpar(h,'OBSERVER')
sxaddpar,header,'DIR',dir
sxaddpar,header,'FOCUS',sxpar(h,'FOCUS')
sxaddpar,header,'BINNING',sxpar(h,'BINNING')
sxaddpar,header,'RATE',sxpar(h,'RATE')
sxaddpar,header,'GAIN',sxpar(h,'GAIN')
sxaddpar,header,'NODE',sxpar(h,'NODE')
sxaddpar,header,'IMSCALE',sxpar(h,'IMSCALE')
sxaddpar,header,'CAMFOCUS',sxpar(h,'CAMFOCUS')
sxaddpar,header,'COLFOCUS',sxpar(h,'COLFOCUS')
sxaddpar,header,'MODE',sxpar(h,'MODE')
sxaddpar,header,'DISPERSE',sxpar(h,'DISPERSE')
sxaddpar,header,'SLITWID',sxpar(h,'SLITWID')
sxaddpar,header,'SLITMASK',sxpar(h,'SLITMASK')
sxaddpar,header,'FILTERS',sxpar(h,'FILTERS')
sxaddpar,header,'FILTPOS1',sxpar(h,'FILTPOS1')
sxaddpar,header,'FILTPOS2',sxpar(h,'FILTPOS2')
	ENDIF
hist=[hist,'TARGET '+cub(k)]
cube_obj(*,*,e,k)=SHIFT(tmp-def_bias(tmp),0,20)
endfor
endfor

writefits,dir+'obj_i.fts',cube_obj,header
cont:
N_y=131
obj=readfits(dir+'obj_i.fts',h)
Nx=sxpar(h,'NAXIS1')
Ny=sxpar(h,'NAXIS2')
Nexp=sxpar(h,'NAXIS3')
Ncube=sxpar(h,'NAXIS4')

Npol=4
cube_obj=fltarr(Nx,N_y,Npol,Nexp,Ncube)
for j=0,Npol-1 do cube_obj(*,*,j,*,*)=obj(*,N_y*j:N_y*(j+1)-1,*,*)
K=0
tmp=fltarr(Nx,N_y)
for j=0,Nexp-1 do begin
window,2,xsize=Nx/2,ysize=N_y*Npol,title='exposure'+string(j+1)
for i=0,Npol-1 do begin
tmp(*,*)=cube_obj(*,*,i,j,K)
tv,255-bytscl(congrid(tmp,Nx/2,N_y),0,500),0,N_y*i
endfor
endfor
sxaddpar,h,'NAXIS',5
sxaddpar,h,'NAXIS2',N_y
sxaddpar,h,'NAXIS3',Npol,' Number of position polarization plane'
sxaddpar,h,'NAXIS4',Nexp,' Number of exposure object'
sxaddpar,h,'NAXIS5',Ncube, ' Number of data cube', after='NAXIS4'
writefits,dir+'cube_obj.fts',cube_obj,h




;
end