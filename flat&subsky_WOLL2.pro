;create_flat_WOLL-2
;dir='h:\red_data.pol\3C390_140324\'
!P.multi=[0,1,1]
;LOGFILE=DIALOG_PICKFILE(/read,path='h:\red_data.pol\AGN\LOGS\',FILTER='*.txt')
LOGFILE=DIALOG_PICKFILE(/read,path='h:\red_data.pol\LOGS\',FILTER='*.txt')
print,systime()
dir=def_wdir(LOGFILE)
flat=readfits(dir+'flat_lin.fts',h)
a=size(flat) & Nx=a(1)  & Ny=a(2)  & Npol=a(3)
;goto,cont
map=fltarr(Nx,Ny)
norm=fltarr(Nx,Npol)

cross=fltarr(Ny,Npol)
w=5
;формирование нормировки плоского поля
;Window,2,xsize=Nx/2,ysize=Ny*4
for k=0,Npol-1 do begin
map(*,*)=flat(*,*,k)
for ky=0,Ny-1 do map(*,ky)=smooth(map(*,ky),50,/edge_truncate)
norm(*,k)=total(map(*,Ny/2-w:Ny/2),2)/(2*w+1)
cross(*,k)=total(map(Nx/2-w:Nx/2,*),1)/(2*w+1)
;tv,255-bytscl(congrid(map,Nx/2,Ny),0,3e4),0,Ny*k
endfor


x=findgen(Nx)
for j=0,Npol-1 do norm(*,j)=LOWESS(x,norm(*,j),Nx/32,2,2)
norm=total(norm,2)/4
for k=0,Npol-1 do begin
for ky=0,Ny-1 do begin
flat(*,ky,k)=flat(*,ky,k)/norm
flat(*,ky,k)=LOWESS(x,flat(*,ky,k),Nx/16,2,2)
print,k,ky
wait,0.02
end
endfor
;Window,3,xsize=Nx/2,ysize=Ny*4
;for k=0,Npol-1 do begin
;map(*,*)=flat(*,*,k)
;tv,255-bytscl(congrid(map,Nx/2,Ny)),0,Ny*k
;endfor
writefits,dir+'flat_norm.fts',flat
;cont:
flat=readfits(dir+'flat_norm.fts')
w=3
norm=fltarr(Npol)
for k=0,Npol-1  do begin
norm=total(flat(Nx/2-w:Nx/2+w,Ny/2-w:Ny/2+w,k))/(2*w+1)/(2*w+1)
flat(*,*,k)=flat(*,*,k)/norm
endfor
;Window,2,xsize=Nx/2,ysize=Ny*4
;for k=0,Npol-1 do begin
;map(*,*)=flat(*,*,k)
;tv,255-bytscl(congrid(map,Nx/2,Ny),0.5,1.5),0,Ny*k
;endfor
writefits,dir+'avg_flat.fts',flat,h
;window,3
;plot,[0,Ny],[0.5,1.5],xst=1,yst=1,/nodata
;y=findgen(Ny)
;w=10
;for j=0,3 do oplot,y,total(flat(Nx/2-w:Nx/2+w,*,j),1)/(2*w+1)

subsky_WOLL2,LOGFILE,/plot
print,systime()
ViewPol_2
end
