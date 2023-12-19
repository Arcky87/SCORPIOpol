;reduction_geometry_WOLL2
;формирование среднего спектра
wdir='h:\red_data.pol\3C390_140324\'
goto,cont1
eta=readfits(wdir+'eta_i.fts',h)
a=size(eta)
Ntra=3
N_deg=3
tra=fltarr(a(1),Ntra,a(3))
for k=0,a(3)-1 do begin
tra(*,*,k) =create_traectory_WOLL2(eta(*,*,k),NP=Ntra,WX=20,NDEG=N_deg,X_beg=250,/plot)
wait,1
endfor
sxaddpar,h,'OBJECT','TRAECTORY'
sxaddpar,h,'IMAGETYPE','eta'
sxaddpar,h,'Ndeg',FIX(N_deg),' Degree of polinomial approximation'
writefits,wdir+'tra.fit',tra,h
cont1:
neon=readfits(wdir+'neon_i.fts',h)
 tra=readfits(wdir+'tra.fit')
b=size(tra) & Ntra=b(2) & Nray=b(3)
a=size(neon) & Nx=a(1)  & Ny=a(2)
x=findgen(Nx)
map=fltarr(Nx,Ny)
;window,2,xsize=1000,ysize=140
;k=3

;map(*,*)=neon(*,*,k)& map=map+100 & map=ALOG10(map)
;tv,255-bytscl(congrid(map,1000,140),2,2.2)
;plot,[0,Nx-1],[0,Ny-1],/nodata,/noerase,position=[0,0,1,1],/norm,xst=1,yst=1
;for j=0,Ntra-1 do oplot,x,tra(*,j,k),color=1e7

;

geometry_2D,neon(*,*,k),tra(*,*,k),scale=2,X0,Y0,X1,Y1,/plot

end