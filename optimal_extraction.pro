;test ext
function  optimal_extraction,OBJ,DPOS=dpos,PLOT=plot,FWHM=fwhm,NDEG=Ndeg
;оптимальная экстракция спектра из изображения
;Ndeg=1
tresh=2
a=size(obj) & Nx=a(1)  & Ny=a(2)
y=findgen(Ny) & x=findgen(Nx)
;d_lambda=sxpar(h,'CDELT1')
;lambda=findgen(Nx)*d_lambda+sxpar(h,'CRVAL1')

Nobj=N_elements(dpos) ; print,Nobj
spectra=fltarr(Nx,Nobj)
w=fltarr(Nobj)+FWHM/2


print,'OPTIMAL EXTRACTION'
;create traectory
Npos=20 & wx=Nx/(Npos+1)
xpos=findgen(Npos)*wx+wx
if xpos(Npos-1)+wx gt Nx-1 then xpos(Npos-1)=Nx-wx-1
ypos=fltarr(Npos,Nobj) & FWHMpos=ypos
vector=fltarr(Ny)
for k=0,Npos-1 do begin
vector(*)=total(obj(xpos(k)-wx:xpos(k)+wx,*),1)
res=MULTIGAUS(y,vector,dpos,FWHM=w)
ypos(k,*)=res.center
FWHMpos(k,*)=res.FWHM
endfor
posfit=fltarr(Nx,Nobj)
FWHMfit=fltarr(Nx,Nobj)
;polynomial approximation
for k=0,Nobj-1 do begin
	f=goodpoly(xpos,ypos(*,k),Ndeg,tresh)
		for j=0,Ndeg do posfit(*,k)=posfit(*,k)+f(j)*x^j
	f=goodpoly(xpos,FWHMpos(*,k),Ndeg,tresh)
		for j=0,Ndeg do FWHMfit(*,k)=FWHMfit(*,k)+f(j)*x^j
endfor
;extraction spectra
eps=0
	for k=eps,Nx-1-eps do begin
	vector=total(obj(k-eps:k+eps,*),1)/(2*eps+1)
	res=MULTIGAUS(y,vector,posfit(k,*),FWHM=FWHMfit(k,*),FIXpos=0,FIXfwhm=1)
	spectra(k,*)=res.flux
		endfor
	if keyword_set(plot) then begin
window,11
!P.multi=[0,1,2]
plot,[0,Nx],[min(dpos)-10,max(dpos)+10],$
	xst=1,yst=1,/nodata,title='traectory'
for k=0,Nobj-1 do begin
oplot,xpos,ypos(*,k),psym=k+5
oplot,x,posfit(*,k)
endfor
plot,[0,Nx],[0,20],xst=1,yst=1,/nodata,title='variation FWHM'
for k=0,Nobj-1 do begin
oplot,xpos,FWHMpos(*,k),psym=k+5
oplot,x,FWHMfit(*,k)
endfor
	endif
cont:
!P.multi=[0,1,1]
return,spectra
end
!P.multi=[0,1,1]
LOGFILE='h:\spectraPOL.log\GRW+70_100813.txt'
wdir=sxpar(read_table(LOGFILE),'w_dir')
cube=readfits(wdir+'obj-sky.fts',h)
a=size(cube)

k=0
j=1
;spectra=optimal_extraction(cube(*,*,k,j),DPOS=82,FWHM=3,NDEG=2,/plot)
window,2
ima=cube(*,*,k,j)
y=findgen(a(2))
plot,y(50:150),ima(2200,50:150),xst=1
gau=gaussfit(y(50:150),ima(2200,50:150),G,Nterms=4)
oplot,y(50:150),gau-G(3)
end
