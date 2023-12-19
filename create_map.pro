;CREATE_MAP
DIR='h:\red_data.pol\AGN\NGC5548_160307\maps\' & m1=1 & m2=4
DIR='h:\red_data.pol\AGN\Mkn744_160307\maps\' & m1=1 & m2=4
DIR='h:\red_data.pol\AGN\MCG+08-11-011_151106\maps\' & m1=0.5 & m2=3
DIR='h:\red_data.pol\AGN\Mkn79_151207\maps\' & m1=0.5 & m2=3.5
 ;DIR='h:\red_data.pol\AGN\Mkn110_151207\maps\' & m1=1 & m2=5
;DIR='h:\red_data.pol\AGN\NGC3227_151210\maps\' & m1=0.5 & m2=3
;DIR='h:\red_data.pol\AGN\NGC4151_160307\maps\' & m1=0.5 & m2=4
;DIR='h:\red_data.pol\AGN\NGC4051_160307\maps\' & m1=1 & m2=3.5
;DIR='h:\red_data.pol\AGN\NGC4593_160308\maps\' & m1=1 & m2=4
;DIR='h:\red_data.pol\AGN\NGC7469_151210\maps\' & m1=0.5 & m2=4.5
;;xc=335 & yc=968  & wc=10
;xc=515 & yc=518  & wc=10
;xg=516 & yg=518  & wc=10
files=FILE_BASENAME(FILE_SEARCH(dir+'s*.fts'))
print,files
h=headfits(dir+files(0))
Nx=sxpar(h,'NAXIS1')  & Ny=sxpar(h,'NAXIS2')
N=N_elements(files)
print,N
cube=fltarr(Nx,Ny,N)
for k=0,N-1 do cube(*,*,k)=FLOAT(readfits(dir+files(k),/silent))
; вычитание фона неба
for k=0,N-1 do begin
robomean,cube(*,*,k),2,0.5,avg_ima,rms_ima
bin=FIX(rms_ima) & M=10
xcross=findgen(2*M+1)*bin+avg_ima-M*bin
ycross=histogram(cube(*,*,k),min=avg_ima-M*bin,max=avg_ima+M*bin,binsize=bin)
gau=gaussfit(xcross,ycross,G,Nterms=3)
NS=G(1)
;
;определение смешения изображения

cube(*,*,k)=cube(*,*,k)-NS
cube(*,*,k)=cube(*,*,k)/total(cube(*,*,k))*1E7
endfor

ima=fltarr(Nx,Ny)
; убирание частиц
for x=0,Nx-1 do begin
for y=0,Ny-1 do ima(x,y)=median(cube(x,y,*))
endfor

ws=30
vector=ALOG10(total(ima(*,Ny/2-ws:Ny/2+ws),2))
x=findgen(Nx)

max_value=max(vector(Nx/2-ws:Nx/2+ws),Nmax)
xs=x(Nx/2-ws+Nmax)+1
ys=Ny/2+10
print
W=Nx/8
map=FLOAT(ima(xs-w:xs+w-1,ys-w:ys+w-1))
dw=2 & d_w=2*dw
set_plot,'PS'
device,file=dir+'map.ps',/portrait,xsize=16,ysize=24,xoffset=4,yoffset=2,bits_per_pixel=24
;window,2,xsize=Nx/2,ysize=Ny/2
tv,255-bytscl(ALOG10(ROTATE(map,3 )),m1,m2),0,8,xsize=16 ,/centimeter
plot,[0,2*w],[0,2*w],/noerase,xst=1,yst=1,position=[0,0.3333,1,1],/norm,$
	/nodata,xcharsize=1e-5,charsize=1
	oplot,[0,2*w],[1,1]*(w-dw)
	oplot,[0,2*w],[1,1]*(w+dw)
	oplot,[0,2*w],[1,1]*(w-dw-d_w),linestyle=2
	oplot,[0,2*w],[1,1]*(w+dw-d_w),linestyle=2
	oplot,[0,2*w],[1,1]*(w-dw+d_w),linestyle=2
	oplot,[0,2*w],[1,1]*(w+dw+d_w),linestyle=2
center=TOTAL(map(w-dw:w+dw,*),1)/(2*dw+1)
galaxy=(TOTAL(map(w-dw+d_w:w+dw+d_w,*),1)/(2*dw+1)+TOTAL(map(w-dw-d_w:w+dw-d_w,*),1)/(2*dw+1))/2
nucleus=center-galaxy
;fit = MPFITPEAK(findgen(Nx/2), nucleus, A,/MOFFAT)
fit = MPFITPEAK(findgen(Nx/2), nucleus, A,/Gauss)
scale=0.357
print,A(2)
;
plot,center,xst=1,yrange=[0.1,max(center)*1.2],yst=1,/ylog ,xrange=[0,2*w],$
	subtitle='galaxy='+string(total(galaxy),format='(I6)')+'  nucleus='+$
	string(total(nucleus),format='(I6)')+'  error='+$
	string(total((nucleus-fit)^2)/Nx,format='(I4)')+'  seeing='+$
	string(A(2)*scale,format='(f5.2)'),$
	position=[0,0,1,0.3333],/norm,/noerase,charsize=1
oplot,galaxy,color=150,thick=2
oplot,center
oplot,nucleus
oplot,fit,thick=2,color=150
xyouts,0.5,1.01,dir,/norm,align=0.5,charsize=1.5
device,/close
set_plot,'WIN'
window,3,xsize=2*w+1,ysize=2*w+1
;tv,255-bytscl(ALOG10(ROTATE(map,3 )),m1,m2)
print,min(map),max(map)
tv,255-bytscl(map,0,10000)
sxaddpar,h,'BZERO',0
writefits,dir+'map.fts',map,h
;print,total(galaxy),total(nuclea),total((nuclea-fit)^2)/Nx
;oplot,nuclea-fit
end