;Stocks analyz
function extraction,ima,tra,WY=wy,S=S

a=size(ima)
vector=fltarr(a(1))
V=fltarr(a(2))
for k=0,a(1)-1 do begin
 V(*)=ima(k,*) & Vy=congrid(V,a(2)*S)
 yc=tra(k)*S
;	print,yc-wy,yc+wy,n_elements(Vy)
 vector(k)=total(Vy(yc-wy:yc+wy))
endfor
return,vector

end


dir='/hdd/Glagol/2019/WOLL-2/20191216/HD5797/'
cube=readfits(dir+'obj-sky.fts',h)
biny=2.0
apert=[2.0,2.0,2.0];in arcsec [5.0,5.0,5.0]
a=size(cube)
Nx=a(1) & Ny=a(2) & Npol=a(3) & Ncub=a(5)
Nexp=[sxpar(h,'NUMEXP1'),sxpar(h,'NUMEXP2'),sxpar(h,'NUMEXP3')]
lam=indgen(Nx)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
lam_st=5550 & lam_fin=7800
names=[sxpar(h,'NAME1'),sxpar(h,'NAME2'),sxpar(h,'NAME3')]
z=0.0
;emlines=[4861,4341]
;objpos=[49,Ny/2,Ny/2];*Ny/2
objpos=[Ny/2,Ny/2,Ny/2]

; построение траекторий
x=findgen(Nx) & wx=50
y=findgen(Ny) & wy=25/biny  ;2x1:  50;25  ;100 - когда далеко от центра кадра
avg_spectra=fltarr(Nx,Npol,Ncub)

sp=fltarr(Nx,Npol,Nexp(0),Ncub)

;FOR C=0,0 do begin;Ncub-1 DO BEGIN
FOR C=0,Ncub-1 do begin;Ncub-1 DO BEGIN
 w=7 ;7
 Npos=(Nx-100)/wx-1
 ypos=fltarr(Npos,Npol,Nexp(c))
 xpos=findgen(Npos)*wx+wx/2+100
  for j=0,Nexp(c)-1 do begin
   for i=0,Npol-1 do begin
    for k=0,Npos-1 do begin
     Vy=total(cube(xpos(k)-wx/2:xpos(k)+wx/2-1,*,i,j,c),1)/wx
     ;определение положения максимума
   ;		max_value=max(Vy,Nmax)
     max_value=max(Vy(objpos(c)-wy:objpos(c)+wy),Nmax) & Nmax=Nmax+objpos(c)-wy
     f=goodpoly(y(Nmax-w:Nmax+w),Vy(Nmax-w:Nmax+w),2,2,Yfit)
     ypos(k,i,j)=-f(1)/f(2)/2
    endfor
   endfor
  endfor

 tra=fltarr(Nx,i,j)

  for j=0,Nexp(c)-1 do begin
  for i=0,Npol-1 do begin
   f=goodpoly(xpos,ypos(*,i,j),2,2)
   tra(*,i,j)=f(0)+f(1)*x+f(2)*x^2;+f(3)*x^3
  endfor
 endfor

 avg_tra=fltarr(Nx,Npol)

 window,3,xsize=900
 !P.multi=[0,1,4]
 for i=0,Npol-1 do begin
  yc=tra(Nx/2,i,Nexp(c)/2)
  plot,[0,Nx-1],[-20,20]+yc,/nodata,xst=1,yst=1,charsize=2
  ;plot,[0,Nx-1],[100,200],/nodata,xst=1,yst=1,charsize=2
   for j=0,Nexp(c)-1 do oplot,x,tra(*,i,j),color=3e5
   for j=0,Nexp(c)-1 do cgoplot,xpos,ypos(*,i,j), psym=16
  avg_tra=total(tra,3)/Nexp(c)
  oplot,x,avg_tra(*,i),thick=2
 endfor

 ;экстракция спектров
 sc=1.05
 V=fltarr(Ny)
 spectra=fltarr(Nx,Npol,Nexp(c))

 ;optimal aperture
 slice=fltarr(Ny)
 for ky=0,Ny-1 do begin
  robomean,cube(Nx/2-5:Nx/2+5,ky,0,0,c),1,0.5,av
  slice(ky)=av
 endfor
 cgdisplay,wid=12
 !p.multi=[0,1,1]
 cgplot,slice
  ;lorentz
  lp=40/biny
  res=mpfitpeak(indgen(Ny),slice,a,/lorentz) & wait, 2
;		res=mpfitpeak(indgen(2*lp+1),slice(Ny/2-lp:Ny/2+lp-1),a,/lorentz) ;& wait, 2
  cgoplot,indgen(2*lp+1)+Ny/2-lp,res, color='grn6'
  stop
  cgoplot,[a(1)-a(2),a(1)+a(2),a(1)],[a(0)/2.0+a(3),a(0)/2.0+a(3),a(0)+a(3)]
  stop
  cgoplot,[a(1)-a(2),a(1)+a(2)],[a(3),a(3)]
  stop
   cgoplot,[a(1)-apert(c)/(0.1785*biny),a(1)-apert(c)/(0.1785*biny)],[-1e4,1e4],thick=2,color='gold'
   stop
   cgoplot,[a(1)+apert(c)/(0.1785*biny),a(1)+apert(c)/(0.1785*biny)],[-1e4,1e4],thick=2,color='gold'
   stop
  print, 'Optimal aperture: ', 9*a(2)
  wait,1
  stop

 for j=0,Nexp(c)-1 do begin
  for i=0,Npol-1 do begin
   ww=apert(c)/(0.1785*biny) & dt=0
   print, 'ww', ww, 'sc', sc
   spectra(*,i,j)=extraction(cube(*,*,i,j,(c)),tra(*,i,j)+dt,wy=ww,S=sc) ;5*a(2)
;				if c eq 0 then begin
;					ww=8.6/0.1785 & dt=14.1/0.1785
;					spectra(*,i,j)=extraction(cube(*,*,i,j,(c)),tra(*,i,j)+dt,wy=ww,S=sc)+extraction(cube(*,*,i,j,(c)),tra(*,i,j)-dt,wy=ww,S=sc)
;				endif
;				if c eq 1 then begin
;					ww=4.0/(0.1785*biny) & dt=0.0/(0.1785*biny)
;					spectra(*,i,j)=extraction(cube(*,*,i,j,(c)),tra(*,i,j)+dt,wy=ww,S=sc)
;				endif
  endfor
 endfor
 stop

 sp(*,*,0:Nexp(c)-1,c)=spectra(*,*,*)

 for j=0,Npol-1 do begin
  for kx=0,Nx-1 do begin
   avg_spectra(kx,j,c)=median(spectra(kx,j,*))
  endfor
 endfor

 for j=0,Npol-1 do begin
  for kx=30,Nx-5 do begin
   robomean,spectra(kx,j,*),1,0.5,avs
;			avs=total(spectra(kx,j,*))/N_elements(spectra(kx,j,*))
   avs=median(spectra(kx,j,*))
   avg_spectra(kx,j,c)=avs
  endfor
 endfor

 window,0,xsize=800,ysize=800
 !P.multi=[0,1,1]
 for j=0,Npol-1 do begin
  plot,avg_spectra(*,j,c),xst=1,/ylog,yrange=[5e3,5e4],/noerase

   for k=0,Nexp(c)-1 do oplot,spectra(*,j,k),color=3e5
  oplot,avg_spectra(*,j,c),thick=2
 endfor
stop
 ;деполяризация
 x0=700 & M=128
 avg_ratio=fltarr(2,Nexp(c))
 rms_ratio=fltarr(2,Nexp(c))

 window,2
 !P.multi=[0,1,2]

 num=findgen(Nexp(c))+1

 for k=0,Nexp(c)-1 do begin
  robomean,spectra(x0:x0+M-1,0,k)/spectra(x0:x0+M-1,1,k),2,0.5,mean,rms
  avg_ratio(0,k)=mean & rms_ratio(0,k)=rms
  robomean,spectra(x0:x0+M-1,2,k)/spectra(x0:x0+M-1,3,k),2,0.5,mean,rms
  avg_ratio(1,k)=mean & rms_ratio(1,k)=rms
 endfor

 for j=0,1 do begin
  robomean,avg_ratio(j,*),3,0.5,mean ;& print, mean
  avg_ratio(j,*)=avg_ratio(j,*)/mean
  rms_ratio(j,*)=rms_ratio(j,*)/mean
  plot,num,avg_ratio(j,*) ,psym=6,yrange=[0.95,1.05],yst=1,xrange=[0,Nexp(c)+1],xst=1
  oploterr,num,avg_ratio(j,*) ,rms_ratio(j,*),psym=6
  ;print, c, avg_ratio(j,*)
 endfor
ENDFOR

writefits,dir+'spectra.fit',sp,h
writefits,dir+'avg_spectra.fit',avg_spectra,h
;writefits,dir+'spectra_host.fit',sp,h
;writefits,dir+'avg_spectra_host.fit',avg_spectra,h


;end

Q=fltarr(Nx,3) & U=Q & F=Q

for c=0,2 do begin
 F(*,c)=total(avg_spectra(*,*,c),2)
 dy=0;0.2
 Q(*,c)=(avg_spectra(*,0,c)-shift_s(avg_spectra(*,1,c),dy))/(avg_spectra(*,0,c)+shift_s(avg_spectra(*,1,c),dy))*100
 dy=0
 U(*,c)=(avg_spectra(*,2,c)-shift_s(avg_spectra(*,3,c),dy))/(avg_spectra(*,2,c)+shift_s(avg_spectra(*,3,c),dy))*100
endfor
;U=shift(U,2,0)
;Q_0=LOWESS(findgen(Nx),Q(*,1),Nx/2,2,2)
;U_0=LOWESS(findgen(Nx),U(*,1),Nx/8,2,2)
P=fltarr(Nx,3)
FI=P
for c=0,2 do begin
 Q(*,c)=Q(*,c);-Q_0
 U(*,c)=U(*,c);-U_0
 P(*,c)=sqrt(Q(*,c)^2+U(*,c)^2)
 FI(*,c)=calc_atan(Q(*,c),u(*,c))
 R=where(FI(*,c) gt 180,ind) & if ind gt 1 then FI(R,c)=FI(R,c)-360
endfor

start_ps, dir+'stokes_ext.eps'
for ob=0,2 do begin
print, ob

 cgdisplay, wid=1, xsize=1200, ysize=1500
 !p.multi=[0,1,5]
 !p.charsize=1.8
 cgplot, lam, F(*,ob),xrange=[lam_st,lam_fin], ytitle='ADU, counts', $
  title=names(ob)+' 12/10/2020   BTA+SCORPIO-2'
  cgoplot, [1,1]*6563*(1+z), [-1e5,1.2*max(F(*,ob))], color='red'
  cgoplot, [1,1]*4861*(1+z), [-1e5,1.2*max(F(*,ob))], color='red'
 cgplot, lam, Q(*,ob), yrange=[median(Q(*,ob))-5,median(Q(*,ob))+5], xrange=[lam_st,lam_fin], ytitle='Q, %', psym=10
  cgoplot, [1,1]*6563*(1+z), [-1e5,1e5], color='red'
  cgoplot, [1,1]*4861*(1+z), [-1e5,1e5], color='red'
 cgplot, lam, U(*,ob), yrange=[median(U(*,ob))-5,median(U(*,ob))+5], xrange=[lam_st,lam_fin], ytitle='U, %', psym=10
  cgoplot, [1,1]*6563*(1+z), [-1e5,1e5], color='red'
  cgoplot, [1,1]*4861*(1+z), [-1e5,1e5], color='red'
 cgplot, lam, P(*,ob), yrange=[median(P(*,ob))-5,median(P(*,ob))+5], xrange=[lam_st,lam_fin], ytitle='P, %', psym=10
  cgoplot, [1,1]*6563*(1+z), [-1e5,1e5], color='red'
  cgoplot, [1,1]*4861*(1+z), [-1e5,1e5], color='red'
 cgplot, lam, FI(*,ob), xrange=[lam_st,lam_fin], charsize=1.8, ytitle='phi, deg', psym=10, $
   yrange=[median(FI(*,ob))-45,median(FI(*,ob))+45]
  cgoplot, [1,1]*4861*(1+z), [-1e5,1e5], color='red'
  cgoplot, [1,1]*6563*(1+z), [-1e5,1e5], color='red'

 if ob ne 2 then  ERASE

endfor

stop_ps

end