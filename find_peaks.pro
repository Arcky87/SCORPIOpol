function find_peaks,vector,W=w,TRESH=tresh,PLOT=plot
;w - width of peaks
if not(keyword_set(W)) then W=4 ; IY default 10
if not(keyword_set(tresh)) then tresh=8 ; IY default 5
Ny=N_elements(vector)
y=findgen(Ny)

;trend subtraction
vector=vector;-LOWESS(y,vector,Ny/4,2)

;estimation of the tresh level
robomean,vector,3,0.5,avg_y,rms_y & print, avg_y
rms_y=45 ;!!!!
fi_peak,y,vector,avg_y+rms_y*tresh,ipix,xpk,ypk,bkpk,ipk
;print, '!!!', xpk

;find exact position of peak
R=where(xpk gt w/2 and xpk lt Ny-1-w/2,ind);(xpk gt 0 and xpk lt Ny-1-w/2,ind) ;(xpk gt w/2 and xpk lt Ny-1-w/2,ind)
if ind gt 1 then begin
 xpk=xpk(R)
  for j=0,ind-1 do begin
   f=goodpoly(y(xpk(j)-w/2:xpk(j)+w/2),vector(xpk(j)-w/2:xpk(j)+w/2),2,5)
   xpk(j)=-f(1)/f(2)/2
  endfor
endif

if keyword_set(plot) then begin
 cgdisplay,wid=2
 !p.multi=[0,1,1]
 cgplot,y,vector, xst=1
 cgoplot,[0,Ny],[1,1]*rms_y*tresh,color='red'
 cgoplot,xpk,ypk,psym=16
 wait, 0.5
endif

return,xpk
end
