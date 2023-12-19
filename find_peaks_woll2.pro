
function find_peaks_WOLL2,vector,W=w,TRESH=tresh,PLOT=plot,MAX_val=max_val
;w - width of peaks
if not(keyword_set(W)) then W=10
if not(keyword_set(tresh)) then tresh=8
Ny=N_elements(vector)
y=findgen(Ny)
;вычитание тренда
;vector=vector-LOWESS(y,vector,Ny/4,2)

;оценка уровня отрезания
robomean,vector,3,0.5,avg_y,rms_y
;fi_peak,y,vector,avg_y+rms_y*tresh,ipix,xpk,ypk,bkpk,ipk
if keyword_set(max_val) then fi_peak,y,vector,max(vector)/tresh,ipix,xpk,ypk,bkpk,ipk
;определение точной позиции пика

R=where(xpk gt w/2 and xpk lt Ny-1-w/2,ind)
if ind gt 1 then begin
xpk=xpk(R)
for j=0,ind-1 do begin
f=goodpoly(y(xpk(j)-w/2:xpk(j)+w/2),vector(xpk(j)-w/2:xpk(j)+w/2),2,3)
xpk(j)=-f(1)/f(2)/2
endfor
endif
if keyword_set(plot) then begin
window,2
plot,y,vector, xst=1
oplot,[0,Ny],[1,1]*max(vector)/tresh
oplot,xpk,ypk,psym=6
endif
return,xpk
end
