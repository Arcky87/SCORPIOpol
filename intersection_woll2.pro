function intersection_WOLL2,lines,traectory,W,PLOT=plot
ext=70  ;70 - 2x1
Nx=N_elements(traectory)
Ny=N_elements(lines)
x=findgen(Nx)-ext  & y=findgen(Ny)-ext
ROBOMEAN,TRAECTORY,3,0.5,YC
xc=total(lines(yc-w+ext:yc+w+ext))/(2*w+1)

T=goodpoly(x(xc-w+ext:xc+w+ext),traectory(xc-w+ext:xc+w+ext),1,2,fit)
L=goodpoly(lines(yc-w+ext:yc+w+ext),y(yc-w+ext:yc+w+ext),1,1,fit)
;print,'Here at L (intersection)'
xc=-(L(0)-T(0))/(L(1)-T(1))
yc=T(0)+T(1)*xc
if keyword_set(plot) then begin
window,0
plot,lines,y,xst=1,yst=1
oplot,x,traectory
oplot,[xc],[yc],color=10,psym=6
wait,0.5
endif
return,[xc,yc]
end
