function intersection,lines,traectory,W,PLOT=plot
ext=20
Nx=N_elements(traectory)
Ny=N_elements(lines)
x=findgen(Nx)-ext  & y=findgen(Ny)-ext
yc=total(traectory)/Nx
xc=total(lines(yc-w+ext:yc+w+ext))/(2*w+1)
T=goodpoly(x(xc-w+ext:xc+w+ext),traectory(xc-w+ext:xc+w+ext),1,2,fit)
L=goodpoly(lines(yc-w+ext:yc+w+ext),y(yc-w+ext:yc+w+ext),1,2,fit)
xc=-(L(0)-T(0))/(L(1)-T(1))
yc=T(0)+T(1)*xc
if keyword_set(plot) then begin
window,0
plot,lines,y,xst=1,yst=1
oplot,x,traectory
oplot,[xc],[yc],psym=6
print,yc,xc
endif
return,[xc,yc]
end