; traectory WOLL2

pro traectory_WOLL2,WDIR,WX=wx,NDEG=Ndeg,Xbeg=Xbeg,plot=plot

if not(keyword_set(Ndeg)) then Ndeg=3
;wdir=def_wdir(LOGFILE)
eta=readfits(wdir+'eta_i.fts',h)
a=size(eta)
Ntra=3
tra=findgen(a(1),Ntra,a(3))
Npos=FIX(a(1)/wx)-1;-3

dNpos=Xbeg/wx
Window,0,xsize=800,ysize=400
!P.multi=[0,1,1]
FOR J=0,a(3)-1 DO BEGIN
xpos=findgen(Npos)*wx+wx/2
ypos=fltarr(Npos,3)
;анализ центрального разреза
xcrsz=[1,1e-5,1e-5,1e-5]
y=findgen(a(2))
Vy=findgen(Npos,a(2))
tmp=findgen(3)
;формирование  узлов траектории
for k=0,Npos-1 do Vy(k,*)=total(eta(xpos(k)-wx/2:xpos(k)+wx/2,*,j),1)/(wx+1)
;справа от центра
ypos(Npos/2,*)=find_peaks_WOLL2(Vy(Npos/2,*),W=10,tresh=2,/max_val);,/plot
tmp(*)=ypos(Npos/2,*)
for k=1,Npos/2-1 do begin
for i=0,2 do tmp(i)=peak_position(Vy(Npos/2+k,*),tmp(i),5)
ypos(Npos/2+k,*)=tmp
endfor
;слева от центра
tmp(*)=ypos(Npos/2,*)
for k=1,Npos/2-dNpos do begin
for i=0,2 do tmp(i)=peak_position(Vy(Npos/2-k,*),tmp(i),5)
ypos(Npos/2-k,*)=tmp
endfor
plot,[0,a(1)-1],[1,a(2)-1],/nodata,xst=1,yst=1,position=[0.05,0.07+0.23*j,0.85,0.07+0.23*(j+1)],/norm,$
	noerase=j,xcharsize=xcrsz(j),ytickinterval=50
x=findgen(a(1))
for k=0,2 do begin
oplot,xpos(dNpos:Npos-1),ypos(dNpos:Npos-1,k),psym=6,symsize=0.5
f=goodpoly(xpos(dNpos:Npos-1),ypos(dNpos:Npos-1,k),Ndeg,1,fit)
robomean,ypos(dNpos:Npos-1,k)-fit,2,0.5,avg_fit,rms_fit
Yfit=0
for i=0,Ndeg do Yfit=Yfit+f(i)*x^i
tra(*,k,j)=Yfit
oplot,x,Yfit
xyouts,a(1),Yfit(a(1)-1),' rms='+string(rms_fit,format='(F5.2)')+' px'
endfor
ENDFOR
writefits,wdir+'tra.fit',tra
end
LOGFILE=DIALOG_PICKFILE(/read,path='h:\red_data.pol\AGN\LOGS\',FILTER='*.txt')

traectory_WOLL2,LOGFILE,WX=40,NDEG=3,Xbeg=100
END
