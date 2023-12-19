function remove_trend,vector
N=N_elements(vector)
f=goodpoly(findgen(N),vector,1,3,Yfit)
vector=(vector-Yfit)+Yfit(N/2)
return,vector
end
;********************************************
pro reduction_plot,x_in,y_in,y_out
dx=x_in(1)-x_in(0)

start:
Window,0,xsize=1900
plot,x_in,y_in,xst=1,yst=1,position=[0,0,1,1] ,/noerase,/norm,psym=10,color=3e5

;manual identefication spectra
sw=1
x_rep=0 & y_rep=0 & W_out=0
xpre=0 & ypre=0
be:
x=0 & y=0
if sw eq 4 then goto,final
cursor,x,y,1,/data
wait,0.1

a=!MOUSE
sw=LONG(a.button)
if xpre ne x and ypre ne y then begin
;
if sw eq 1 then begin
;NN=N_elements(x_out)
;if NN lt 3 then oplot,[x],[y],psym=6,color=2.2e8,symsize=2
;if NN gt 2 then

oplot,[x],[y],psym=7,color=2.2e8,symsize=2

x_rep=[x_rep,x]
y_rep=[y_rep,y]

endif
xpre=x & ypre=y
endif
if sw eq 1 then goto,be

final:
N=N_elements(x_rep)
for k=1,N-1 do begin
R=where(ABS(x_in-x_rep(k)) lt dx/2, ind)

if ind gt 0 then y_in(R)=y_rep(k)
endfor
oplot,x_in,y_in ,thick=1,psym=10
;
input=''
read,input,prompt='REPEAT? (y/n)'
if input eq 'y' then goto,start
y_out=y_in
END
;obj='MCG+08-11-011'
obj='Mkn110'
obj='Mkn231'
obj='3C273'
obj='PG0844+349'
obj='NGC4051'
file='h:\red_data.pol\Sy1\result\'+obj+'.txt'
set_plot,'WIN'


N=numlines(file) & print,N
res=fltarr(6,N)
openr,1,file
readf,1,res
close,1
;window,0
;plot,res(4,*),xst=1
;f=goodpoly(findgen(N),res(4,*),1,3,Yfit)
;oplot,Yfit
;oplot,res(4,*)/Yfit,color=3E5
;res(4,*)=res(4,*)/Yfit
;end
;reduction_plot,res(0,*),res(1,*),y_out &res(1,*)=y_out
;reduction_plot,res(0,*),res(5,*),y_out &res(5,*)=y_out
;res(5,*)=shift(res(5,*),-17)
;res(5,*)=remove_trend(res(5,*))
;******************************************************
;reduction_plot,res(0,*),res(4,*),y_out &res(4,*)=y_out
;******************************************************

dFI=0
res(5,*)=res(5,*)+dFI
res(2,*)=RES(4,*)*COS(res(5,*)*!PI/90)
res(3,*)=RES(4,*)*SIN(res(5,*)*!PI/90)

openw,1,'h:\red_data.pol\Sy1\result\'+obj+'.txt'
for k=0,N-1 do printf,1,res(*,k),format='(I4,E12.3,3F7.2,F8.1)'
close,1
END
