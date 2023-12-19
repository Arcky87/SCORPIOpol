function corr_2order,spectra,PLOT=plot,LINEAR=linear,WTITLE=wtitle,NUMWIN=numwin
if not(keyword_set(wtitle)) then wtitle=''
if not(keyword_set(numwin)) then numwin=0
bobo=3
Ndeg=2
Nbeg=linear(0)  & Nend=linear(1)
a=size(spectra) & N=a(1)
x=findgen(N) & out=fltarr(N,2)
corr=fltarr(N,2)
kQ=spectra(*,0)/spectra(*,1)
f=goodpoly(x(Nbeg:Nend),kQ(Nbeg:Nend),Ndeg,3)
kQlin=0 & for j=0,Ndeg do kQlin=kQlin+f(j)*x^j
kQfit=LOWESS(x,kQ,N/bobo,3,3)
Q=(spectra(*,0)-spectra(*,1)*kQfit/kQlin)/(spectra(*,0)+spectra(*,1)*kQfit/kQlin)

kU=spectra(*,2)/spectra(*,3)
f=goodpoly(x(Nbeg:Nend),kU(Nbeg:Nend),Ndeg,3)
kUlin=0 & for j=0,Ndeg do kUlin=kUlin+f(j)*x^j

kU=median(kU,5)
kUfit=LOWESS(x,kU,N/bobo,3,1)


U=(spectra(*,2)-spectra(*,3)*kUfit/kUlin)/(spectra(*,2)+spectra(*,3)*kUfit/kUlin)
if keyword_set(plot) then begin
window,20+numwin,title=wtitle,xsize=500,ysize=270,xpos=1400,ypos=310*numwin
!P.multi=[0,1,2]
plot,kQ/kQlin,xst=1,yrange=[0.6,1.4],yst=1,$
	position=[0.05,0.52,0.99,0.99],/norm,$
	charsize=1,xcharsize=1e-5
	xyouts,N/2,1.3,'I!D0!N(0)/I!D1!N(90)',align=0.5,charsize=1.5
oplot,kQfit/kQlin,color=3e5,thick=2
;oplot,kQfit,color=3e5,thick=2
;oplot,kQlin,thick=2,color=1.3e5
plot,kU/kUlin,xst=1,yst=1,$;yrange=[0.6,1.4]
	position=[0.05,0.05,0.99,0.52],/norm,$
	charsize=1,/noerase
	xyouts,N/2,0.65,'I!D2!N(45)/I!D3!N(135)',align=0.5,charsize=1.5

oplot,kUfit/kUlin,color=3e5,thick=2
;oplot,kUlin,thick=2,color=1.3e5
endif
corr(*,0)=kQfit/kQlin & corr(*,1)=kUfit/kUlin

return,corr
!P.multi=[0,1,1]
end
;dir='d:\Sy1\NGC4151_160307\'
dir='d:\Sy1\Mkn817_140529\'
avg_spectra=readfits(dir+'avg_spectra.fit',h,/silent)
a=size(avg_spectra) & Nx=a(1) & x=findgen(Nx)
R=where(avg_spectra eq 0 ,ind) & if ind gt 1 then avg_spectra(R)=1
slope=slope_flat(dir)  & print,'slope',slope
if slope eq 1 then avg_spectra=tmp(*,reverse(num),*)
;

j=1
kQ=avg_spectra(*,0,j)/avg_spectra(*,1,j)
kU=avg_spectra(*,2,j)/avg_spectra(*,3,j)
window,0
!P.multi=[0,1,2]
bobo=2
plot,kQ,xst=1,yrange=[0.,2],yst=1
oplot,LOWESS(x,kQ,N/bobo,3,1),color=3e5
plot,kU,xst=1,yrange=[0.,2],yst=1
oplot,LOWESS(x,kU,N/bobo,3,1),color=3e5
stoks=calc_stoks(avg_spectra(*,*,j));,/plot)
star_0=fltarr(N,2)
for k=0,1 do star_0(*,k)=lowess(x,stoks(*,k+1),N/2,2,2)
stoks=calc_stoks(avg_spectra(*,*,0),null=star_0,/plot);,/plot
window,1
!P.multi=[0,1,3]
plot,stoks(*,0),xst=1
plot,stoks(*,1),xst=1
plot,stoks(*,2),xst=1
end

;