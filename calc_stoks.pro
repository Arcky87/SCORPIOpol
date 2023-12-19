function calc_stoks,spectra,PLOT=plot,CORR=corr,Null=null

a=size(spectra) & N=a(1)
if not(keyword_set(null)) then null=fltarr(N,2)
if not(keyword_set(corr)) then corr=fltarr(N,2)+1
;print,size(corr)
kQ=corr(*,0) & kU=corr(*,1)
Stoks=fltarr(N,3)
Stoks(*,0)=spectra(*,0)+spectra(*,1)*kQ+spectra(*,2)+spectra(*,3)*kU
Stoks(*,1)=(spectra(*,0)-spectra(*,1)*kQ)/(spectra(*,0)+spectra(*,1)*kQ)-null(*,0)
Stoks(*,2)=(spectra(*,2)-spectra(*,3)*kU)/(spectra(*,2)+spectra(*,3)*kU)-null(*,1)
if keyword_set(plot) then begin
window,5,xsize=380,ysize=600
!P.multi=[0,1,3]
;
plot,stoks(*,0),xst=1
plot,stoks(*,1),xst=1;,position=[0.05,0.53,0.99,0.76],/norm,$
	;/noerase,charsize=1,xcharsize=1e-5
;xyouts,N/2,min(Q*100),'Q'

plot,stoks(*,2),xst=1;,position=[0.05,0.07,0.99,0.30],/norm,$
	;/noerase,charsize=1

endif
return,stoks
!P.multi=[0,1,1]
end