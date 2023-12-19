;Rc_VS_UV
dir='h:\red_data.pol\AGN\'
tab=read_table(dir+'Rc_VS_UV.txt')
N=N_elements(tab)
flux=fltarr(N) & Rc=flux & obj=strarr(N)
for k=0,N-1 do begin
tmp=str_sep(STRCOMPRESS(tab(k)),' ')
 obj(k)=tmp(0)
flux(k)=FLOAT(tmp(1))
  Rc(k)=FLOAT(tmp(3))
  print,obj(k), flux(k), Rc(k)
  endfor
  print,flux
  print,Rc
Rc=Rc(sort(Flux))
obj=obj(sort(Flux))
Flux=Flux(sort(flux))

Window,0
!P.multi=[0,1,1]

f=goodpoly(ALOG10(flux),ALOG10(Rc),1,3,fit)
print,f
plot,flux,Rc,psym=6,/ylog,/xlog
for k=0,N-1 do xyouts,flux(k),Rc(k)*1.1,obj(k),align=0.5
oplot,flux,10^(f(0)+f(1)*ALOG10(flux))
end