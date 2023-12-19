;calibration_UV_Rc
name=['Mkn335','Akn120','Mkn79','Mkn110','NGC3227','NGC4151','NGC4051','3C273','NGC4593','IRAS13349','NGC5548','Mkn817','Mkn509','Mkn1502','Mkn231','Mkn6']
  dx=[      1,      1.5,      1,       1,        1,        1,      1.5,      1,        1,          1,        1,       1,       1,        1,    0.7,     1]
  UV=[ 9.5896,  58.8224, 0.6166, 20.3210,   0.2302,   0.8336,    0.1223,   1918,    1.748,   1105.78,   1.4394, 19.2434, 35.1134,  10.5919, 49.3674,   42.]
  UV=DOUBLE(UV)*10.^38
  UV=UV*10^5.
  Rc=[    142,     452,      55,    100,        25,       36,      44,     963,       43,        1130,      114,     180,      130,      90,     404,  214]
N=N_elements(name)
window,0
plot,UV,Rc,yst=1,/ylog ,/xlog ,yrange=[10,2E3], psym=6,title='log(Rc)=-14.866+0.386*log(UV)',$
	xtitle='UV luminosity, erg/sec',ytitle='Rc, ligth days'
for k=0,N-1 do xyouts,UV(k)*dx(k),Rc(k)*1.1,name(k),align=0.5
f=poly_fit(ALOG10(UV),ALOG10(Rc),1,sigma=err)
UV_lim=DOUBLE([1,1.0e5])*10.0^35*10.^7
print,UV_lim
R_lim=10.^f(0)*UV_lim^f(1)
print,R_lim
print,f
print,err
oplot,UV_lim,R_lim,linestyle=1
UV=[34.7511,7.0658,36.7829,20.7672,36.4078,4.1518,2.5204,0.883,19.2432,44.2022,53.7836,2782.85,43.7617,435.337,24.3]
UV_new=DOUBLE(UV)*10.^38
UV_new=UV_new*10^5.
N=N_elements(UV)
Rc=10.^f(0)*UV_new^f(1)
for k=0,N-1 do print,UV(k),Rc(k)
;oplot,UV_new,Rc,psym=1
 end