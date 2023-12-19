;test_param
dir='h:\red_data.pol\tinatin\2MASXJ060210+2829_141121\
N=numlines(dir+'test.txt')
tab=fltarr(2,N)
close,1
openr,1,dir+'test.txt'
readf,1,tab
close,1
Window,0
!P.multi=[0,1,2]
z=0.033  & wav=[4861,5500,6563]
xrng=[4700,7500]
R=where(tab(0,*) gt xrng(0) AND tab(0,*) lt xrng(1))
plot,ALOG10(tab(0,R)),ALOG10(tab(1,R)),xst=1
for k=0,2 do oplot,[1,1]*ALOG10(wav(k)*(1+z)),[0,2],linestyle=2
;f=goodpoly(ALOG10(tab(0,R)),ALOG(tab(1,R)),1,3,FIT)
f=poly_fit(ALOG10(tab(0,R)),ALOG10(tab(1,R)),1,FIT,SIGMA=rms)
print,f(1),rms(1)

;print,f(1),rms(1)
;window,1


plot,tab(0,R),tab(1,R),xst=1
oplot,tab(0,R),10^FIT
for k=0,2 do oplot,[1,1]*wav(k)*(1+z),[0,10],linestyle=2
d_wave=20
for k=0,2 do begin

RW=where(tab(0,*) gt wav(k)*(1+z)-d_wave AND  tab(1,*) lt wav(k)*(1+z)+d_wave)
robomean,tab(1,RW),2,0.5,avg_val,rms_val
robomean,10^FIT(RW)+.1,2,0.5,avg_fit,rms_fit
print,wav(k),avg_val,rms_val,avg_fit,rms_fit
endfor
end