;check_stoks_WOLL2
;проверка спектров перед образование векторов Сещкса
;goto,fin
log_dir='h:\red_data.pol\Sy1\'
;log_dir='h:\red_data.pol\AGN\'
LOGFILE=DIALOG_PICKFILE(/read,path=log_dir+'LOGS\',FILTER='*.txt')
;LOGFILE=log_dir+'LOGS\'+'Mkn876_160408.txt
	name_out=str_sep(FILE_BASENAME(LOGFILE),'.')
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
	wdir=log_dir+wdir(N_elements(wdir)-2)+'\'
;проверка шкалы длин волн
neon=readfits(wdir+'neon_lin.fts',h)

Nx=sxpar(h,'NAXIS1')
Ny=sxpar(h,'NAXIS2')
yc=Ny/2  & wy=20
vector=total(neon(*,yc-wy:yc+wy,*),2)
x_shift=fltarr(4)
M=10 & xcross=findgen(2*M+1)-M
for k=0,3 do begin
ycross=cross_norm(vector(*,k),vector(*,0),M)
gau=gaussfit(xcross,ycross,G)
x_shift(k)=G(1)
endfor
print,x_shift,format='(4F7.2)'

;vector(*,k)=shift_s(vector(*,k),-x_shift(k))


;формирование среднего спектра
;

;s=readfits(wdir+'avg_spectra.fit',h)

spectra=readfits(wdir+'spectra.fit',hs)
Nx=sxpar(hs,'NAXIS1') & Npol=sxpar(hs,'NAXIS2') & Ntarget=sxpar(hs,'NAXIS4')
Nexp=intarr(Ntarget) & wx=50  & Name=strarr(Ntarget) & PA=fltarr(Ntarget)
for k=0,Ntarget-1 do begin
Nexp(k)=sxpar(hs,'NUMEXP'+string(k+1,format='(I1)'))
Name(k)=sxpar(hs,'NAME'+string(k+1,format='(I1)'))
PA(k)=sxpar(hs,'PA'+string(k+1,format='(I1)'))
endfor
star_name=name(2)
PA_slit=PA(2)
mean_spectra=fltarr(Nx,Npol,Ntarget)
;проверка совместимости спектров
norm=fltarr(Npol,max(Nexp),Ntarget)
for k=0,Ntarget-1 do begin
for j=0,Npol-1 do begin
for i=0,Nexp(k)-1 do begin
robomean,spectra(Nx/2-wx:Nx/2+wx,j,i,k),3,0.5,avg_spectra,rms_spectra
norm(j,i,k)=avg_spectra
spectra(*,j,i,k)=spectra(*,j,i,k)/norm(j,i,k)
spectra(*,j,i,k)=shift_s(spectra(*,j,i,k),-x_shift(j))
endfor
spectra(*,j,*,k)=spectra(*,j,*,k)*max(norm(j,*,k))
for x=0,Nx-1 do mean_spectra(x,j,k)=MEDIAN(spectra(x,j,0:Nexp(k)-1,k))
endfor
endfor

title_target=['object','unpolarized star','polarized star']
angle=['0','90','45','135']+' deg'
set_plot,'PS'
device,file=wdir+'spectra.ps',/landscape
;window,4,xsize=1200,ysize=600
!P.multi=[0,1,1]
for k=0,Ntarget-1 do begin
y_max=max(mean_spectra(*,*,k))
plot,[0,1],[0,1],/nodata,position=[0.005+0.33*k,0.01,0.005+0.33*(k+1),0.97],/norm, noerase=k,$
   TICKLEN=0,xcharsize=1e-5,ycharsize=1e-5,$
   title=title_target(k)+'  '+name(k)+' '+'Nexp='+string(Nexp(k),format='(I2)')+' PA='+STRING(PA(k),format='(F7.1)'),charsize=0.6
for j=0,Npol-1 do begin
plot,[0,Nx-1],[0,y_max],xst=1,yst=1,$
	position=[0.005+0.33*k,0.01+0.24*j,0.005+0.33*(k+1),0.01+0.24*(j+1)],$
	/norm,/noerase,/nodata,TICKLEN=0,charsize=1e-5
for i=0,Nexp(k)-1 do oplot,spectra(*,j,i,k),color=120
oplot,mean_spectra(*,j,k)
xyouts,Nx*0.85,y_max*0.85,angle(j),charsize=0.6
endfor
endfor
xyouts,0,-0.01,wdir,/norm,charsize=0.6
device,/close
set_plot,'WIN'
;калибровка нулевого стандакта
s=mean_spectra
window,1,xsize=600,ysize=427
!P.multi=[0,1,2]
x=findgen(Nx)
dx=100
beg=dx
cut=Nx-dx
;& deg=5
R01=s(*,0,1)/s(*,1,1)
R23=s(*,2,1)/s(*,3,1)
plot,R01,xst=1,charsize=1,color=1e5,yrange=[0,2]
fit=LOWESS(x(beg:cut),R01(beg:cut),Nx/2 ,3,1)
R01=INTERPOL( fit,x(beg:cut),x)
oplot,x,R01
plot,R23,xst=1,charsize=1,color=1e5,yrange=[0,2]
fit=LOWESS(x(beg:cut),R23(beg:cut),Nx/2 ,3,1)
R23=INTERPOL( fit,x(beg:cut),x)
oplot,x,R23
;END





name=strarr(3) & PA=fltarr(3)
for j=0,2 do begin
name(j)=sxpar(h,'NAME'+string(j+1,format='(I1)'))


endfor

window,1,xsize=600,ysize=427
!P.multi=[0,1,2]
x=findgen(Nx)

;cut=4200
;& deg=5
R01=s(*,0,1)/s(*,1,1)
R23=s(*,2,1)/s(*,3,1)
plot,R01,xst=1,charsize=1,color=1e5,yrange=[0,3]
fit=LOWESS(x(beg:cut),R01(beg:cut),Nx/2 ,3,1)
R01=INTERPOL( fit,x(beg:cut),x)
oplot,x,R01
plot,R23,xst=1,charsize=1,color=1e5,yrange=[0,3]
fit=LOWESS(x(beg:cut),R23(beg:cut),Nx/2 ,3,1)
R23=INTERPOL( fit,x(beg:cut),x)
oplot,x,R23

;
stoks=fltarr(Nx,4,3)
cube_stoks=fltarr(Nx,4,Nexp(0),3)
;Nx=2291
;s=s(200:Nx-200,*,*)
a=size(s) & Nx=a(1) & x=findgen(Nx)
wave=findgen(Nx)*sxpar(h,'CDELT1') +sxpar(h,'CRVAL1')
Q=fltarr(Nx,3) & U=Q
;***********************
mode=slope_flat(wdir)
;**********************
;R01=1 & R23=1
for k=0,Ntarget-1 do begin
stoks(*,0,k)=s(*,0,k)+R01*s(*,1,k)+s(*,2,k)+R23*s(*,3,k)
stoks(*,1,k)=(s(*,0,k)-R01*s(*,1,k))/(s(*,0,k)+R01*s(*,1,k))
stoks(*,2,k)=(s(*,2,k)-R23*s(*,3,k))/(s(*,2,k)+R23*s(*,3,k))
for j=0,Nexp(k)-1 do begin
cube_stoks(*,0,j,k)= spectra(*,0,j,k)+R01*spectra(*,1,j,k)+spectra(*,2,j,k)+R23*spectra(*,3,j,k)
cube_stoks(*,1,j,k)=(spectra(*,0,j,k)-R01*spectra(*,1,j,k))/(spectra(*,0,j,k)+R01*spectra(*,1,j,k))
cube_stoks(*,2,j,k)=(spectra(*,2,j,k)-R23*spectra(*,3,j,k))/(spectra(*,2,j,k)+R23*spectra(*,3,j,k))
endfor
endfor
Q_null=Smooth(stoks(*,1,1),wx,/edge_truncate)
U_null=SMooth(stoks(*,2,1),wx,/edge_truncate)
Q_null=0 & U_null=0

for k=0,Ntarget-1 do begin
stoks(*,1,k)=stoks(*,1,k)-Q_null
stoks(*,2,k)=stoks(*,2,k)-U_null
for j=0,Nexp(k)-1 do begin
cube_stoks(*,1,j,k)=cube_stoks(*,1,j,k)-Q_null
cube_stoks(*,2,j,k)=cube_stoks(*,2,j,k)-U_null

if mode eq 1 then begin
Q_tmp=cube_stoks(*,2,j,k)
U_tmp=-cube_stoks(*,1,j,k)
cube_stoks(*,1,j,k)=Q_tmp
cube_stoks(*,2,j,k)=U_tmp
endif
endfor
endfor
if mode eq 1 then begin
Q_tmp=stoks(*,1,*) & U_tmp=stoks(*,2,*)
stoks(*,1,*)=U_tmp & stoks(*,2,*)=-Q_tmp
endif
writefits,wdir+'avg_stoks.fit',stoks,h
writefits,wdir+'stoks.fit',cube_stoks,hs
wx=50


window,0,xsize=600,ysize=427
T=1
robomean,stoks(Nx/2-Nx/4:Nx/2+Nx/4,1,T)*100,2,0.5,avg_val,rms_val
plot,stoks(*,1,T)*100,yrange=[-5,5],xst=1,charsize=1,$
	title='Q='+string(avg_val,format='(F7.2)')+$
	'!9 +!3'+string(rms_val,format='(F5.2)')

oplot,[0,Nx],[0,0],linestyle=2
robomean,stoks(Nx/2-Nx/4:Nx/2+Nx/4,2,T)*100,2,0.5,avg_val,rms_val
plot,stoks(*,2,T)*100,yrange=[-5,5],xst=1,charsize=1,$
	title='U='+string(avg_val,format='(F7.2)')+$
	'!9 +!3'+string(rms_val,format='(F5.2)')

oplot,[0,Nx],[0,0],linestyle=2
;***********
FLUX=fltarr(Nx,3) & avg_Q=fltarr(NX)  & avg_U=avg_Q
FLUX(*,*)=stoks(*,0,*)

T=2

;box=5
;avg_Q=LOWESS(findgen(Nx),Q(*,T),box,2,2)

avg_Q(*)=stoks(*,1,T)-Q_null
;avg_U=LOWESS(findgen(Nx),U(*,T),box,2,2)
avg_U(*)=stoks(*,2,T)-U_null
avg_P=sqrt(avg_Q^2+avg_U^2)




avg_FI=fltarr(Nx)
for x=0,Nx-1 do  avg_FI(x)=calc_atan(avg_Q(x),avg_U(x))/2; & R=where(avg_FI gt 178) & avg_FI(R)=avg_FI(R)-180

Window,2,xsize=600,ysize=900
!P.background=2^16-1 & !P.color=0
!P.multi=[0,1,5]
!P.charsize=2
R=where(wave gt 5000 and wave lt 6000)
RR=where(wave gt 4500 and wave lt 7500,IND); & PRINT,IND
RR=findgen(Nx)
;PRINT,AVG_p
robomean,avg_Q(R)*100,3,0.5,QQ,err_QQ
robomean,avg_U(R)*100,3,0.5,UU,err_UU
robomean,avg_P(R)*100,3,0.5,P,err_P
robomean,avg_FI(R) ,3,0.5,F,err_F
plot,wave(RR),FLUX(RR,T),xst=1,$
title=name(T)+'   PAslit='+string(PA(T),format='(F7.1)')
plot,wave(RR),avg_Q(RR)*100,xst=1,yrange=[-5,5],yst=1
plot,wave(RR),avg_U(RR)*100,xst=1,yrange=[-5,5],yst=1
plot,wave(RR),avg_P(RR)*100,xst=1,yrange=[0,10],yst=1,$
	title='P(V)='+string(P,format='(F5.2)')+'!9 +!3'+string(err_P,format='(F5.2)')
plot,wave(RR),avg_FI(RR),xst=1,title='!9P!3(V)='+string(F,format='(F7.1)')+'!9 +!3'+string(err_F,format='(F4.1)'),$
subtitle=wdir

print,QQ,err_QQ
print,UU,err_UU
print,P,err_P
print,F,err_F
;print,star_name
;input=''
;read,input,prompt='tabulated angle='
;zero=[calc_atan(qq,uu)/2, calc_atan(qq,-uu)/2,calc_atan(-qq,uu)/2,calc_atan(uu,qq)/2,calc_atan(-uu,qq)/2,calc_atan(uu,-qq)/2];,format='(6f8.1)'
;zero=zero-PA_slit+FLOAT(input)
;print,'PA_tab',FLOAT(input)
;print,'PA_slit',PA_slit
;if mode eq 0 then PA_0=177.9 else PA_0=219.3
;print,name_out(0),star_name,zero(0),-F+PA_0+PA_slit
;END

fin:
END



;нуль пункт старого WOLL2
zero=[182.7,180.2,179.0,177.1,176.9,177.9,179.2,178.3,178.3,177.7,175.6,177.1,176.7,176.3,$
	 177.43,177.2,179.1,179.2,178.3,178,176.9,173.7,180.6,175.6,182.4,176.3,176.5,177.8]
	 robomean,zero,3,0.5,mean_zero,tms_zero
print,mean_zero,rms_zero
xhist=findgen(21)*2+160
window,5
!P.multi=[0,1,1]
plot,xhist,histogram(zero,min=160,max=200,bin=2),psym=10
end
;нуль пункт yjdjuj WOLL2
zero=[226.4,220.4,218.4,213.2,205.3,232.8,222.9,230.4,211.6,212.1,212.1,224,223]
	; 177.43,177.2,179.1,179.2,178.3,178,176.9,173.7,180.6,175.6,182.4,176.3,176.5,177.8]
	 robomean,zero,3,0.5,mean_zero,tms_zero
print,mean_zero,rms_zero
xhist=findgen(11)*10+160
window,5
!P.multi=[0,1,1]
plot,xhist,histogram(zero,min=160,max=260,bin=10),psym=10
END
T=0
ISM=[-0.35,-0.14]/100.
ISM=[0,0]
avg_Q(*)=stoks(*,1,T)-Q_null-ISM(0)
;avg_U=LOWESS(findgen(Nx),U(*,T),box,2,2)
avg_U(*)=stoks(*,2,T)-U_null-ISM(1)
avg_P=sqrt(avg_Q^2+avg_U^2)
avg_FI=fltarr(Nx)
for x=0,Nx-1 do  avg_FI(x)=calc_atan(avg_Q(x),avg_U(x))/2; & R=where(avg_FI gt 178) & avg_FI(R)=avg_FI(R)-180

Window,4,xsize=600,ysize=900,xpos=600,ypos=0
!P.background=2^16-1 & !P.color=0
!P.multi=[0,1,5]
!P.charsize=2
R=where(wave gt 5000 and wave lt 6000)
RR=where(wave gt 4500 and wave lt 8000)
robomean,avg_Q(R)*100,1,0.5,Q,err_Q
robomean,avg_U(R)*100,1,0.5,U,err_U
robomean,avg_P(R)*100,1,0.5,P,err_P
robomean,avg_FI(R) ,1,0.5,F,err_F
plot,wave(RR),FLUX(RR,T),xst=1,title=name(T)+'   PAslit='+string(PA(T),format='(F7.1)')
plot,wave(RR),avg_Q(RR)*100,xst=1,yrange=[-5,5],yst=1
plot,wave(RR),avg_U(RR)*100,xst=1,yrange=[-5,5],yst=1
plot,wave(RR),avg_P(RR)*100,xst=1,yrange=[0,10],yst=1,$
	title='P(V)='+string(P,format='(F5.2)')+'!9 +!3'+string(err_P,format='(F5.2)')
plot,wave(RR),avg_FI(RR),xst=1,title='!9P!3(V)='+string(F,format='(F7.1)')+'!9 +!3'+string(err_F,format='(F4.1)'),$
subtitle=wdir
print,Q,err_Q
print,U,err_U
print,P,err_P
print,F,err_F
END
window,4
!P.multi=[0,1,1]
;avg_vector=robust_estimation(avg_P,wx=10,tresh=1 )
avg_P=median(avg_P,10)
plot,wave(RR),avg_P(RR),xst=1
oplot,wave(RR),LOWESS(wave(RR),avg_P(RR),50,3,3),color=1e5,thick=2

end
