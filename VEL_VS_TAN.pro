pro kepler_motion,VEL,FI,PLOT=plot,PATH=path







END
;*****************************************
;Rsc=452 & titl=' Akn 120, Rsc=452 ligth days'  ;

;Rsc=963 & titl=' 3C273, Rsc=963 ligth days'  ;
;Rsc=38 & titl=' NGC4051, Rsc=38 ligth days'
;Rsc=1130 & titl='IRAS13349+2438 , Rsc=1130 ligth days'
;Rsc=142 & titl='Mkn335 , Rsc=142 ligth days'
;Rsc=60& titl=' NGC5548, Rsc=60 ligth days'
;Rsc=44 & titl=' NGC4151, Rsc=44 ligth days'
;Rsc=180 & titl=' Mkn817, Rsc=180 ligth days'
 ;
;*****************************************

crsz=1.2
;obj='Mkn110_151207' & Rsc=100 & titl=' Mkn110, Rsc=100 ligth days' & z=0.035291 & dv=0 & dFI=0 & sign=-1
;obj='Mkn79_141020' & Rsc=55 & titl=' Mkn79, Rsc=55 ligth days' & z=0.022189 & dv=0 & sign=-1
;obj='NGC4593_150318' & Rsc=43 & titl=' NGC4593, Rsc=43 ligth days' & z=0.009001 & dv=0 & sign=-1
;obj='PG1700+518_160405' & Rsc=687 & titl='PG1700+518_160405, Rsc=687 ligth days' & z=0.2890 & dv=0 & sign=1 &  Vmax=15e3 & dFI=-15 & Vmin=0 & TRESH=2 & V_min=3e3 &  V_max=10e3
obj='MCG+08-11-011_151106' & Rsc=687 & titl='MCG+08-11-011, Rsc=90 ligth days' & z=0.020484 & dv=0
;obj='3c390_140529' & Rsc=119 & titl='3c390.3, Rsc=119 ligth days' & z=0.0561 & dv=0
obj='PG0844+349_141121' & Rsc=189 & titl='PG0844+349, Rsc=189 ligth days' & z=0.0660 & dv=0 & sign=1 & Vmax=10e3 & dFI=0 & Vmin=0 & TRESH=2 & V_min=1e3 &  V_max=6e3 & TETA=0
;obj='3C120_131106' & Rsc=226 & titl='3C120, Rsc=226 ligth days' & z=0.0330 & dv=0
;obj='E1841+643' & Rsc=2161 & titl='E1841+643, Rsc=2161 ligth days' & z=0.297


;������ ������� RES
dir='h:\red_data.pol\Sy1\'

tab=read_table(dir+obj+'\'+obj+'.res')
Ntab=N_elements(tab)-1
wave=fltarr(Ntab) & avg_FI=wave & avg_Q=wave  & avg_U=wave & avg_F=wave  & P=wave &
for  j=0,Ntab-1 do begin

wave(j)=FLOAT(strmid(tab(j+1),0,4))
avg_F(j)=FLOAT(strmid(tab(j),5,12))/FLOAT(strmid(tab(j),40,6))*100
avg_Q(j)=FLOAT(strmid(tab(j+1),17,5))/100
avg_U(j)=FLOAT(strmid(tab(j+1),29,5))/100
avg_FI(j)=strmid(tab(j+1),53,7)
avg_FI(j)=calc_atan(avg_Q(j),avg_U(j))/2
;print,wave(j),avg_Q(j),avg_U(j),avg_FI(j)
endfor


;������ ������� TXT
;dir='h:\red_data.pol\AGN\publish_AGN\'
;Ntab=numlines(dir+obj+'.txt')
;tab=fltarr(2,Ntab)r
;wave=fltarr(Ntab) & avg_FI=wave
;openr,1,dir+obj+'.txt'
;readf,1,tab
;close,1
;wave(*)=tab(0,*) & avg_FI(*)=tab(1,*)

;����������� ������ ���������� �� 5100
cont=5100*(1+z)  & d_cont=50*(1+z)
Rc=where(ABS(wave-cont) lt d_cont)
robomean,avg_F(Rc),3,0.5,avg_cont,rms_cont
print,avg_cont,rms_cont
;avg_F=avg_F/avg_cont
Ha=(6563)*(1+z)
d_wave=Vmax*Ha/3e5

R=where(ABS(wave-Ha) lt d_wave,Nw)
wave=wave(R)  & avg_FI=avg_FI(R) & avg_Q=avg_Q(R) & avg_U=avg_U(R) & avg_F=avg_F(R)
robomean,avg_Q*avg_F,3,0.5,Q_cont
robomean,avg_U*avg_F,3,0.5,U_cont
;Q_cont=0 & U_cont=0
Vel=((wave-Ha)/Ha*3e5)+dV



set_plot,'WIN'
crsz=1
window,0,xsize=700,ysize=1100
!p.multi=[0,1,1]
;����������� ������ ����������
N=N_elements(vel)
CF=goodpoly([vel(0:N/4),vel(N-N/4:N-1)],[avg_F(0:N/4),avg_F(N-N/4:N-1)],1,3)
cont_F=CF(0)+CF(1)*vel
;avg_F=avg_F-cont_F
;Flux
order_F=FIX(ALOG10(max(avg_F)))
plot,vel,avg_F/10^order_F,xst=1,charsize=crsz,psym=10,$
	yrange=[0.001,1]*max(avg_F)*1.1/10^order_F,yst=1,$
	position=[0.08,0.78,0.99,0.96],xcharsize=1e-5
	xyouts,0.025,0.87,'Total Flux, 10!U'+string(order_F,format='(I1)')+'!N, ADU',/norm,align=0.5,orientation=90
	oplot,vel,cont_F/10^order_F,linestyle=2
	oplot,[0,0],[-1,1]*1e3,linestyle=1
;Q-Stokes
order_Q=FIX(ALOG10(max(ABS(avg_Q*avg_F)))) & print,order_Q

plot,vel,avg_Q*avg_F/10.^order_Q,xst=1,charsize=crsz,psym=10,$
	position=[0.08,0.60,0.99,0.78],xcharsize=1e-5,/noerase
	xyouts,0.025,0.69,'Q-Stoks Flux, 10!U'+string(order_Q,format='(I1)')+'!N, ADU',/norm,align=0.5,orientation=90
	oplot,[min(vel),max(vel)],[1,1]*Q_cont/10^order_Q,linestyle=2
	oplot,[0,0],[-1,1]*1e3,linestyle=1
;U-Stokes
order_U=FIX(ALOG10(max(ABS(avg_U*avg_F))))
plot,vel,avg_U*avg_F/10^order_U,xst=1,charsize=crsz,psym=10,$
	yrange=[-1,1]*max(ABS(avg_U))*avg_F*1.3/10^order_U,yst=1,$
   	position=[0.08,0.42,0.99,0.60],xcharsize=1e-5,/noerase
	oplot,[min(vel),max(vel)],[1,1]*U_cont/10^order_U,linestyle=2
	oplot,[0,0],[-1,1]*1e3,linestyle=1
	xyouts,0.025,0.51,'U-Stoks Flux, 10!U'+string(order_U,format='(I1)')+'!N, ADU',/norm,align=0.5,orientation=90
print,'QU cont', Q_cont,U_cont
;Polarized Flux
Q=(avg_Q*avg_F-Q_cont) & U=(avg_U*avg_F-U_cont)
order_P=FIX(ALOG10(max(sqrt(Q^2+U^2))))
plot,vel,sqrt(Q^2+U^2)/10^order_P,xst=1,charsize=crsz,psym=10,$
	yrange=[min(sqrt(Q^2+U^2)),max(sqrt(Q^2+U^2))]*1.1/10^order_P,yst=1,$
	position=[0.08,0.24,0.99,0.42],xcharsize=1e-5,/noerase
	oplot,[0,0],[-1,1]*1e3,linestyle=1
	xyouts,0.025,0.33,'Polarized Flux, 10!U'+string(order_P,format='(I1)')+'!N, ADU',/norm,align=0.5,orientation=90
;Angle of Polarization
plot,vel,avg_FI,xst=1,charsize=crsz,psym=10,$
	xtitle='Velocity, km/s',$
	position=[0.08,0.06,0.99,0.24] ,/noerase
	robomean,avg_FI,3,0.5,level_FI
	oplot,[min(vel),max(vel)],[1,1]*level_FI,linestyle=2
	oplot,[0,0],[-1,1]*1e3,linestyle=1
	xyouts,0.025,0.15,'Angle of Polarization, deg',/norm,align=0.5,orientation=90

;******************************************************
window,1,xsize=550,ysize=550,xpos=620,ypos=0
!P.multi=[0,1,1]
Ampl=300
plot,[-1,1]*Ampl,[-1,1]*Ampl,xst=1,yst=1,/nodata,position=[0.1,0.1,0.995,0.995],/norm
;R=where(vel gt 2e3 and U lt -0.2 ) & Q(R)=0 & U(R)=0

RPL=where(vel gt 0 and vel lt  V_min) & RPH= where(vel gt  V_min and vel lt  V_max)
RNL=where(vel lt 0 and vel gt -V_min) & RNH= where(vel lt -V_min and vel gt -V_max)
FI=findgen(32)*!PI/16 & FI=[FI,FI(0)]
USERSYM,cos(FI),sin(FI)

oplot,Q(RPL),U(RPL),psym=8
oplot,Q(RPH),U(RPH),psym=8,symsize=1.5,color=1e5

USERSYM,cos(FI),sin(FI),/fill
oplot,Q(RNL),U(RNL),psym=8
oplot,Q(RNH),U(RNH),psym=8 ,symsize=1.5,color=1e5
oplot,[-1,1]*Ampl,[0,0],linestyle=2
oplot,[0,0],[-1,1]*Ampl,linestyle=2
PA=2*level_FI-90
oplot,[-sin((PA)/180.*!PI),sin((PA)/180.*!PI)]*ampl,[1,-1]*ampl,linestyle=1
R_e=1.1* ampl  & a=6
x_e=cos(FI)*R_e & y_e=sin(FI)*R_e/a-50
ROT_XY, X_e, Y_e,PA , 0, 0, X_e, Y_e,  /degrees
oplot,x_e,y_e,linestyle=1
;end





robomean,avg_FI,1,0.5,FI_cont
Vel=((wave-Ha)/Ha*3e5)+dV
avg_FI=avg_FI-FI_cont

f=goodpoly(vel,avg_FI,1,2,fit)
;avg_FI=avg_FI-fit
avg_FI=avg_FI+dFI
;R=where(avg_fi lt -10) & avg_FI(R)=-3
;Avg_FI(R)=Avg_FI(R)+10


set_plot,'WIN'
;VE__VS_TAN

key_plot=0
if key_plot eq 0 then Window,2,xsize=600,ysize=900
!P.multi=[0,1,2]
if key_plot eq 1 then begin
set_plot,'PS'
device,file=dir+obj+'\Angle_VS_Vel.ps',xsize=16,ysize=23,xoffset=2,yoffset=4
endif
if key_plot eq 2 then begin
set_plot,'PRINTER'
device,scale_factor=1
endif

!P.multi=[0,1,2]
;goto,fin

;����������� LOG(v/c))_VS_log(tan(FI))
RL=where(Vel lt -Vmin,ind_RL)  & RR=where(Vel gt Vmin,ind_RR)
print,ind_RL,ind_RR
avg_FI=avg_FI*sign

T=TAN(avg_FI*!PI/180)

;rms_T=rms_FI/180.*!PI/(cos(rms_FI/180.*!PI))^2
VR=VEL(RR)/3e5  & TR=-T(RR) & ;rms_TR=rms_T(RR)
R=where(TR gt 0.0) & TR=TR(R) & VR=VR(R)  & ;rms_TR=rms_TR(R)
VL=-VEL(RL)/3e5  & TL=T(RL) & ;rms_TL=rms_T(RR)
R=where(TL gt 0,ind)
if ind gt 0 then begin
TL=TL(R)  & VL=VL(R) & ;rms_TL=rms_TL(R)
endif

A=-2.5 & B=-0.5

;f=goodpoly([ALOG10(TL),ALOG10(TR)],[ALOG10(VL),ALOG10(VR)],1,2,Fit)
;print,'f: ',f(0),f(1)

LogT=[ALOG10(TL),ALOG10(TR)] & LogV=[ALOG10(VL),ALOG10(VR)]
LogV_Kepler=A+B*LogT

robomean,LogV-LogV_Kepler,tresh,0.5,avg_A,rms_A

avg_A=A+avg_A
print,avg_A,rms_A




if Rsc ne 0 then begin
;Log_M_BH=ALOG10(1.78)+10+2*avg_A+ALOG10(Rsc)
Log_M_BH=ALOG10(0.58)+10+2*avg_A +ALOG10(Rsc)
err_Log_M_BH=2*rms_A
print,Log_M_BH,err_Log_M_BH
endif

LogV_Kepler=avg_A+B*LogT
;
plot,Vel,avg_FI,xst=1 ,charsize=crsz,yrange=[-10,10]*8,$
	title=titl,xtickformat='(I6)',psym=10,xtitle='Velocity, km/s',ytitle='!7D!9P!3'

xyouts,7000,80,'broad H!7a!3',charsize=2,align=0.5
oplot,[0,0] ,[-1,1]*100,linestyle=2
oplot,[-1,1]*2e4,[0,0],linestyle=2
print,'sign ',sign
if key_plot eq  0 then col=4e4 else col=200
RPOS=WHERE(Vel gt Vmin)  & RNEG=WHERE(Vel lt -Vmin)
M_BH=Log_M_BH

FI_Kepler=(ATAN(10.0^(-2*avg_A)*Vel(RPOS)^2/9E10/2)*180/!PI-90)/2
oplot,Vel(RPOS),-FI_Kepler*sign,color=col,thick=10
FI_Kepler=(90-ATAN(10.0^(-2*avg_A)*Vel(RNEG)^2/9E10)*180/!PI)/2

oplot,Vel(RNEG),-FI_Kepler*sign ,color=col,thick=10
oplot,Vel,avg_FI,psym=10 & oplot,[-1,1]*2e4,[0,0],linestyle=2
oplot,Vel(RR),avg_FI(RR),psym=10,thick=2
oplot,Vel(RL),avg_FI(RL),psym=10,thick=2

TETA=findgen(32)*!PI/16 & teta=[teta,teta(0)]
USERSYM,cos(teta),sin(teta)
R=where(ABS(ALOG10(VR)-avg_A-B*ALOG10(TR)) lt tresh*rms_A)

plot,ALOG10(TR(R)),ALOG10(VR(R)),psym=8,yrange=[-3,0],charsize=crsz,$
	title='Log(V/c) = a + 0.5 Log(tan!9P!3)',$
	xtickinterval=1,xtitle='Log(tan!9P!3)',ytitle='Log(V/c)',xrange=[-3,0.5],xst=1
oplot,[-2.8],[-2.5],psym=8
xyouts,[-2.85],[-2.52],'( ) - Velocity > 0'
R=where(ABS(ALOG10(VL)-avg_A-B*ALOG10(TL)) lt tresh*rms_A)
USERSYM,cos(teta),sin(teta),/fill
oplot,ALOG10(TL(R)),ALOG10(VL(R)),psym=8
oplot,[-2.8],[-2.7],psym=8
xyouts,[-2.85],[-2.72],'( ) - Velocity < 0'
index=sort(LogT)
oplot,LogT(index),LogV_Kepler(index),linestyle=2
xyouts,-1.5,-0.3,'a = '+string(avg_A,format='(F5.2)')+'!9+!3'+string(rms_A,format='(F4.2)'),$
	charsize=crsz*1.5,align=0.5
if Rsc ne 0 then xyouts,-1.5,-0.7,'Log(M!DBH!N/M!D!9n!N!3)='$
	+string(Log_M_BH,format='(F5.2)')+'!9+!3'+string(err_Log_M_BH,format='(F5.2)'),$
	charsize=crsz*1.5,align=0.5
if key_plot gt 0  then device,/close
set_plot,'WIN'
end