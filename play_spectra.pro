pro show_pol,wave,obj,Q,U,P,PA,z,WIN=win,WTITLE=wtitle,XRNG=xrng
if not(keyword_set(wtitle)) then wtitle=''
if not(keyword_set(win)) then win=1
window,win,xsize=600,ysize=750,title=wtitle,xpos=win*20,ypos=win*20
!P.multi=[0,1,5]
!P.charsize=1
N=N_elements(wave)
if not(keyword_set(xrng)) then xrng=[wave(0),wave(N-1)]
;xrng=[4500,7200]

dz=0
plot,wave,obj/max(obj)*10,xst=1,xrange=xrng,position=[0.05,0.80,0.99,0.99],$
	/norm,xcharsize=1e-5,yrange=[0,11],yst=1,psym=10
	oplot,[1,1]*6563*(1+z),[0,1]*1e6,color=3e5
plot,wave,Q,xst=1,xrange=xrng,position=[0.05,0.61,0.99,0.80],$
	/norm,/noerase,xcharsize=1e-5,psym=10
	oplot,[1,1]*6563*(1+z),[-1,1]*1e6,color=3e5
plot,wave,U,xst=1,xrange=xrng,position=[0.05,0.42,0.99,0.61],$
	/norm,/noerase,xcharsize=1e-5,psym=10
	oplot,[1,1]*6563*(1+z),[-1,1]*1e6,color=3e5
plot,wave,P,xst=1,xrange=xrng,position=[0.05,0.23,0.99,0.42],$
	/norm,/noerase,xcharsize=1e-5,psym=10
	oplot,[1,1]*6563*(1+z),[0,1]*1e6,color=3e5
RR=where(wave gt xrng(0) and wave lt xrng(1))
robomean,PA(RR),3,0.5,avg_PA,rms_PA
plot,wave(RR),PA(RR),xst=1 ,yst=1,position=[0.05,0.04,0.99,0.23],$
	/norm,/noerase,psym=10,yrange=5*[-1,1]*rms_PA+avg_PA
oplot,[1,1]*6563*(1+z-dz) ,[0,1]*1e6,color=3e5

RW=where(wave gt 5000 and wave lt 6000)
robomean,P(RW),3,0.5,avg_P,rms_P
print, avg_P,rms_P
robomean,PA(RW),3,0.5,  avg_PA,rms_PA
print, avg_PA,rms_PA
end
;*************************
function  angle_calculation,Q,U
Teta=ATAN(U/Q)*180/!PI
if Q gt 0 and U ge 0 then TETA_0=0
if Q gt 0 and U lt 0 then TETA_0=360
if Q lt 0 then TETA_0=180
if Q eq 0 and U gt 0 then TETA_0=90
if Q eq 0 and U lt 0 then TETA_0=270
return, (TETA+TETA_0)/2
END
;******************************
function  limit_integration,w
;вычисление пределов интегрирования в окне w
M=41
ind=indgen(M)
ww=intarr(2,M) & for j=0,1 do ww(j,*)=ind & ww=reform(ww,2*M)
;print,ww,format='(40I3)'
wr=ww(0:M-1)
wl=[0,shift(ww,-1)] & wr=[0,wr(0:M-1)]
;print,wl,format='('+string(M,format='(I2)')+'I3)'
;print,wr,format='('+string(M,format='(I2)')+'I3)'
;print,w,format='('+string(M,format='(I2)')+'I3)'
;print,wl+wr,format='('+string(M,format='(I2)')+'I3)'
return,[wl(w),wr(w)]
end
;**************************************
function remove_trend,vector
N=N_elements(vector)
f=goodpoly(findgen(N),vector,1,3,Yfit)
vector=(vector-Yfit)+Yfit(N/2)
return,vector
end
;***********************************************************
;play spectra]
set_plot,'WIN'
dz=0
ISM=[0,0]
;dir='h:\red_data.pol\Sy1\MCG+08-11-011_151106\' & ISM=[-0.0,-0.3] & kp=1
;dir='h:\red_data.pol\Sy1\Akn120_140324\' 	& ISM=[0,0] & kp=1
;dir='h:\red_data.pol\Sy1\Mkn1095_141123\' 	& ISM=[0,0] & kp=1
;dir='h:\red_data.pol\Sy1\3C120_131106\'	 	& ISM=[-0.3,-0.3] & kp=1
;dir='h:\red_data.pol\Sy1\Mrk1502_141122\'	& ISM=[-0.6,-0.25] & dz=0.001
;dir='h:\red_data.pol\Sy1\1Zw1_131103\'	& ISM=[0,0.] & dz=0.0605 & kp=0.4
;dir='h:\red_data.pol\Sy1\Mkn1148_141120\'	& ISM=[-0.2,-0.3]  & dz=0.001 &kP=1
;dir='h:\red_data.pol\Sy1\Mkn335_131109\'	& ISM=[-0.,0.]  & dz=0.00 &kP=1
;
;dir='h:\red_data.pol\Sy1\NGC3227_140325\'	& ISM=[-0.0,0.0]  & dz=0.00
;dir='h:\red_data.pol\Sy1\NGC3227_140325\'	& ISM=[-0.0,0.0]  & dz=0.00
;dir='h:\red_data.pol\Sy1\Mkn110_151207\'	& ISM=[0.0,0.0]  & dz=0.00
;dir='h:\red_data.pol\Sy1\Mkn704_150325\'	& ISM=[-0.2,0.2]  & dz=0.00
;dir='h:\red_data.pol\Sy1\PG0844+349_141121\'
;dir='h:\red_data.pol\Sy1\Mkn79_141020\'		&  ISM=[-0.0,0.0]  & dz=0.00 & kp=1
;dir='h:\red_data.pol\Sy1\Mkn79_151106\'		&  ISM=[-0.0,0.0]  & dz=0.00 & kp=1
;dir='h:\red_data.pol\Sy1\Mkn79_151207\'		&  ISM=[-0.0,0.0]  & dz=0.00 & kp=1
;dir='h:\red_data.pol\Sy1\Mkn6_131104\'		&  ISM=[-0.3,0]  & dz=0.00 & kp=1
;
;dir='h:\red_data.pol\Sy1\Mkn231_150318\'		&  ISM=[-0.0,0.0]  & dz=0.00 & kp=1
;dir='h:\red_data.pol\Sy1\NGC4593_150318\'		& ISM=[0.0,0.0]  & dz=-0.001
;;dir='h:\red_data.pol\Sy1\NGC4593_160308\'		& ISM=[0.0,0.0]  & dz=-0.0
;dir='h:\red_data.pol\Sy1\3C273_140325\'		& ISM=[-1.0,0.0]  & dz=-0.0
;dir='h:\red_data.pol\Sy1\NGC4051_140324\'	& ISM=[0.0,0.0] & dz=-0.0 &kp=1
;dir='h:\red_data.pol\Sy1\NGC4151_140523\'	& ISM=[0.0,0.0] & dz=-0.0
;dir='h:\red_data.pol\Sy1\NGC4151_160307\'	& ISM=[0.5,-1] & dz=-0.0

;dir='h:\red_data.pol\Sy1\Mkn744_150329\'	& ISM=[0.0,0.0] & dz=-0.0
;
;dir='h:\red_data.pol\Sy1\Mkn876_150326\'	& ISM=[0,0] & dz=-0.0 & kp=1
;dir='h:\red_data.pol\Sy1\Mkn876_160408\'	& ISM=[0,0] & dz=-0.0
;dir='h:\red_data.pol\Sy1\Mkn841_140529\'	& ISM=[0.6,0.8] & dz=-0.0
dir='d:\Sy1\Mkn817_140529\'	& ISM=[0,0] & dz=-0.0 & kp=1
;dir='D:\Sy1\sbs1419+538_190216\' 	& ISM=[0,0] & dz=-0.0 & kp=1
;dir='d:\Sy1\Mkn335_131109\'	& ISM=[0,0] & dz=-0.0 & kp=1
;dir='D:\SCORPIO\red_data.pol\NGC4151_160307\'	& ISM=[0,0] & dz=-0.0 & kp=1
;dir='d:\Sy1\NGC4151_160307\'	& ISM=[0,0] & dz=-0.0
;dir='h:\red_data.pol\Sy1\NGC5548_140325\'	& ISM=[0,0] & dz=-0.0
;dir='h:\red_data.pol\Sy1\NGC5548_150325\'	& ISM=[0,0] & dz=-0.0
;dir='h:\red_data.pol\Sy1\Mkn668_150324\'	& ISM=[0,0] & dz=-0.0
;dir='h:\red_data.pol\Sy1\IRAS13349+2438_140523\' & ISM=[0,5] & dz=-0.0
;
;dir='h:\red_data.pol\Sy1\IRAS03450_141020\'	& ISM=[-0.00,-0.0] & dz=-0.0 & kp=1
;dir='h:\red_data.pol\Sy1\3C390_140529\'		& ISM=[0,0] & dz=-0.0
;dir='h:\red_data.pol\Sy1\3C445_141113\'	& ISM=[0,02] & dz=-0.0 & kP=1
;dir='h:\red_data.pol\Sy1\3C445_131105\'	& ISM=[0,02] & dz=-0.0 & kP=1
;dir='d:\Sy1\Mkn304_131109\'	& ISM=[0,-2] & dz=-0.0 & kP=0.3
;dir='d:\Sy1\Mkn509_141021\'	& ISM=[0,0] & dz=-0.0 & kP=1
;dir='h:\red_data.pol\Sy1\PG1700+518_160405\'	& ISM=[0,0] & dz=-0.0 & kP=1
;dir='d:\Sy1\IRAS13349+2438_140523\'  	& ISM=[0,0] & dz=-0.0 & kP=1
;dir='h:\red_data.pol\Sy1\Mkn304_131109\'
;dir='h:\red_data.pol\Sy1\NGC7469_151210\'& ISM=[0,-2] & dz=-0.0 & kP=0.3
;dir='h:\red_data.pol\Sy1\E1841+643_140324\'& ISM=[0,0] & dz=-0.0 & kP=1
;dir='h:\red_data.pol\Sy1\E1841+643_161124\'& ISM=[-1,0] & dz=-0.0 & kP=1
;dir='h:\red_data.pol\Sy1\E1841+643_17012\'& ISM=[0,0] & dz=-0.0 & kP=1
;dir='h:\red_data.pol\Sy1\Akn564_131109\'& ISM=[0,0] & dz=-0.0 & kP=0.5
;dir='h:\red_data.pol\Sy1\Mkn1501_141120\'& ISM=[0,0] & dz=-0.0 & kP=1
;
;dir='h:\red_data.pol\Sy1\IRAS13349+2438_140523\'

avg_spectra=readfits(dir+'avg_spectra.fit',h,/silent)
;
PA=fltarr(3) & for j=0,2 do PA(j)=sxpar(h,'PA'+string(j+1,format='(I1)'))
name=strarr(3) & for j=0,2 do name(j)=sxpar(h,'NAME'+string(j+1,format='(I1)'))
Nexp=intarr(3) & for j=0,2 do Nexp(j)=sxpar(h,'NUMEXP'+string(j+1,format='(I1)'))
Nx=sxpar(h,'NAXIS1') & N=Nx
wave=findgen(Nx)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
z=sxpar(h,'Z') & print,'Object redshift:  ', z & z=z+dz

RW=where(wave gt 6000 and wave lt 8000)

;определение относительного смещения поляризационных каналов
dx=fltarr(4)
neon=readfits(dir+'neon_lin.fts',hn,/silent) & a=size(neon)
neon=total(neon(*,*,*),2)
Rn=where(wave gt 6550 and wave lt 6650)
x=findgen(N)
for k=0,3 do begin
fit=gaussfit(x(Rn),neon(Rn,k),G)
dx(k)=g(1)
endfor
dx=dx-total(dx)/4
print,'Shifts from neon:  ', dx
;dx=[-0.0,0.0,0, 0]
for k=0,3 do avg_spectra(*,k,0)=shift_s(avg_spectra(*,k,0),-dx(k))

;проверка типа анализатора по плоскому полю
tmp=avg_spectra
num=indgen(4)
slope=slope_flat(dir)  & print,'slope',slope
if slope eq 1 then avg_spectra=tmp(*,reverse(num),*)

;кривая чувствительности (?)
if file_test(dir+'sent.fts') eq 1 then sent=readfits(dir+'sent.fts',/silent) else sent=1

;чтение параметров из хэдера
PA=fltarr(3) & for j=0,2 do PA(j)=sxpar(h,'PA'+string(j+1,format='(I1)'))
name=strarr(3) & for j=0,2 do name(j)=sxpar(h,'NAME'+string(j+1,format='(I1)'))
Nexp=intarr(3) & for j=0,2 do Nexp(j)=sxpar(h,'NUMEXP'+string(j+1,format='(I1)'))
N=sxpar(h,'NAXIS1')
wave=findgen(N)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
RW=where(wave gt 6000 and wave lt 8000)
!P.charsize=1.5
LIN=[0,Nx-1]
Npol=4
x=findgen(N)
type=['obj ',' unpolarized star ',' polarized star ']
avg_stoks=fltarr(N,5,3) ; 0=> flux, 1=>Q, 2=>U, 3=>P, 4=>PA

;1 - определение нуль-пункта и деполяризации?
;cor=corr_2order(avg_spectra(*,*,1),LINEAR=LIN,/plot);,WTITLE=dir+type(1))
;cor - коэф-ты деполяризации
cor=fltarr(N,2)
cor(*,0)=avg_spectra(*,0,1)/avg_spectra(*,1,1) & cor(*,1)=avg_spectra(*,2,1)/avg_spectra(*,3,1)
kQ=LOWESS(x,cor(*,0),N/4,3,3) & kU=LOWESS(x,cor(*,1),N/3,3,3)
	cgdisplay, wid=0
	!p.multi=[0,1,2]
	cgplot, x, cor(*,0), xrange=[0,Nx-1], yrange=[median(cor(*,0))-0.11,median(cor(*,0))+0.11]
		cgoplot, x, kQ, color='orange'
	cgplot, x, cor(*,1), xrange=[0,Nx-1], yrange=[median(cor(*,1))-0.11,median(cor(*,1))+0.11]
		cgoplot, x, kU, color='orange'
cor(*,0)=kQ & cor(*,1)=kU
;cor(*,0)=replicate(1,N) & cor(*,1)=replicate(1,N)

;проверка по нулевому стандарту
S0=calc_stoks(avg_spectra(*,*,1),/plot,corr=cor)
print, 'Rough Q,U of unpol std, %:  ', median(S0(*,1))*100, median(S0(*,2))*100
;нуль-пункт поляризации (вообще он =0)
avg_null=fltarr(N,2) & for k=0,1 do avg_null(*,k)=LOWESS(x,S0(*,k+1),N/2,1,1)

;2 - вычисление параметров Стокса по среднему спектру
avg_corr=fltarr(N,2,3)
for j=0,2 do begin
;		;коррекция 2-го порядка
	;	avg_corr(*,*,j)=corr_2order(avg_spectra(*,*,j),LINEAR=[LIN(0),LIN(1)])
	;вычисление параметров Стокса
	avg_stoks(*,0:2,j)=calc_stoks(avg_spectra(*,*,j),corr=cor)
	;учет нуль-пункта
	avg_stoks(*,1:2,j)=avg_stoks(*,1:2,j)-avg_null(*,0:1)
	;учет межзвездной поляризации
		if j eq 0 then begin
			avg_stoks(*,1,j)=(avg_stoks(*,1,j)-ISM(0)/100.)*Kp
			avg_stoks(*,2,j)=(avg_stoks(*,2,j)-ISM(1)/100.)*Kp
		endif
	;вычисление степени поляризации
	avg_stoks(*,3,j)=SQRT(avg_stoks(*,1,j)^2+avg_stoks(*,2,j)^2)
	;вычисление угла плоскости поляризации
		for i=0,N-1 do begin
			;avg_stoks(i,4,j)=PA(j)-angle_calculation(avg_stoks(i,1,j),avg_stoks(i,2,j))+180
			avg_stoks(i,4,j)=PA(j)-angle_calculation(avg_stoks(i,1,j),avg_stoks(i,2,j))+179
			R=where(avg_stoks(*,4,j) gt 180,ind) & if ind gt 1 then avg_stoks(R,4,j)=avg_stoks(R,4,j)-180
			R=where(avg_stoks(*,4,j) lt 45,ind) & if ind gt 1 then avg_stoks(R,4,j)=avg_stoks(R,4,j)+90
		endfor
	avg_stoks(*,0,j)=avg_stoks(*,0,j);/sent
	;подавление выбросов в угле наклона
	f=goodpoly (x(*),avg_stoks(*,4,j),1,1,fit)
	RA=where(ABS(avg_stoks(*,4,j)-fit) gt 40, ind)
	if ind gt 1 then avg_stoks(RA,4,j)=fit(RA)
endfor

RR=where(wave gt 6400 and wave lt 6500) ;& avg_stoks(RR,4,0)=avg_stoks(RR,4,0)+10
Ha=6562.6
xrng=[-1,1]*500+Ha*(1+z)
;xrng=[wave(0),wave(N-1)]
window,2,xsize=1000,ysize=900,title=dir
!P.multi=[0,3,5]
	for k=0,4 do begin
for j=0,2 do begin
	amp=1
	if k gt 0 and k lt 4 then amp=100
		robomean,avg_stoks(RW,k,j)*amp,1,0.5,avg_s,rms_s


	if k gt 0 then titl=string(avg_s)+string(rms_s) ELSE titl=name(j)
	plot,wave,avg_stoks(*,k,j)*amp,xst=1,title=titl,xrange=xrng
	if j eq 0 then oplot,[1,1]*Ha*(1+z),[-1,1]*1e6,color=3e5
	;correction angle
	if k eq 4 and j eq 0 then begin
	angle=LOWESS(wave,avg_stoks(*,k,j),20,3,3)
	noise=avg_stoks(*,k,j)-angle
	oplot,wave,angle,color=3e5,thick=2
;	oplot,[wave(0),wave(N-1)],[1,1]*avg_s
;	RR=where(wave gt 6800 and wave lt 7200)
;	angle(RR)=(angle(RR)-avg_s+1)*2+avg_s
;	oplot,wave,angle,thick=2
;	oplot,wave,angle+noise
	;avg_stoks(*,4,0)=angle+noise/2
	oplot,wave,angle+noise/2,color=3e5,thick=2
	avg_stoks(*,4,0)=angle+noise/2
	endif
	endfor
endfor
;формирование выходного файла объекта
RA=where(ABS(wave-Ha*(1+z)) lt 500, ind)
wave=wave(RA)  & Na=N_elements(wave)
flux=avg_stoks(RA,0,0)
P=avg_stoks(RA,3,0)*100
;P(124:ind-1)=P(124:ind-1)+0.5
PA=avg_stoks(RA,4,0)
;P=P-LOWESS(findgen(ind),P,ind/2,2,2)+0.7
;RS=where(P lt 0) & P(RS)=-P(RS)
;PA(212:ind-1)=PA(212:ind-1)+20
;PA(0:162)=PA(0:162)+15
;PA=PA-LOWESS(findgen(ind),PA,ind/2,2,2)+98
;PA(ind/2:ind-1)=shift(PA(ind/2:ind-1),-5)
;RS=where(PA lt 60) & PA(RS)=100
;d_PA=10
;Q=P*cos(PA*!PI/90)-0.5 & Q=shift(Q,-d_PA)& U=shift(U,d_PA)
;U=P*sin(PA*!PI/90) & U=remove_trend(U)
;;
;
;P=sqrt(Q^2+U^2)/4
;for j=0,Na-1 do PA(j)=angle_calculation(Q(j),U(j)) & R=where(PA gt 90) & PA(R)=PA(R)-90
;PA=PA-180
;P=shift(P,-10) & PA=shift(PA,-10)
Q=P*cos(PA*!PI/90)
U=P*sin(PA*!PI/90)
;U=remove_trend(U)
P=sqrt(Q^2+U^2)
;
;END
;Pa=shift(PA,10)
window,1,xsize=400 ,ysize=1100
!P.multi=[0,1,5]
plot,wave,flux,xst=1 	&	oplot,[1,1]*(1+z)*Ha,[-1,1]*1e6,linestyle=2
plot,wave,Q,xst=1, yrange=[-2,2]		&	oplot,[1,1]*(1+z)*Ha,[-1,1]*1e6,linestyle=2
plot,wave,U,xst=1, yrange=[-2,2]		&	oplot,[1,1]*(1+z)*Ha,[-1,1]*1e6,linestyle=2
;R=where(P gt 1) & P(R)=P(R)-0.3
plot,wave,P,xst=1, yrange=[0,5]		& oplot,[1,1]*(1+z)*Ha,[-1,1]*1e6,linestyle=2
;PA=PA-180
PA=remove_trend(PA); & PA=median(PA,3)
plot,wave,PA,xst=1		& oplot,[1,1]*(1+z)*Ha,[-1,1]*1e6,linestyle=2
;oplot,wave,LOWESS(findgen(Na),PA,20,2,2),color=3e5,thick=3


;запись таблицы результата
obj=str_sep(dir,'\')
obj=str_sep(obj(N_elements(obj)-2),'_')
obj=obj(0)
openw,1,'D:\SCORPIO\red_data.pol\result\'+obj+'.txt'
for k=0,Na-1 do printf,1,wave(k),flux(k),Q(k),U(k),P(k),PA(k),format='(I4,E12.3,3F7.2,F8.1)'
close,1
END

;OLD BLOCK!!
;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;;вычисление параметров Стокса по среднему спектру
;avg_corr=fltarr(N,2,3)
;for j=0,2 do begin
;;коррекция 2-го порядка
;avg_corr(*,*,j)=corr_2order(avg_spectra(*,*,j),LINEAR=[LIN(0),LIN(1)])
;;вычисление параметров Стокса
;avg_stoks(*,0:2,j)=calc_stoks(avg_spectra(*,*,j),corr=avg_corr(*,*,j))
;
;;учет нуль-пункта
;avg_stoks(*,1:2,j)=avg_stoks(*,1:2,j)-avg_null(*,0:1)
;;учет межзвездной поляризации
;if j eq 0 then begin
;avg_stoks(*,1,j)=(avg_stoks(*,1,j)-ISM(0)/100.)*Kp
;
;avg_stoks(*,2,j)=(avg_stoks(*,2,j)-ISM(1)/100.)*Kp
;endif
;;вычисление степени поляризации
;avg_stoks(*,3,j)=SQRT(avg_stoks(*,1,j)^2+avg_stoks(*,2,j)^2)
;
;;вычисление угла плоскости поляризации
;for i=0,N-1 do begin
;
;;avg_stoks(i,4,j)=PA(j)-angle_calculation(avg_stoks(i,1,j),avg_stoks(i,2,j))+180
;avg_stoks(i,4,j)=PA(j)+calc_atan(avg_stoks(i,1,j),avg_stoks(i,2,j))/2+180-317.3;219
;;if avg_stoks(i,4,j) gt 360 then avg_stoks(i,4,j)=avg_stoks(i,4,j)-360
;;if avg_stoks(i,4,j) gt 180 then avg_stoks(i,4,j)=avg_stoks(i,4,j)-180
;;if avg_stoks(i,4,j) lt 0 then avg_stoks(i,4,j)=avg_stoks(i,4,j)+180
;;if avg_stoks(i,4,j) lt 90 then avg_stoks(i,4,j)=avg_stoks(i,4,j)+90
;;;avg_stoks(i,4,j)=avg_stoks(i,4,j)-160
;;avg_stoks(i,4,j)=-avg_stoks(i,4,j)+160-150
;
;;if avg_stoks(i,4,j) lt 45 then avg_stoks(i,4,j)=avg_stoks(i,4,j)+45
;;if avg_stoks(i,4,j) gt 90 then avg_stoks(i,4,j)=avg_stoks(i,4,j)-45
;;avg_stoks(i,4,j)=avg_stoks(i,4,j)+120
;endfor

;avg_stoks(*,0,j)=avg_stoks(*,0,j);/sent
;;подавление выбросов в угле наклона
;f=goodpoly (x(*),avg_stoks(*,4,j),1,1,fit)
;RA=where(ABS(avg_stoks(*,4,j)-fit) gt 40, ind)
;if ind gt 1 then avg_stoks(RA,4,j)=fit(RA)
;
;endfor
;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



