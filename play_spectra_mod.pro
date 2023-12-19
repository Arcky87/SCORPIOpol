pro show_pol,wave,obj,Q,U,P,PA,z,WIN=win,WTITLE=wtitle,XRNG=xrng
if not(keyword_set(wtitle)) then wtitle=''
if not(keyword_set(win)) then win=1
window,win,xsize=600,ysize=750,title=wtitle,xpos=win*20,ypos=win*20
!P.multi=[0,1,5]
!P.charsize=2
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
;вычисление пределов интегрировани€ в окне w
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

;***********************************************************
;play spectra

;dir='h:\red_data.pol\AGN\Mkn335_131109	\'
;dir='h:\red_data.pol\Sy1\Mkn509_141021\'
;dir='h:\red_data.pol\Sy1\MCG+08-11-011_151106\'
;dir='h:\red_data.pol\Sy1\Mkn231_150318\'
;dir='h:\red_data.pol\Sy1\Mkn876_160408\'
dir='d:\Sy1\Mkn817_140529\'
;dir='h:\red_data.pol\Sy1\IRAS03450_141020\'
;dir='h:\red_data.pol\Sy1\Akn120_140324\'
;dir='h:\red_data.pol\Sy1\Mkn110_151207\
;dir='h:\red_data.pol\Sy1\3C120_131106\'
;dir='h:\red_data.pol\Sy1\3C445_141113\'
;dir='h:\red_data.pol\Sy1\Mkn704_141123\'
;dir='h:\red_data.pol\Sy1\NGC4593_150318\'
;dir='h:\red_data.pol\Sy1\Mkn841_140529\'
;dir='h:\red_data.pol\Sy1\PG0844+349_141121\'
;dir='h:\red_data.pol\Sy1\PG1700+518_160405\'
;dir='h:\red_data.pol\Sy1\NGC4151_140523\
;dir='h:\red_data.pol\Sy1\NGC4051_140324\
;dir='h:\red_data.pol\Sy1\Mkn1148_141120\'
;dir='h:\red_data.pol\Sy1\Mrk1502_141122\'
;dir='h:\red_data.pol\Sy1\IRAS13349+2438_140523\'
spectra=readfits(dir+'spectra.fit',h)
;avg_spectra=readfits(dir+'avg_spectra.fit',h,/silent)
tmp=spectra
num=indgen(4)
slope=slope_flat(dir)  & print,'slope',slope
if slope eq 1 then spectra=tmp(*,reverse(num),*,*)
z=sxpar(h,'Z')
sent=readfits(dir+'sent.fts',/silent)
 N=sxpar(h,'NAXIS1') & Npol=sxpar(h,'NAXIS2') & Ntarget=sxpar(h,'NAXIS4')
avg_spectra=fltarr(N,Npol,Ntarget)
PA=fltarr(Ntarget)   & for j=0,Ntarget-1 do   PA(j)=sxpar(h,'PA'+string(j+1,format='(I1)'))
name=strarr(Ntarget) & for j=0,Ntarget-1 do name(j)=sxpar(h,'NAME'+string(j+1,format='(I1)'))
Nexp=intarr(Ntarget) & for j=0,Ntarget-1 do Nexp(j)=sxpar(h,'NUMEXP'+string(j+1,format='(I1)'))
;формирование среднего спектра
		for i=0,Ntarget-1 do begin
	for j=0,Npol-1 do begin
for k=0,N-1 do begin
robomean,spectra(k,j,0:Nexp(i)-1,i),1,0.5,mean
avg_spectra(k,j,i)=mean;median(spectra(k,j,0:Nexp(i)-1,i))
endfor
	endfor
			endfor

wave=findgen(N)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
RW=where(wave gt 6000 and wave lt 8000)
!P.charsize=1.5
;***************************
;goto,cont
;***************************
Nmax=1800
Npol=4
x=findgen(N)
type=['obj ',' unpolarized star ',' polarized star ']
avg_corr=fltarr(N,2,3)
;определение средней коррекции по объектам
for j=0,2 do begin
;коррекци€ 2-го пор€дка
avg_corr(*,*,j)=corr_2order(avg_spectra(*,*,j),LINEAR=[0,Nmax],WTITLE=dir+type(j),numwin=j);,/plot)

endfor

;определение нуль- пункта
S0=calc_stoks(avg_spectra(*,*,1),corr=avg_corr(*,*,1));,/plot)
avg_null=fltarr(N,2) & for k=0,1 do avg_null(*,k)=LOWESS(x,S0(*,k+1),N/2,1,1)

;вычисление робастных оценок параметров —токса
avg_stoks=fltarr(N,5,3) ; 0=> flux, 1=>Q, 2=>U, 3=>P, 4=>PA
rms_stoks=avg_stoks
stoks=fltarr(N,3,Nexp(0),Ntarget)
;**************
wl=1 & wr=1 & tresh=1
;**************
		FOR j=0,2 DO BEGIN
   	for k=0,Nexp(j)-1 do begin
stoks(*,*,k,j)=calc_stoks(spectra(*,*,k,j),corr=avg_corr(*,*,j))
	endfor

for k=wl,N-1-wr do begin
	for i=0,2 do begin
		robomean,stoks(k,i,0:Nexp(j)-1,j),tresh,0.5,val,rms
		avg_stoks(k,i,j)=val
		rms_stoks(k,i,j)=rms
	endfor
endfor
	avg_stoks(*,1:2,j)=avg_stoks(*,1:2,j)-avg_null(*,0:1)
		ENDFOR
window,0,xsize=500,ysize=900
!P.multi=[0,1,3]
j=0
for k=0,2 do plot,wave,avg_stoks(*,k,j),xst=1

END

;avg_stoks(*,0:2,j)=calc_stoks(avg_spectra(*,*,j),corr=avg_corr(*,*,j))
;
;;учет нуль-пункта
;avg_stoks(*,1:2,j)=avg_stoks(*,1:2,j)-avg_null(*,0:1)
;;вычисление степени пол€ризации
;avg_stoks(*,3,j)=SQRT(avg_stoks(*,1,j)^2+avg_stoks(*,2,j)^2)
;
;;вычисление угла плоскости пол€ризации
;for i=0,N-1 do begin
;
;avg_stoks(i,4,j)=PA(j)-angle_calculation(avg_stoks(i,1,j),avg_stoks(i,2,j))+180
;;avg_stoks(i,4,j)=PA(j)-angle_calculation(-avg_stoks(i,2,j),-avg_stoks(i,1,j))+180
;if avg_stoks(i,4,j) gt 360 then avg_stoks(i,4,j)=avg_stoks(i,4,j)-360
;if avg_stoks(i,4,j) gt 180 then avg_stoks(i,4,j)=avg_stoks(i,4,j)-180
;
;
;;endfor
;
;avg_stoks(*,0,j)=avg_stoks(*,0,j)/sent
;;;подавление выбросов в угле наклона
;;f=goodpoly (x(*),avg_stoks(*,4,j),1,1,fit)
;;RA=where(ABS(avg_stoks(*,4,j)-fit) gt 30, ind)
;;if ind gt 1 then avg_stoks(RA,4,j)=fit(RA)
;
;;endfor
;Ha=6562.6
;window,2,xsize=1500,ysize=900,title=dir
;!P.multi=[0,3,5]
;	for k=0,4 do begin
;for j=0,2 do begin
;	amp=1
;	if k gt 0 and k lt 4 then amp=100
;		robomean,avg_stoks(RW,k,j)*amp,1,0.5,avg_s,rms_s
;
;
;	if k gt 0 then titl=string(avg_s)+string(rms_s) ELSE titl=name(j)
;	plot,wave,avg_stoks(*,k,j)*amp,xst=1,title=titl
;	if j eq 0 then oplot,[1,1]*Ha*(1+z),[-1,1]*1e6,color=3e5
;
;	endfor
;endfor
;;cont:
;;вычисление параметров —токса по дл€ каждой экспозиции
;
;spectra=readfits(dir+'spectra.fit',h,/silent)
;a=size(spectra)
;N=a(1) & Npol=a(2) & Ntarget=a(4)
;stoks=fltarr(N,5,Nexp(0),Ntarget)
;
;
;
;		for j=0,Ntarget-1 do begin
;	for k=0,Nexp(j)-1 do begin
;;коррекци€ 2-го пор€дка
;
;stoks(*,0:2,k,j)=calc_stoks(spectra(*,*,k,j),corr=avg_corr(*,*,j))
;
;;учет нуль-пункта
;stoks(*,1:2,k,j)=stoks(*,1:2,k,j)-avg_null(*,0:1)
;;вычисление степени пол€ризации
;stoks(*,3,k,j)=SQRT(stoks(*,1,k,j)^2+stoks(*,2,k,j)^2)
;
;;вычисление угла плоскости пол€ризации
;for i=0,N-1 do begin
;stoks(i,4,k,j)=PA(j)-angle_calculation(stoks(i,1,k,j),stoks(i,2,k,j))+180
;;avg_stoks(i,4,j)=PA(j)-angle_calculation(-avg_stoks(i,2,j),-avg_stoks(i,1,j))+180
;if stoks(i,4,k,j) gt 360 then stoks(i,4,k,j)=stoks(i,4,k,j)-360
;if stoks(i,4,k,j) gt 180 then stoks(i,4,k,j)=stoks(i,4,k,j)-180
;endfor
;stoks(*,0,k,j)=stoks(*,0,k,j)/sent
; 	endfor
; 	endfor
;writefits,dir+'stoks.fts',stoks,h
;;cont:
;
;
;;робастные оценки пол€ризационных параметров
;
;stoks=readfits(dir+'stoks.fts',h)
;
;N=sxpar(h,'NAXIS1') & Npol=sxpar(h,'NAXIS2') & Nobj=sxpar(h,'NAXIS4')
;z=sxpar(h,'Z')
;  PA=fltarr(Nobj) & for j=0,Nobj-1 do   PA(j)=sxpar(h,'PA'+string(j+1,format='(I1)'))
;name=strarr(Nobj) & for j=0,Nobj-1 do name(j)=sxpar(h,'NAME'+string(j+1,format='(I1)'))
;Nexp=intarr(Nobj) & for j=0,Nobj-1 do Nexp(j)=sxpar(h,'NUMEXP'+string(j+1,format='(I1)'))
;
;wave=findgen(N)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
;RW=where(wave gt 6000 and wave lt 8000)
;!P.charsize=1.5
;avg_stoks=fltarr(N,Npol,Nobj)  & rms_stoks=avg_stoks
;;определение пределов интегрировани€
;w=3 & ww=limit_integration(w)
;tresh=1
;		for i=0,Nobj-1 do begin
;	for j=0,Npol-1 do begin
;for k=ww(0),N-1-ww(1) do begin
;robomean,stoks(k-ww(0):k+ww(1),j,0:Nexp(i)-1,i),tresh,0.5,avg_val,rms_val
;avg_stoks(k,j,i)=avg_val
;rms_stoks(k,j,i)=rms_val
;endfor
;	endfor
;		endfor
;print,'ok'
;Ha=6562.6
;window,3,xsize=1500,ysize=900,title=dir
;!P.multi=[0,3,5]
;	for k=0,Npol-1 do begin
;for j=0,Nobj-1 do begin
;	amp=1
;	if k gt 0 and k lt 4 then amp=100
;		robomean,avg_stoks(RW,k,j)*amp,1,0.5,avg_s,rms_s
;
;
;	if k gt 0 then titl=string(avg_s)+string(rms_s) ELSE titl=name(j)
;	plot,wave,avg_stoks(*,k,j)*amp,xst=1,title=titl
;	if j eq 0 then oplot,[1,1]*Ha*(1+z),[-1,1]*1e6,color=3e5
;
;	endfor
;endfor
;
;END
