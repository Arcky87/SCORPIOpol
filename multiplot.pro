;test mult]i plot
pro multiplot,NUM=num,TARGET=target,XRANGE=xrange,ERR=ERR
COMMON GRAPH,spectra,stoks,wave,dir
if not(keyword_set(num)) then num=[1,1,1,1,1,1]
RN=where(num eq 1,indN)

Nplot= indN+1

if not(keyword_set(err)) then err=0
if not(keyword_set(target))then T=0 ELSE T=TARGET
S=N_elements(wave)
R=where(spectra.data(*,0) gt 0, F)
if keyword_set(xrange) then $
R= where(spectra.data(*,0) gt xrange(0) AND spectra.data(*,0) lt xrange(1))
;
!P.charsize=1
if not(keyword_set(Nplot)) then Nplot=5
if not(keyword_set(xtitle)) then xtitle=''
if not(keyword_set(ytitle)) then ytitle=strarr(Nplot)

X_margin=[0.1,0.995]
Y_margin=[0.05,0.97]
dY=(Y_margin(1)-Y_margin(0))/Nplot

for k=0,Nplot-1 do begin
if k eq Nplot-1 then xchar=1 ELSE xchar=1e-5
if k eq 0 then noerase=0  ELSE noerase=1

if k eq 0 then begin
;определение пределов вывода по ординате
max_spectra=max(spectra.data(0:F-1,T+1))
order=RFIX(ALOG10(max_spectra))

plot,spectra.data(R,0),spectra.data(R-1,T+1)/10.^order, $
	yrange=[-0.02,1.2]*max_spectra/10.^order,yst=1,psym=10,$
	title=STRCOMPRESS('Date '+spectra.date+', object '+spectra.name(T)+$
	', exposure '+string(spectra.exptime(T),format='(I4)')+$
	' s, PA '+string(spectra.PA(T),format='(F6.1)')+' deg'),$
	noerase=noerase,xcharsize=xchar,xst=1,$
	position=[X_margin(0),Y_margin(1)-dY*(k+1),X_margin(1),Y_margin(1)-dY*k],/norm
xyouts,0.05,Y_margin(1)-dY*(k+1)+dY/2,'Flux, 10!U'+string(order-2,format='(I1)')+'!N, ADU/px',/norm,orientation=90,align=0.5
endif

if k gt 0 then begin
if strmid(stoks.ytitle(RN(k-1)),0,1) eq 'F' then begin
order=RFIX(ALOG10(stoks.yrange(1,k-1,T)))
stoks.ytitle(0)='Flux, 10!U'+string(order,format='(I1)')+'!N, ADU/px'
endif  ELSE order=0
plot,wave(0:S-1),stoks.value(0:S-1,RN(k-1),T)/10.^order,position=[X_margin(0),Y_margin(1)-dY*(k+1),X_margin(1),Y_margin(1)-dY*k],/norm,$
		noerase=noerase,xcharsize=xchar,xtitle=xtitle,xrange=xrange,xst=1,$
	yrange=stoks.yrange(*,RN(k-1),T)/10.^order,yst=1,psym=10

if err gt 0  then OPLOTERR, wave(0:S-1),stoks.value(0:S-1,RN(k-1),T)/10.^order,stoks.error(0:S-1,RN(k-1),T)/10.^order,psym=10
xyouts,total(xrange)/2,(stoks.yrange(1,RN(k-1),T)-stoks.yrange(0,RN(k-1),T))*0.8+stoks.yrange(0,RN(k-1),T),stoks.param(RN(k-1),T),align=0.5
xyouts,0.05,Y_margin(1)-dY*(k+1)+dY/2,stoks.ytitle(RN(k-1)),/norm,orientation=90,align=0.5
;ytitle=stoks.ytitle(k-1)
endif

endfor
xyouts,total(X_margin)/2,0.02,'Wavelength, '+STRING(197B),/norm,align=0.5
xyouts,0,0,dir,/norm
fin:
end


;******************************************************
pro robust_estimate_stoks,wdir,WOLL,WX=wx,TRESH=tresh,TARGET=target,XRNG=xrng,DPA=dPA,ISM=ISM,CORR=corr
COMMON GRAPH,spectra,stoks,wave,dir
dir=wdir
if not(keyword_set(dPA)) then dPA=0
if not(keyword_set(ISM)) then ISM=[0,0]
;******************************************************
set_plot,'WIN'

stoks_in=readfits(wdir+'stoks.fit' ,head)
Nx=sxpar(head,'NAXIS1') & Ntarget=3
lambda=findgen(Nx)*sxpar(head,'CDELT1')+sxpar(head,'CRVAL1')

;—“–” “”–A ѕ≈–≈ћ≈ЌЌќ…  SPECTRA
Nmax=5000  & Ntarget=3
spectra={spectra,data:fltarr(Nmax,Ntarget+1),$
     		title:strarr(3),name:strarr(Ntarget),date:strarr(1),$
     		Nexp:intarr(3),exptime:fltarr(Ntarget),PA:fltarr(Ntarget)}  ; j=2
;—“–” “”–A ѕ≈–≈ћ≈ЌЌќ…  STOKS
Smax=2000  & Ntarget=3
stoks={stoks,value:fltarr(Smax,6,Ntarget),error:fltarr(Smax,6,Ntarget),yrange:fltarr(2,6,Ntarget),$
     		mean:fltarr(6,Ntarget),rms:fltarr(6,Ntarget),param:strarr(6,Ntarget),ytitle:strarr(6)}
stoks.ytitle=['Flux, ADU','Q-Stoks, %','U-stoks, %','V-Stoks','Degree polarization, %','Angle polarization, deg']
Npol=4
;чтение параметров наблюдений

spectra.data(0:Nx-1,0)=lambda
spectra.date=sxpar(head,'DATE-OBS')

for k=0,Ntarget-1 do begin
spectra.name(k)=sxpar(head,'NAME'+string(k+1,format='(I1)'))
spectra.exptime(k)=sxpar(head,'EXPTIME'+string(k+1,format='(I1)'))
spectra.PA(k)=sxpar(head,'PA'+string(k+1,format='(I1)'))
spectra.Nexp(k)=sxpar(head,'NUMEXP'+string(k+1,format='(I1)'))
endfor

Npos=float(Nx)/wx-1
xpos=findgen(Npos)*wx+wx/2
wave=lambda(xpos)
V=fltarr(Ntarget)
CASE WOLL OF
'WOLL1':begin
		spectra.data(0:Nx-1,1:3)=stoks_in(*,0,0:2)
		Nexp=sxpar(head,'NAXIS3')-Ntarget
		V(*)=total(stoks_in(*,3,0:2),1)/Nx; & print,V
		for k=0,3 do begin
 		for j=0,Npos-1 do begin
 		for i=0,Ntarget-1 do begin
 		IF K EQ 3 and V(i) EQ 0 THEN GOTO,OUT
		if i eq 0 then robomean,stoks_in(xpos(j)-wx/2:xpos(j)+wx/2,k,3:Nexp-1),tresh,0.5,avg_val,rms_val  ELSE $
		robomean,stoks_in(xpos(j)-wx/2:xpos(j)+wx/2,k,i),tresh,0.5,avg_val,rms_val
		stoks.value(j,k,i)=avg_val
		stoks.error(j,k,i)=rms_val
		OUT:
		endfor & endfor & endfor
		;вычисление степени пол€ризации
		stoks.value(*,4,*)=sqrt(stoks.value(*,1,*)^2+stoks.value(*,2,*)^2)
		stoks.error(*,4,*)=sqrt(stoks.error(*,1,*)^2+stoks.error(*,2,*)^2)/sqrt(2)
		stoks.value(*,0,*)=stoks.value(*,0,*)*stoks.value(*,4,*)
		stoks.error(*,0,*)=stoks.value(*,0,*)*stoks.error(*,4,*)
		;вычисление угла плоскости пол€ризации
		;stoks.value(*,5,*)=calc_atan(stoks.value(*,1,*),stoks.value(*,2,*))/2
		stoks.value(*,5,*)=atan(stoks.value(*,2,*)/stoks.value(*,1,*))/2*180/!PI
		for k=0,2 do begin
		;RF=WHERE(stoks.value(*,5,k) gt 180, ind) & if ind gt 0 then stoks.value(RF,5,k)=stoks.value(RF,5,k)-180
		endfor
		stoks.error(*,5,*)=28.4*stoks.error(*,4,*)


		END
'WOLL2': begin
		for j=0,Ntarget-1 do begin
		for k=0,Nx-1 do begin
		robomean,stoks_in(k,0,0:spectra.Nexp(j)-1,j),3,0.5,avg_val,rms_val
		spectra.data(k,j+1)=median(stoks_in(k,0,0:spectra.Nexp(j)-1,j))
		;spectra.data(k,j+1)=avg_val
		endfor & endfor
;анализ нулевого и спектрофотометрического стандарта
	  	;window,2
	  	;!P.multi=[0,1,4]
	  	null=TOTAL(stoks_in(*,*,0:spectra.Nexp(1)-1,1),3)/spectra.Nexp(1)
	  	for k=0,Npol-1 do begin
	  ;	plot,null(*,k),xst=1
	  	null(*,k)=LOWESS(findgen(Nx),null(*,k),Nx/8,2,3)
	  ;	oplot,null(*,k),thick=1,color=1e5
	  	endfor
		;учет инструментальной пол€ризации
		for i=0,Ntarget-1 do begin
		for j=0,spectra.Nexp(i)-1 do begin
		for k=1,Npol-1 do stoks_in(*,k,j,i)=stoks_in(*,k,j,i)-null(*,k)
		endfor & endfor


;робастные оценки параметров —токса в окне
		for k=0,Npol-1 do begin
 		for j=0,Npos-1 do begin
 		for i=0,Ntarget-1 do begin
 		IF k EQ 3 and V(i) EQ 0 THEN GOTO,OUT2
		robomean,stoks_in(xpos(j)-wx/2:xpos(j)+wx/2,k,0:spectra.Nexp(i)-1,i),tresh,0.5,avg_val,rms_val
		stoks.value(j,k,i)=avg_val
		stoks.error(j,k,i)=rms_val
		OUT2:
		endfor & endfor & endfor
			if keyword_set(corr) then begin
		;исправление 2-го пор€дка в объекте
		for k=1,2 do begin
		f=goodpoly(ALOG10(wave(0:Npos/2)),ALOG10(stoks.value(0:Npos/2,k,0)),1,3,FIT)
        fit=f(0)+ALOG10(wave)*f(1)
        fit=10.^fit
        stoks.value(0:Npos-1,k,0)=stoks.value(0:Npos-1,k,0)-LOWESS(findgen(Npos),stoks.value(0:Npos-1,k,0),Npos/4,2,3)+fit
        endfor
        	endif
        ;поворот плоскости пол€ризации
		PA_0=317.3

		for j=0,Ntarget-1 do begin
		TETA=2*(spectra.PA(j)-PA_0+dPA)*!PI/180;-!PI
	;	TETA=!PI/2
		Q=stoks.value(*,1,j)*cos(TETA)-stoks.value(*,2,j)*sin(TETA)
		U=stoks.value(*,1,j)*sin(TETA)+stoks.value(*,2,j)*cos(TETA)
		stoks.value(*,1,j)=Q & stoks.value(*,2,j)=U
		endfor
        ;вычисление степени пол€ризации
		stoks.value(*,4,*)=sqrt(stoks.value(*,1,*)^2+stoks.value(*,2,*)^2)
		stoks.error(*,4,*)=sqrt(stoks.error(*,1,*)^2+stoks.error(*,2,*)^2)
		stoks.value(*,0,*)=stoks.value(*,0,*)*stoks.value(*,4,*)
		stoks.error(*,0,*)=stoks.value(*,0,*)*stoks.error(*,4,*)
		;вычисление угла плоскости пол€ризации
		stoks.value(*,5,*)=calc_atan(stoks.value(*,1,*),stoks.value(*,2,*))/2

	;	RF=WHERE(stoks.value(*,5,0) lt 90, ind) & if ind gt 0 then stoks.value(RF,5,0)=stoks.value(RF,5,0)+90
		;stoks.value(*,5,*)=atan(stoks.value(*,2,*)/stoks.value(*,1,*))/2*180/!PI
		for k=0,2 do begin
	;	RF=WHERE(stoks.value(*,5,k) gt 180, ind) & if ind gt 0 then stoks.value(RF,5,k)=stoks.value(RF,5,k)-90
		endfor
		stoks.error(*,5,*)=28.4*stoks.error(*,4,*)
		end
ENDCASE
;вычисление средних значений  в широкой полосе
stoks_str=['','Q','U','V','P','!9P!3']
stoks_amp=[1,100,100,100,100,1]
min_lim=[0,-2,-2,-2, 0,-90]
max_lim=[1.2, 2, 2, 2, 2, 90]
form_str=['','(F5.2)','(F5.2)','(F5.2)','(F5.2)','(F7.1)']

		Wc=5500  & dW=500 & Rc=where(wave gt Wc-dW AND wave lt Wc+dW)
		for j=0,Ntarget-1 do begin

		for k=1,5 do begin
		stoks.value(*,k,j)=stoks_amp(k)*stoks.value(*,k,j)
		if total(stoks.value(Rc,k,j)) ne 0 then begin
	   	robomean, stoks.value(Rc,k,j),3,0.5,avg_val,rms_val
		stoks.mean(k,j)=avg_val
		stoks.rms(k,j)=rms_val
		stoks.param(k,j)=STRCOMPRESS(stoks_str(k)+'('+string(Wc,format='(I4)')+')='$
			+string(stoks.mean(k,j),format=form_str(k))+'!9 +!3'$
			+string(stoks.rms(k,j),format=form_str(k)))

		endif
		stoks.yrange(*,k,j)=stoks.mean(k,j)+[min_lim(k),max_lim(k)]
		if k eq 4 then 		stoks.yrange(0,k,j)=0
		endfor
		;определение пор€дка
		RW=where(wave gt xrng(0) and wave lt xrng(1))

		order=ALOG(max(stoks.value(RW,0,j)))
	    rms=stdev(stoks.value(RW,0,j),mean)
	    stoks.yrange(*,0,j)=[0,mean+6*rms]

		endfor




cont:
end
path='h:\red_data.pol\TINATIN\'
LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path=path+'LOGS')
print,LOGFILE
WOLL=sxpar(read_table(LOGFILE),'ANALYZER')
WOLL=str_sep(WOLL,'-')
WOLL='WOLL'+WOLL(1)
tmp=str_sep(LOGFILE,'\')
path=tmp(0)
for j=1,N_elements(tmp)-3 do begin
path=path+'\'+tmp(j)
endfor
path=path+'\'
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
wdir=path+wdir(N_elements(wdir)-2)+'\'
print,wdir
print,WOLL
;проверка типа эталона
if WOLL eq 'WOLL2' then begin
flat=readfits(wdir+'flat_i.fts',h)
a=size(flat)
flat=total(flat(a(1)/2-a(1)/4:a(1)/2+a(1)/4,*,*),1)
out=[flat(*,0),flat(*,1),flat(*,2),flat(*,3)]
window,0
plot,out,xst=1,title=sxpar(h,'FILTERS')
endif
XRNG=[4500,8000]
robust_estimate_stoks,wdir,WOLL,WX=10,TRESH=2,TARGET=0,xrng=xrng,dPA=0;,/corr
window,2,xsize=800,ysize=1100
!P.charsize=1.5
multiplot,num=[1,1,1,0,1,1],target=target,xrange=xrng;,/err

end