;dispersion
function dispersion_WOLL2,neon,table,N_DEG=N_deg,TRESH=tresh,PLOT=plot,DX=dx

dpos=40
plotsym,0
if not(keyword_set(dx)) then dx=0
a=size(neon)  & Nx=a(1)  & Ny=a(2)
x=findgen(Nx)
if not(keyword_set(N_deg)) then N_deg=3
if not(keyword_set(tresh)) then tresh=3

D=dblarr(Ny,N_deg+1)

;чтение таблицы линий
b=size(table) & N_lines=b(2)
pos=fltarr(N_lines) & wave=pos
wave(*)=table(0,*) & pos(*)=table(1,*)
xpos=fltarr(Ny,N_lines)
wx=3


	FOR y=0,Ny-1 do begin
spectra=neon(*,y)
;точное определение позиций линий
	for j=0,N_lines-1 do begin
	ypos=max(spectra(pos(j)-2*wx:pos(j)+2*wx),Nmax)
	xpos(y,j)=pos(j)-2*wx+Nmax
	p=goodpoly(x(xpos(y,j)-wx:xpos(y,j)+wx),spectra(xpos(y,j)-wx:xpos(y,j)+wx),2,2,fit)
	xpos(y,j)=-p(1)/p(2)/2
	endfor
;построение дисперсионной кривой по высоте щели
D(y,*)=goodpoly(xpos(y,*),wave,N_deg,tresh,Yfit,newx,newy)
	ENDFOR
;формирование модели  дисперсионной кривой

	for k=0,N_deg do begin
f=goodpoly(findgen(Ny),D(*,k),2,2,Yfit)
D(*,k)=Yfit
	endfor
;определение эффективных длин волн
ray=['angle = 0 deg','angle = 90 deg','angle = 45 deg','angle = 135 deg']
if keyword_set(plot)  then begin
!P.background=2^24-1 & !P.color=0
	cgdisplay,wid=plot+10,xsize=1400,ysize=600,title=ray(plot-1),xpos=dpos*(plot-1),ypos=dpos*(plot-1)
		plot,[0,Ny],[0,N_lines+1],/nodata,ycharsize=1e-5,$
			position=[0.01,0.05,0.4,0.99],/norm,xst=1,yst=1
	 			endif
	for k=0,N_lines-1 do begin
obs_wave=fltarr(Ny)
	for y=0,Ny-1 do begin
		fit=0 & for j=0,N_deg do fit=fit+D(y,j)*xpos(y,k)^j
			obs_wave(y)=fit
				endfor
robomean,obs_wave,3,0.5,avg_wave,rms_wave

if keyword_set(plot) then begin
	oplot,(obs_wave-wave(k))/10+k+1,psym=8,symsize=0.6
		oplot,[0,Ny-1],[1,1]*(k+1),linestyle=2
			endif
robomean,obs_wave-wave(k),5,0.5,avg_err,rms_err
;rms_err=stdev(obs_wave-wave(k),avg_err)
if ABS(avg_err) gt 3*rms_err  $
 then error=string(avg_err,format='(F6.2)') else error=''
if keyword_set(plot) then $
xyouts,Ny*1.02,(k+1)-0.1,string(wave(k),format='(F8.2)')$
	+'!9 +!3'+string(rms_err,format='(F5.2)')+error,/data
endfor
;график дисперсии и отождествления линий
D_mean=total(D,1)/Ny
disp=0 & for j=1,N_deg do disp=disp+D_mean(j)*j*findgen(Nx)^(j-1)

spectra=neon(*,Ny/2) & spectra=sqrt(spectra-min(spectra))
if keyword_set(plot) then begin
plot,spectra,xst=1,yrange=[0,max(spectra)*1.5],yst=1,$
position=[0.53,0.5,0.99,0.99],/norm,/noerase
endif
for j=0,N_lines-1 do begin
ypos=max(spectra(pos(j)-wx:pos(j)+wx),Nmax)
xpos=pos(j)-wx+Nmax
if keyword_set(plot) then xyouts,xpos+10,ypos,$
	string(wave(j),format='(F9.2)'),/data,orientation=90,charsize=0.75
endfor
if keyword_set(plot) then plot,disp,xst=1,yst=1,$
position=[0.53,0.051,0.99,0.45],/norm,/noerase ;,yrange=[FIX(min(disp)),RFIX(max(disp))]

return,D
end
