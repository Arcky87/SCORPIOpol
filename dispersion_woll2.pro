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
wx=20; default 6

FOR y=0,Ny-1 do begin
  spectra=neon(*,y)
;точное определение позиций линий (переписал для гауссианы)
    for j=0,N_lines-1 do begin
      x_start = pos(j) - wx > 0
      x_end = pos(j) + wx < (Nx - 1)

    ;  x_region = [x_start:x_end]

      ; if (x_end - x_start + 1) lt 5 then begin ;(старый метод)
      ;   ypos=max(spectra(pos(j)-2*wx:pos(j)+2*wx),Nmax)
      ;   xpos(y,j)=pos(j)-2*wx+Nmax
      ;   p=goodpoly(x(xpos(y,j)-wx:xpos(y,j)+wx),spectra(xpos(y,j)-wx:xpos(y,j)+wx),2,2,fit)
      ;   xpos(y,j)=-p(1)/p(2)/2
      ;   continue
      ; endif

    ; Извлекаем данные для фитирования
      x_fit = x[x_start:x_end]
      y_fit = spectra[x_start:x_end]
     ; if j eq 0 then  cgplot,x_fit, y_fit

      ; Оценка начальных параметров
      y_max = max(y_fit, i_max)
      y_min = min(y_fit)

      x0_guess = x_fit[i_max]
      amp_guess = y_max - y_min
      sigma_guess = 2.0 ; начальное приближение для ширины
      const_guess = y_min
      estimates = [amp_guess, x0_guess, sigma_guess, const_guess]

      success = 0
      a_params = 0

      catch, error_status

      if (error_status eq 0) then begin
      y_fitted = GAUSSFIT(x_fit, y_fit, a_params, $
                                ESTIMATES=estimates, $
                                NTERMS=4)
   ;  if j eq 0 then begin
   ;   cgoplot,x_fit, y_fitted, color='red'
   ;   cgtext,10,10, 'line ' + string(y,format='(I4)'),/device
   ;    wait,0.1
   ;  endif

   ;   print, 'a_params num',n_elements(a_params),' a_params 2', a_params[2]
      ;stop
        if (n_elements(a_params) ge 4) and (a_params[2] gt 0) then begin
        success = 1
        xpos(y, j) = a_params[1]
        endif
      endif
    if (success eq 0) then begin
          ; ; Если фитирование неудачно, используем параболу
          start_pix1 = round(pos(j)) - 2
          end_pix1 = round(pos(j)) + 2
          if start_pix1 lt 0 then start_pix1 = 0
          if end_pix1 ge Nx then end_pix1 = Nx - 1

          ypos = max(spectra(start_pix1:end_pix1), Nmax)
          x_rough = start_pix1 + Nmax
          ;ypos = max(spectra(pos(j)-2:pos(j)+2), Nmax)
          ;x_rough = pos(j) - 2 + Nmax
          start_pix2 = x_rough - 2
          end_pix2 = x_rough + 2
          if start_pix2 lt 0 then start_pix2 = 0
          if end_pix2 ge Nx then end_pix2 = Nx - 1
          p = goodpoly(x(start_pix2:end_pix2), spectra(start_pix2:end_pix2), 2, 2, fit)
          ;p = goodpoly(x(x_rough-2:x_rough+2), spectra(x_rough-2:x_rough+2), 2, 2, fit)
          if (p[2] ne 0) then begin
          xpos(y, j) = -p(1)/p(2)/2
          endif else begin
            xpos(y, j) = x_rough
          endelse
          ; 
          print,'No success!', success
    endif

    catch, /cancel
  endfor

 ;       ypos = max(spectra(pos(j)-2:pos(j)+2), Nmax)
 ;       x_rough = pos(j) - 2 + Nmax
 ;       p = goodpoly(x(x_rough-2:x_rough+2), spectra(x_rough-2:x_rough+2), 2, 2, fit)
 ;       xpos(y, j) = -p(1)/p(2)/2
 ;     endelse
 ;   endfor

;построение дисперсионной кривой по высоте щели
D(y,*)=goodpoly(xpos(y,*),wave,N_deg,tresh,Yfit,newx,newy)
ENDFOR

;формирование модели  дисперсионной кривой
D_original = D ; для тестов алгоритма сглаживания
for k=0,N_deg do begin
f=goodpoly(findgen(Ny),D(*,k),2,2,Yfit) ; --- original
D(*,k)=Yfit
endfor

;определение эффективных длин волн
ray=['angle = 0 deg','angle = 90 deg','angle = 45 deg','angle = 135 deg']
if keyword_set(plot)  then begin
!P.background=2^24-1 & !P.color=0
 cgdisplay,wid=plot+10,xsize=1800,ysize=1000,title=ray(plot-1),xpos=10+dpos*(plot-1),ypos=10+dpos*(plot-1)
 cgplot,[0,Ny],[0,N_lines+1],/nodata,ycharsize=1e-5,$
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
 cgoplot,(obs_wave-wave(k))/10+k+1,psym=8,symsize=0.6,COLOR='red'
  cgoplot,[0,Ny-1],[1,1]*(k+1),linestyle=2
   endif
robomean,obs_wave-wave(k),3,0.1,avg_err,rms_err
;rms_err=stdev(obs_wave-wave(k),avg_err)
if ABS(avg_err) gt 3*rms_err  $
 then error=string(avg_err,format='(F6.2)') else error=''
if keyword_set(plot) then $
;xyouts,Ny*1.02,(k+1)-0.1,string(wave(k),format='(F8.2)')$
cgtext,Ny*1.02,(k+1)-0.1,string(wave(k),format='(F8.2)')$
+'!9 +!3'+string(rms_err,format='(F5.2)')+error,color='black',/data
endfor
;график дисперсии и отождествления линий
D_mean=total(D,1)/Ny
disp=0 & for j=1,N_deg do disp=disp+D_mean(j)*j*findgen(Nx)^(j-1)

spectra=neon(*,Ny/2) & spectra=sqrt(spectra-min(spectra))
if keyword_set(plot) then begin
cgplot,spectra,xst=1,yrange=[0,max(spectra)*1.5],yst=1,$
position=[0.53,0.5,0.99,0.99],/norm,/noerase
endif
for j=0,N_lines-1 do begin
  start_pix = round(pos(j)) - wx
  end_pix = round(pos(j)) + wx
  if start_pix lt 0 then start_pix = 0
  if end_pix ge Nx then end_pix = Nx - 1
  ypos=max(spectra(start_pix:end_pix),Nmax)
;ypos=max(spectra(pos(j)-wx:pos(j)+wx),Nmax)
xpos=start_pix+Nmax
;xpos=pos(j)-wx+Nmax
if keyword_set(plot) then cgtext,xpos+10,ypos,$
 string(wave(j),format='(F9.2)'),/data,orientation=90,charsize=0.75
endfor
if keyword_set(plot) then cgplot,disp,xst=1,yst=1,COLOR='blue',$
position=[0.53,0.051,0.99,0.45],/norm,/noerase ;,yrange=[FIX(min(disp)),RFIX(max(disp))]

return,D
end
