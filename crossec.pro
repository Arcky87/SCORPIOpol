; cross-section through image
function CROSSEC,im,x,y,r=r,xout=xout,yout=yout,dh=dh,median=median,mask=mask,sigma=sigma

if n_params() lt 3 then begin
 print,'Uses CROSSEC,im,x,y[,r=r,xout=xout,yout=yout,dh=dh,sigma=sigma,/median,mask=mask]'
 return,-1
endif

if not(keyword_set(dh)) then dh=0
; cross-section  direction
 dx=float(x(1)-x(0))
 dy=float(y(1)-y(0))
 pm=max(abs([dx,dy]))
 ax=dx/pm
 ay=dy/pm

ss=size(im)
 ; width
 pa=atan(dy,dx)+!pi/2
 xrec=(findgen(2*dh+1)-dh)*cos(pa)
 yrec=(findgen(2*dh+1)-dh)*sin(pa)

pv=fltarr(pm)
sigma=fltarr(pm)
xout=fltarr(pm)
yout=fltarr(pm)

FOR p=0,pm-1 DO BEGIN
 xc=x(0)+ax*p
 yc=y(0)+ay*p
IF xc le ss(1) and yc le ss(2) and xc ge 0 and yc ge 0THEN BEGIN
 xout(p)=xc
 yout(p)=yc
 imcut=im(round(xc+xrec),round(yc+yrec))
 if keyword_set(mask) then begin
  recmask=where (imcut gt mask,num)
  if num gt 0 then begin
         if not(keyword_set(median)) then pv(p)=total(imcut(recmask))/num  else pv(p)=median(imcut(recmask))
  endif else  pv(p)=mask
   ;RMS
  if num gt 1 then sigma(p)=total( (imcut(recmask)-pv(p))^2 )/num

 endif else begin ; if  MASK not set
  if not(keyword_set(median)) then pv(p)=total(imcut)/(1+2*dh)  else pv(p)=median(imcut)
   ; RMS
  if dh gt 0 then sigma(p)=total( (imcut-pv(p))^2 )/(1+2*dh)
 endelse

ENDIF
IF xc lt  0 then xout(p)=0
IF yc lt  0 then yout(p)=0
IF xc gt  ss(1) then xout(p)=ss(1)
IF yc gt  ss(2) then yout(p)=ss(2)
ENDFOR
r=sqrt((xout-x(0))^2+(yout-y(0))^2)
sigma=sqrt(sigma)
return,pv
END


