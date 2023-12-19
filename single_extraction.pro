;test fit_profile


; Gaussian Function
function gauss, x, p, _extra=extra
  sz = size(x)
  if sz(sz(0)+1) EQ 5 then smax = 26D else smax = 13.
 u = (x - p(1))/p(2)
  u=u^2
  mask = u LT (smax^2)  ;; Prevents floating underflow
  if n_elements(p) GE 4 then f = p(3) else f = 0
  if n_elements(p) GE 5 then f = f + p(4)*x
  return,  f + p(0) * mask * exp(-0.5 * temporary(u) * mask)
end

; Lorentzian Function
function lorentz, x, p, _extra=extra
  u = (x - p(1))/p(2)
  u=u^2
  if n_elements(p) GE 4 then f = p(3) else f = 0
  if n_elements(p) GE 5 then f = f + p(4)*x
  return, f + p(0) / (u + 1)
end
; Moffat Function

function moffat, x, p, _extra=extra
  u = (x - p(1))/p(2)
  u=u^2
  if n_elements(p) GE 5 then f = p(4) else f = 0
  if n_elements(p) GE 6 then f = f + p(5)*x
  return, f + p(0) / (u + 1)^p(3)
end

;Sersic Function
function sersic, x, p, _extra=extra

sz = size(x)
  if sz(sz(0)+1) EQ 5 then smax = 26D else smax = 13.
  u = ABS((x - p(1))/p(2))
  tmp=-0.333333+2*p(3)+0.987654E-02/p(3)+1.802861E-3/(p(3)^2)

  u= (u)^(1/p(3))
  ;if n_elements(p) GE 6 then f = p(5) else f=0
  return,  p(0) *exp(-tmp*u)
end
;********************
function fit_profile,x,spectra,POS=pos,WIDTH=width,DEG=deg,YFIT=yfit,PROFILE=profile,$
	PLOT=plot,FIX_POS=fix_pos,FIX_WIDTH=fix_width,FIX_DEG=fix_deg
PROFILE=STRUPCASE(PROFILE)
Np=0
if not(keyword_set(pos)) then pos=N_elements(x)/2
if not(keyword_set(width)) then width=2
if not(keyword_set(deg)) then deg=2
if not(keyword_set(fix_pos)) then fix_pos=0
if not(keyword_set(fix_width)) then fix_width=0
if not(keyword_set(fix_deg)) then fix_deg=0

if PROFILE eq 'GAUSS' OR PROFILE eq 'LORENTZ' then begin
deg=0
fix_deg=1
endif
start_param=FLOAT([max(spectra),pos,width,deg])
sigma=sqrt(ABS(spectra))
parinfo = replicate({fixed:0, limited:[0,0],limits:[0.D,0]}, 4)
parinfo(*).fixed=[0,fix_pos,fix_width,fix_deg]

CASE PROFILE OF
   'GAUSS':	BEGIN
  			start_param(3)=0 & parinfo(3).fixed=1
			param = mpfitfun('GAUSS', x, spectra,sigma,start_param,parinfo=parinfo,Nprint=NP)    ; Fit a function
			Yfit=GAUSS(x,param)
			end
  'MOFFAT':	BEGIN
			param = mpfitfun('MOFFAT', x, spectra,sigma,start_param,parinfo=parinfo,Nprint=NP)    ; Fit a function
			Yfit=MOFFAT(x,param)
			end
 'LORENTZ':	BEGIN
			param = mpfitfun('LORENTZ', x, spectra,sigma,start_param,parinfo=parinfo,Nprint=NP)    ; Fit a function
			Yfit=LORENTZ(x,param)
			end
  'SERSIC': BEGIN
			param = mpfitfun('SERSIC', x, spectra,sigma,start_param,parinfo=parinfo,Nprint=NP)    ; Fit a function
			Yfit=SERSIC(x,param)
			end
ENDCASE
return,param
end
;************************************************************************
function single_extraction,cube,NTYPE=Ntype,PLOT=plot,POS=pos,WX=wx,WY=wy
if not(keyword_set(wx)) then wx=0
if not(keyword_set(wy)) then wy=20
if not(keyword_set(Ntype)) then Ntype=4
a=size(cube)
Nx=a(1) & Ny=a(2) & Npol=a(3) & Nexp=a(4)

print,Nx,Ny,Npol, Nexp
Nray=2
spectra=fltarr(Nx,Nray,Npol,Nexp)
p=fltarr(2)
type=['gauss','moffat','lorentz','sersic','strob']
pol=['0','45','22.5','67.5','0','90']+' deg'
y=findgen(Ny/2) & x=findgen(Nx)
;oпределение начальных параметров
for j=0,Nray-1 do begin
Vy=total(cube(Nx/2-Nx/4:Nx/2+Nx/4,Ny/2*j:Ny/2*(j+1)-1,0,0),1)/(Nx/2+1)
yfit = mpfitpeak(y, Vy, param,/gauss, error=sqrt(Vy))
p(j)=param(1)
endfor

!P.multi=[0,2,Npol]
if keyword_set(pos) then p=[1,1]*pos
;экстракция спектров
Vy=fltarr(Ny/2)
w=20
ray=['  Ordinary ray ','   Extraordinary ray ']

if Ntype eq 4 then str_wy='  halfwidth='+string(wy,format='(I3)')+' px' ELSE str_wy=''
for j=0,Nexp-1 do begin
if keyword_set(plot) then begin
window,2,xsize=1000,ysize=750,retain=1,ypos=200,title='  exp='+string(j)+'   EXTRACTION PROFILE '+type(Ntype)+str_wy+$
', halfwidth window along dispersion ='+string(wx,format='(I3)')+' px'

endif
for i=0,Npol-1 do begin
;for i=0,3 do begin
;fit=fltarr(Nx,Ny/2)
ima=fltarr(Nx,Ny/2)
total_spectra=fltarr(Nx)
err=fltarr(Nx)

for r=0,Nray-1 do begin
flux=fltarr(Nx)
ima(*,*)=cube(*,Ny/2*r:Ny/2*(r+1)-1,i,j)

for k=wx,Nx-1-wx do begin
cursor,xc,yc,0,/data
M=!MOUSE
if M.button eq 4 then begin
wdelete,2
warn= DIALOG_MESSAGE( 'EXTRACTION INTERRUPTED!')
RETURN,-1
endif
Vy(*)=total(ima(k-wx:k+wx,*),1)/(2*wx+1)
if Ntype lt 4 then begin
;yfit = mpfitpeak(y, Vy, param,/moffat, error=sqrt(Vy))
param=fit_profile(y,Vy,POS=p(r),WIDTH=10,YFIT=yfit,PROFILE=type(Ntype));,$
	;PLOT=plot,FIX_POS=fix_pos,FIX_WIDTH=fix_width,FIX_DEG=fix_deg
;fit(k,*)=Yfit
endif
total_spectra(k)=total(Vy(p(r)-wy:p(r)+wy))
if Ntype lt 4 then spectra(k,r,i,j)=total(Yfit) else  spectra(k,r,i,j)=total_spectra(k)
endfor
norm=max(median(spectra(*,r,i,j),20))
if keyword_set(plot) then begin
plot,spectra(*,r,i,j)/norm,xst=1,title='angle='+pol(i),charsize=1e-5,$
	yrange=[-0.2,1.2],yst=1;,color=2^32-1
xyouts,0,1.1,ray(r)+' angle='+pol(i),/data
err(*)=(spectra(*,r,i,j)-total_spectra)/total_spectra
if Ntype lt 4 then begin
robomean,err(Nx/2-Nx/8:Nx/2+Nx/8),3,0.5,avg_err,rms_err
xyouts,Nx/2,-0.13,'avg_err='+string(avg_err*100,format='(F5.2)')+$
	'%, rms_err= '+string(rms_err*100,format='(F4.2)')+'%',/data,align=0.5
	endif
oplot,(spectra(*,r,i,j)-total_spectra)/norm
oplot,[0,Nx],[0,0],linestyle=2
wait,0.1
endif
endfor
endfor
endfor
fin:
return,spectra
end
name='h:\red_data.pol\Arp102b_120417\obj-sky.fts'
;name='h:\red_data.pol\Mkn6\Mkn6_120518\obj-sky.fts'
;name='h:\red_data.pol\3C390.3\3C390_120516\obj-sky.fts'
name='h:\red_data.pol\TINATIN\WD_GJ427_110531\obj-sky.fts'
;name='h:\red_data.pol\TINATIN\WD1756+827_110531\obj-sky.fts'
cube=readfits(name,h)
lambda=sxpar(h,'CRVAL1')+findgen(sxpar(h,'NAXIS1'))*sxpar(h,'CDELT1')
spectra=single_extraction(cube,WY=20,Ntype=4,WX=1);,/plot);,POS=pos
window,2,xsize=1400,ysize=1000
!P.multi=[0,1,2]
j=0
SP=TOTAL(spectra(*,0,*,0),3)
plot,lambda,SP/max(SP),xst=1,title=sxpar(h,'NAME1'),xrange=[4500,7200],charsize=2,yrange=[0,1.2],yst=1,$
	ytitle='Relative intensity', xtitle='Wavelength, A'
V=((spectra(*,0,4,0)-spectra(*,1,4,0))/(spectra(*,0,4,0)+spectra(*,1,4,0)))/2-$
  ((spectra(*,0,5,0)-spectra(*,1,5,0))/(spectra(*,0,5,0)+spectra(*,1,5,0)))/2
  plot,lambda,V*100,xst=1,yrange=[-1,1],ytitle='Circular polarization, %',xtitle='Wavelength, A',xrange=[4500,7200],charsize=2
end