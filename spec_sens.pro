; SCORPIO longslit dqe determination

PRO spec_sens,table_name,filein=filein,$
                w_dir=w_dir,skip=skip

if n_params() eq 0 then begin

print,'spec_sens,table_name,filein=filein,w_dir=w_dir,skip=skip'
print,''
print,'Skip -- skip around broad lines  (default=50 A) and a smoothing willbe performed with window=3*skip'
print,'Set GAUS=1 for gaussian sum instead easy sum in the WIN'
return
endif

help,skip,win
device,dec=0 & loadct,27

if not(keyword_set(filein)) then filein='obj-sky.fts'
if keyword_set(w_dir) then w_dir=slash(w_dir) else w_dir=''
if not(keyword_set(skip)) then skip=50 ; skipping in A around lines
; read file
im=readfits(w_dir+filein,head)
 if n_elements(image) eq 1 then begin
  res=dialog_message('No file '+w_dir+filein(i))
  return
 endif

Xs=sxpar(head,'naxis1')

;Ys=sxpar(head,'naxis2')
;yy=findgen(Ys)
ld0 =sxpar(head,'CRVAL1')
dld0=sxpar(head,'CDELT1')
lambda=ld0+dld0*findgen(Xs); wavelength grid
ADU=sxpar(head,'GAIN')
z_star=float(sxpar(head,'Z'))
T_exp=sxpar(head,'EXPTIME2')

tot=fltarr(xs)
gau=fltarr(xs)

; read standard star table
;table=read_ascii(w_dir+table_name)
;print, (table)
;if n_elements(table) lt 3 then return

file=w_dir+table_name
N=numlines(file)
table=fltarr(3,N)
openr, 1, file
readf, 1, table
close, 1

;  object  integration
;
;; object search
;sum=total(im(xs*0.2:xs*0.8,*),1)
;m=max(sum,mp)
;for x=00,xs-1 do begin
; tot(x)=total(im(x,(mp-Win)>0:(mp+Win)<(ys-1)))
; ; gauss fitting along slit
; if keyword_set(gaus) then begin
;
;  g=multigaus( yy,im(x,*),mp,fwhm=win)
;  gau(x)=g(0).flux
; endif
;if x mod 300 eq 0 then print,'X:',x
;ENDFOR
;
;if keyword_set(gaus) then begin
;window,0
;   plot,lambda,tot,xst=1,xtit='wavelength',ytit='Flux',col=200
;   oplot,lambda,gau,col=20
;   legend,['Total','Gaussian'],col=[200,20],textcolor=[200,20]
;   tot=gau
;endif

tot=im


lambda_tab=reform(table(0,*))
mag_tab=reform(table(1,*))

s=size(table)
; set step
if s(1) ge 3 then w_tab=reform(table(2,*)) else begin ;!!!
  w_tab=lambda_tab-shift(lambda_tab,1)
  w_tab(0)=w_tab(1)

  nn=n_elements(w_tab)
  w_tab(nn-1)=w_tab(nn-2)
endelse


print,'Apertures in the table:',minmax(w_tab)
;calculation extintion
S=2.51E05       ;total square mirror of telescope in cm^2
a=0.012 & c=0.12
extin_tab=(a*1./((lambda_tab/10000.)^4.)+c)*2.5/2.3
mag_tab=mag_tab+extin_tab/cos(z_star*!DTOR)
N_tab=948.*S*(10^(-0.4*mag_tab))/ADU
N_tab=N_tab*(5500./lambda_tab)^2


vector=tot/(T_exp*dld0); counts/sec/A

; shift lambda to true wavelengh scale:
Wmed=Xs*0.2
temp=INTERPOL(mag_tab,lambda_tab,lambda)
obs=vector-med_ext(vector,Wmed)
temp=med_ext(temp,Wmed)-temp
apod=cosin_apod(Xs,20)
m=150 & x_cross=findgen(2*M+1)-M
cross=CROSS_NORM(obs*apod,temp*apod,m)
gau=multigaus(x_cross,cross,0,fwhm=m/10)

if gau(0).flux eq -1 or abs(gau(0).center) gt m/2 then sh=0 else sh=-gau(0).center

vector=vecshift(vector,dx=sh)
print,'Wavelength shift [px]: ',sh


; **********
; smoothed observed data:

rec=where (lambda_tab gt lambda(0) and lambda_tab lt lambda(Xs-1),num_t)
l_obs=lambda_tab(rec)
f_obs=fltarr(num_t)

w=median(w_tab(rec))

if w GT 2*dld0 then begin

for j=0,num_t-1 do begin
 l_c=l_obs(j)
 w=w_tab(rec(j))
 inde=where (lambda ge l_c-W/2 and lambda Le l_c+W/2,nums)
 if nums gt 0 then f_obs(j)=total(vector(inde))/nums

endfor
endif else f_obs=INTERPOL(vector,lambda,l_obs)

;rec=where(f_obs eq 0,ccc)

DQE_tab=100.*f_obs/N_tab(rec)

; DQE smoothing and interpolation

dw=skip/dld0
bad_lambda=[4600,4650,6250,6860]  ;[4340,4861,5420,6562,6860]
 Nbad=N_elements(bad_lambda)
for j=0,Nbad-1 do begin

index=where(abs(l_obs-bad_lambda(j)) gt dw)
l_obs=l_obs(index)
dqe_tab=dqe_tab(index)
endfor

; Smoothing DQE
dqe_smo=median(dqe_tab,3)
DQE_average=INTERPOL(dqe_smo,l_obs,lambda)

old=dqe_average
dqe_average=smooth(dqe_average,2*dw,/edge)
dqe_average(0:1.5*dw)=smooth(old(0:1.5*dw),0.2*dw,/edge)

window,5
plotsym,0,0.2,/fill
plot,lambda,DQE_average,xst=1,xtit='wavelength',ytit='DQE'
oplot,l_obs,DQE_tab,psym=8,col=20


;write sent
sent=3.39E-9*ADU/(9.48*DQE_average*S)
mkhdr,h_sent,sent
sxaddpar,h_sent,'CRVAL1', ld0
sxaddpar,h_sent,'CDELT1', dld0
sxaddpar,h_sent,'BUNIT', 'erg/cm^2/sec/A'
writefits,w_dir+'sent.fts',sent,h_sent


sxdelpar,h_sent,'BUNIT'
writefits,w_dir+'dqe.fts',dqe,h_sent

;FDECOMP, fileout, disk, dir, name
fileps='DQE.ps'
set_plot,'ps'
device,file=w_dir+fileps,xs=18,ys=16,xoff=1,yoff=5
loadct,0,/sil
!p.multi=[0,1,1]
plot,lambda,DQE_average,xst=1,xtit='wavelength',ytit='DQE',tit=sxpar(head,'OBJECT')+'   '+sxpar(head,'DATE')
oplot,l_obs,DQE_tab,psym=8


 if  !VERSION.OS_family eq 'Windows' then set_plot,'win' else set_plot,'X'

loadct,27,/sil
end