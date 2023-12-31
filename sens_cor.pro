; SCORPIO long-slit data reduction
; correction to true flux
pro sens_cor,filein,fileout,w_dir=w_dir

if keyword_set(w_dir) then w_dir=slash(w_dir) else w_dir=''
sent=READFITS(w_dir+'sent.fts',h)
spectr=READFITS(w_dir+filein,h)

N_x=sxpar(h,'NAXIS1')	;number spectral element
;N_y=sxpar(h,'NAXIS2')    ;number of spectra
obj_red=fltarr(N_x)
param=[sxpar(h,'CRVAL1'),sxpar(h,'CDELT1'),N_x]


exp_obj=sxpar(h,'EXPTIME1')
z_obj=float(sxpar(h,'Z'))
object_name=STRCOMPRESS(sxpar(h,'OBJECT'))
date_obs=STRCOMPRESS(sxpar(h,'DATE-OBS'))
;calculation mean value extintion
a=0.01
c=0.14
lambda=findgen(N_x)*param(1)+param(0)
secz=1/cos(z_obj*!DTOR)
ext=(a*1/((lambda/10000.)^4)+c)*2.5/2.3
ext_obj=10.^(0.4*ext*secz)
 print,'Z=',z_obj, 'Ext_obj=',median(ext_obj)

;definition  unit
sent=sent/(exp_obj*param(1))
order=fix(min(ALOG10(sent)))
if order lt 0 then order=order-1
out_units=10.^order		;output units in erg/cm^2/sec/A
;obj_red=fltarr(N_x)

 obj_red=spectr*sent*ext_obj


sxaddpar,h,'BUNIT','erg/cm^2/sec/A',befor='COMMENT'
writefits,w_dir+fileout,obj_red,h

END
