dir='C:\RED_DATA\SCORPIO\B0624+6907_V_170130\'
cu=readfits(dir+'cube.fts')

a=size(cu) & Nexp=a(4)
!p.multi=[0,1,1]

coords=find_objects(cu(*,*,0,1),0.7,0.7)
print, coords
b=size(coords)
if (b(0) eq 1) then begin
	Nobj=1
endif else begin
	Nobj=b(2)
endelse

w=15
fwhm=3
skyap=[10,15]
aps=fltarr(3,Nobj,3,Nexp)

for ex=0,Nexp-1 do begin
	for ob=0,Nobj-1 do begin
		for pol=0,2 do begin
			ima=cu(*,*,pol,ex)
			optimal_aperture_short, ima, coords(*,ob), w, fwhm, sky, aper_data
			aps(*,ob,pol,ex)=aper_data
		endfor
	endfor
endfor

aper_size=fltarr(Nobj)
for ob=0,Nobj-1 do begin
	aper_size(ob)=fwhm*median(aps(2,ob,*,*)) & print, 'Aperture size for object', string(ob), aper_size(ob)
endfor

if Nobj ge 2 then begin
	robomean,aper_size,1,0.5,ap,rms
endif else begin
	ap=aper_size
endelse

ap=15
print, 'Final aperture: ', ap

flux=fltarr(Nobj,4,Nexp) & fluxerr=fltarr(Nobj,4,Nexp)
for ex=0,Nexp-1 do begin
	for ob=0,Nobj-1 do begin
		for pol=0,2 do begin
			ima=cu(*,*,pol,ex)
			aper, ima, aps(0,ob,pol,ex), aps(1,ob,pol,ex), fl, err_fl, sky, err_sky, 1, ap, $
				[ap+5,sqrt(4*ap*ap+10*ap+25)], [-32767,60000], /flux, /silent
			flux(ob,pol,ex)=fl ;- sky
			fluxerr(ob,pol,ex)=err_fl
		endfor
	endfor
endfor

for ob=0,Nobj-1 do begin
	file=dir+'F_'+string(ob,format='(I2.2)')+'.txt'
	print, file
	openw, ob+4, file
		for k=0,Nexp-1 do begin
				printf, ob+4, k, flux(ob,0,k),flux(ob,1,k),flux(ob,2,k)
		endfor
	close, ob+4
endfor

;coordinates
	file=dir+'coordinates.txt'
	print, file
	openw, 4, file
		for ob=0,Nobj-1 do begin
			printf, 4, ob, coords(0,ob), coords(1,ob)
		endfor
	close, 4

end