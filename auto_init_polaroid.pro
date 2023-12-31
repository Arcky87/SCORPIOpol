function flat_norm,flat_cube
; files - list of flat files
; bias - level
bias=1000.0

s=size(flat_cube) & N=s(4) & Nx=s(1) & Ny=s(2)
cube=fltarr(Nx,Ny,3,N)
normal=fltarr(3,N)

for k=0,N-1 do begin
	for i=0,2 do begin
		cube(*,*,i,k)=flat_cube(*,*,i,k)
		robomean,cube(Nx/6.0:5.0*Nx/6.0,Ny/6.0:5.0*Ny/6.0,i,k),3,0.5,avgg,rms
		;robomean,cube(Nx/2-Nx/12:Nx/2+Nx/12,Ny/2-Ny/12:Ny/2+Ny/12,i,k),3,0.5,avgg,rms
		normal(i,k)=avgg ;& print, avgg
	endfor
endfor

mea=fltarr(4)
for i=0,2 do mea(i)=total(normal(i,*),2)/N
for i=0,2 do begin
	for k=0,N-1 do begin
		normal(i,k)=normal(i,k)/mea(i)
	endfor
endfor

ima=fltarr(Nx,Ny,4)
for i=0,2 do begin
	for kx=0,Nx-1 do begin
		for ky=0,Ny-1 do begin
			ima(kx,ky,i)=median(cube(kx,ky,i,*)/normal(i,*))
		endfor
	endfor
endfor
	for i=0,2 do begin
		robomean,ima(Nx/2-Nx/12:Nx/2+Nx/12,Ny/2-Ny/12:Ny/2+Ny/12,i),3,0.5,avgg,rms
		ima(*,*,i)=ima(*,*,i)/avgg
	endfor

;delete zeros
for i=0,2 do begin
	for kx=0,Nx-1 do begin
		for ky=0,Ny-1 do begin
			if ima(kx,ky,i) lt 0.85 then ima(kx,ky,i)=0.99
		endfor
	endfor
endfor

return, ima
end


;initial files for photometry

;PARAMETERS
rdir = '/home/elias/SCORPIO/sppol_pipeline_v2023.8/s191125'
wdir = '/home/elias/SCORPIO/sppol_pipeline_v2023.8/reduced'
fldir = '/home/elias/SCORPIO/sppol_pipeline_v2023.8/s191125'
cube_obj = ['s149508']
cube_flat = 's149508'
filters = ['V']
limexp=0.0
checkdark=0 ;1 or 0
;=================================
Nx=1024 & Ny=1024
bias=1000.0

;----------------------------------

;test wdir & oldies
if FILE_TEST(wdir) eq 0 then SPAWN,'mkdir '+wdir
if FILE_TEST(wdir+'cube_i.fts') eq 1 then FILE_DELETE,FILE_SEARCH(wdir+'cube_i.fts')
if FILE_TEST(wdir+'flat_i.fts') eq 1 then FILE_DELETE,FILE_SEARCH(wdir+'flat_i.fts')
if FILE_TEST(wdir+'flat.fts') eq 1 then FILE_DELETE,FILE_SEARCH(wdir+'flat.fts')
if FILE_TEST(wdir+'cube.fts') eq 1 then FILE_DELETE,FILE_SEARCH(wdir+'cube.fts')

;copy flat cube and unzip
FILE_COPY,fldir+cube_flat+'.zip',wdir,/OVERWRITE, /ALLOW_SAME
unzip,wdir,cube_flat

;find true flats
LIST=FILE_SEARCH(wdir+'s*.fts')
Nlist=N_elements(LIST)
fl_n='' & fl_p='' & fl_m=''
for k=0,Nlist-1 do begin
	header=headfits(list(k))
	type=strmid(SXPAR(HEADER,'MODE'),0,6)
	CASE type OF
	  'IMA':BEGIN
	  		FILE_DELETE,(list(k))
	  		end
	  'UNK':BEGIN
	  		FILE_DELETE,(list(k))
	  		end
	  'POL':BEGIN
	  		FILE_DELETE,(list(k))
	  		end
	  'ImaPol':BEGIN
;	 		file_fla=[file_fla,sxpar(header,'FILE')]
	 			a2=strmid(sxpar(header,'POLAMODE'),9,3) & print, a2
				f2=strmid(sxpar(header,'FILTERS'),2,1) & print, f2
				ty=strmid(sxpar(header,'IMAGETYP'),0,3) & print, ty
				if (a2 eq '-60' and f2 eq filters and ty eq 'neo') then fl_m=[fl_m,sxpar(header,'FILE')]
				if (a2 eq '+60' and f2 eq filters and ty eq 'neo') then fl_p=[fl_p,sxpar(header,'FILE')]
				if (a2 eq '0' and f2 eq filters and ty eq 'neo') then fl_n=[fl_n,sxpar(header,'FILE')]
	 	 	END
	 ENDCASE
ENDFOR

;check line and continuum filters
Nf=N_elements(fl_m) ;suppose the arrays are equal!
file_flat=strarr(3,Nf-1)
file_flat(0,*)=fl_m(1:Nf-1)
file_flat(1,*)=fl_n(1:Nf-1)
file_flat(2,*)=fl_p(1:Nf-1)

;extracted flat
flat=fltarr(Nx,Ny,3,Nf-1)
for kp=0,2 do begin
	for ke=0,Nf-2 do begin
		tmp=readfits(wdir+file_flat(kp,ke),h)
		flat(*,*,kp,ke)=tmp - bias
	endfor
endfor
writefits,wdir+'flat_i.fts',flat,h

;norm flat
flat_n=flat_norm(flat)
writefits,wdir+'flat.fts',flat_n,h

FILE_DELETE,FILE_SEARCH(wdir+'s*.fts')

;+++++++++++++++++++
;OBJECT

;copy cube and unzip
FILE_COPY,rdir+cube_obj+'.zip',wdir,/OVERWRITE, /ALLOW_SAME
ntmp=n_elements(cube_obj)
for j=0,ntmp-1 do unzip,wdir,cube_obj(j)

;find true frames
LIST=FILE_SEARCH(wdir+'s*.fts')
Nlist=N_elements(LIST)
ob_n='' & ob_p='' & ob_m=''
for k=0,Nlist-1 do begin
	header=headfits(list(k))
	type=strmid(SXPAR(HEADER,'MODE'),0,6)
	CASE type OF
	  'IMA':BEGIN
	  		FILE_DELETE,(list(k))
	  		end
	  'UNK':BEGIN
	  		FILE_DELETE,(list(k))
	  		end
	  'POL':BEGIN
	  		FILE_DELETE,(list(k))
	  		end
	  'ImaPol':BEGIN
;	 		file_fla=[file_fla,sxpar(header,'FILE')]
	 			a2=strmid(sxpar(header,'POLAMODE'),9,3) & print, a2
				f2=strmid(sxpar(header,'FILTERS'),2,1) & print, f2
				ty=strmid(sxpar(header,'IMAGETYP'),0,3) & print, ty
				if (a2 eq '-60' and f2 eq filters and ty eq 'obj') then ob_m=[ob_m,sxpar(header,'FILE')]
				if (a2 eq '+60' and f2 eq filters and ty eq 'obj') then ob_p=[ob_p,sxpar(header,'FILE')]
				if (a2 eq '0' and f2 eq filters and ty eq 'obj')   then ob_n=[ob_n,sxpar(header,'FILE')]
	 	 	END
	 ENDCASE
ENDFOR

;check line and continuum filters
No=N_elements(ob_m)
file_obj=strarr(3,No-1)
file_obj(0,*)=ob_m(1:No-1)
file_obj(1,*)=ob_n(1:No-1)
file_obj(2,*)=ob_p(1:No-1)

;extracted flat
obj=fltarr(Nx,Ny,3,No-1)
for kp=0,2 do begin
	for ke=0,No-2 do begin
		tmp=readfits(wdir+file_obj(kp,ke),h)
		obj(*,*,kp,ke)=tmp - bias
	endfor
endfor
sxaddpar,h,'DATAMIN',string(-65535.0)
writefits,wdir+'cube_i.fts',obj,h

FILE_DELETE,FILE_SEARCH(wdir+'s*.fts')
FILE_DELETE,FILE_SEARCH(wdir+'*.zip')

;FINAL
obj=readfits(wdir+'cube_i.fts')
flat=readfits(wdir+'flat.fts')
obj_cor=obj

for j=0,No-2 do begin
	for p=0,2 do begin
		obj_cor(*,*,p,j)=obj(*,*,p,j)/flat(*,*,p)
	endfor
endfor

;delete background
for j=0,No-2 do begin
	for p=0,2 do begin
		robomean,obj_cor(Nx/6.0:5.0*Nx/6.0,Ny/6.0:5.0*Ny/6.0,p,j),3,0.5,avgg,rms
		obj_cor(*,*,p,j)=obj_cor(*,*,p,j) - avgg
	endfor
endfor

hist=''
target=['ANGLE = -60', 'ANGLE = 0', 'ANGLE = +60']
for kr=0,2 do begin
	hist=[hist,target(kr)]
		for ke=0,No-2 do begin
			hist=[hist,'    '+file_obj(kr,ke)]
		endfor
endfor
sxaddhist,hist,h

writefits,wdir+'cube.fts',obj_cor,h


END
