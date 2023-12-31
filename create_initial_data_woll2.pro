pro create_initial_data_WOLL2,LOGFILE,Y_SHIFT=Y_shift
Z=sxpar(read_table(LOGFILE),'REDSHIFT')
if not(keyword_set(Y_shift)) then Y_shift=0
;data analisys in cubes
name=['OBJ','star_0','star']
;create work directory
wdir=def_wdir(LOGFILE) & rdir=def_rdir(LOGFILE) & ext=def_ext(LOGFILE)
if FILE_TEST(wdir) eq 0 then SPAWN,'mkdir '+wdir
;give x limits
;Xstart=20 & Xfinish=1700	;20 - overscan  ;1026 + GS17 ;4290
;Xstart=20 & Xfinish=1400	;20 - overscan  ;940 + LP425 ;261  ??
;Xstart=25 & Xfinish=1620	;20 - overscan  ;940         ;42-90
;Xstart=250 & Xfinish=1900	;20 - overscan  ;940         ;42-90
;Xstart=420 & Xfinish=1620	;20 - overscan  ;940         ;42-90  ;quartz
;Xstart=390 & Xfinish=2150	;20 - overscan  ;940 + GS11  ;42-90
;Xstart=360 & Xfinish=2200	;20 - overscan  ;940 + GS11  ;42-90
;Xstart=250 & Xfinish=2300	;20 - overscan  ;1200 + GS11 ;42-90
;Xstart=24 & Xfinish=2300	;20 - overscan  ;1026 + GS17 ;42-90
;Xstart=90 & Xfinish=2150	;20 - overscan  ;1200        ;42-90
;Xstart=100 & Xfinish=1700	;20 - overscan  ;940         ;42-90
;Xstart=20 & Xfinish=2320	;20 - overscan  ;1800        ;42-90
;Xstart=10 & Xfinish=2055	;20 - overscan  ;1026 + LP525 ;261
;Xstart=50 & Xfinish=1850	;20 - overscan  ;1026 + GS17 ;42-90
;Xstart=24 & Xfinish=4100	;20 - overscan  ;1026 + GS17 ;261    ;1x2
;Xstart=24 & Xfinish=2050	;20 - overscan  ;1026 + GS17 ;261    ;2x4
;Xstart=310 & Xfinish=1900	;20 - overscan  ;940 + LP425 ;261    ;2x4
;Xstart=630 & Xfinish=3800	;20 - overscan  ;940 + LP425 ;261    ;1x2
Xstart=350 & Xfinish=2050	;20 - overscan  ;940 + LP425 ;261    ;1x2


;====================================================================
;'OBJ','star_0','star'
cube_filename=strarr(3)
file_obj=''  & file_nul=''  & file_sta=''
file_neo=''  & file_fla=''  & file_eta=''
;choose file names
for J=0,2 DO BEGIN
	print,'INIT DATA for  ', name(j)
	cube_filename(j)=sxpar(read_table(LOGFILE),name(j))
	cd,wdir,CURRENT=old_dir
	FILE_COPY,rdir+cube_filename(j)+'.zip',wdir,/OVERWRITE
	spawn,'/home/elias/bin/7zzs x '+cube_filename[J]+ext + ' -o'+wdir
	;data type selection
	LIST=FILE_SEARCH(wdir+'*.fts')
	Nlist=N_elements(list)
		for k=0,Nlist-1 do begin
			header=headfits(list(k))
			type=strmid(SXPAR(HEADER,'IMAGETYP'),0,3)
			mode=strmid(SXPAR(HEADER,'MODE'),0,3)
			mask=strmid(SXPAR(HEADER,'SLITMASK'),0,3)
			woll=strmid(SXPAR(HEADER,'FILTERS'),2,6)

			if type eq 'fla' and mode eq 'Spe' then  type='fla'
			if type eq 'fla' and mode eq 'Spe' and mask eq '3 d' then type='eta'

			CASE type OF
			  'bia':BEGIN
			  		end
			  'obj' :BEGIN
			 		if mode eq 'Spe' then begin
			 		if J eq 0  then file_obj=[file_obj,sxpar(header,'FILE')]
			 		if J eq 1  then file_nul=[file_nul,sxpar(header,'FILE')]
			 		if J eq 2  then file_sta=[file_sta,sxpar(header,'FILE')]
			 		endif
			 		END
			 'neo':BEGIN
			 		file_neo=[file_neo,sxpar(header,'FILE')]
			 		END
			 'fla':BEGIN
			 		file_fla=[file_fla,sxpar(header,'FILE')]
			 	 	END
			 'eta':BEGIN
			  		file_eta=[file_eta,sxpar(header,'FILE')]
			  	 	END
			 ENDCASE
	 	endfor
	tmp=FILE_SEARCH(wdir+'*.fts')
		if tmp(0) ne '' then  FILE_DELETE,FILE_SEARCH(wdir+'*.fts')
endfor

;list of name for every type
print,'obj : ',file_obj
print,'null: ',file_nul
print,'star: ',file_sta
print,'neon: ',file_neo
print,'flat: ',file_fla
print,'eta : ',file_eta

;determination of limits of frames
;1) H_frame - y-size of the extracted frames
curr_cube=strmid(file_eta(1),0,STRLEN(file_eta(1))-6)  & print,file_eta(1)
H_frame=125 ; 120 --- default IY
ima=READ_CUBE_WOLL2(wdir,curr_cube+'.zip',file_eta(1),h)
bin=sxpar(h,'BINNING') & bin=str_sep(bin,'x') & bin=FIX(bin(1))
H_frame=H_frame*4.0/bin
Ny=H_frame

ima=shift(ima,0,Y_shift)
y_c=CENTER_FRAMES_WOLL2(ima)
if y_c(0) lt H_frame/2 then d_y=H_frame/2-Y_c(0)+5 ELSE d_y=0

;2) Nx - x-size of the frames within the given limits
Nx=Xfinish - Xstart + 1

;create obj data cube
HIST='INITIAL FILES'
cub=['object','unpolarized star','polarized star'] & Ncube=3
Nobj=N_elements(file_obj)-1  & file_obj=file_obj(1:Nobj)
Nnul=N_elements(file_nul)-1 & file_nul=file_nul(1:Nnul)
Nsta=N_elements(file_sta)-1 & file_sta=file_sta(1:Nsta)
ima=read_cube(wdir,cube_filename(0)+'.zip',file_obj(0),head)

Nxi=sxpar(head,'NAXIS1')  & Nyi=sxpar(head,'NAXIS2')
Ncube=3  & Nray=4

;create header
header=create_header_WOLL_2(Nx,Ny,Nray,Nobj,Ncube,0,1)
if Nnul gt Nobj then Nnul=Nobj
if Nsta gt Nobj then Nsta=Nobj
files=strarr(Nobj,Ncube)
files(0:Nobj-1,0)=file_obj
files(0:Nnul-1,1)=file_nul(0:Nnul-1)
files(0:Nsta-1,2)=file_sta(0:Nsta-1)
cube_obj=fltarr(Nx,H_frame,Nray,Nobj,Ncube)
Nexp=[Nobj,Nnul,Nsta]

TARGET=['************ TARGET: OBJECT ***********',$
        '** TARGET: UNPOLARIZED STANDARD STAR **',$
        '*** TARGET: POLARISED STANDARD STAR ***']

for k=0,Ncube-1 do begin
	;modification of the FITS-header
	sxaddpar,header,'cube'+string(k+1,format='(I1)'),sxpar(read_table(LOGFILE),cub(k))
	FOR J=0,Nexp(k)-1 DO BEGIN
		tmp=FLOAT(READ_CUBE(wdir,cube_filename(k)+'.zip',files(J,k),h))
		tmp=tmp(Xstart:Xfinish,*)
			;������������ ������ ������
;			tmp(*,68)=tmp(*,69)  ; ��� 2x4
;			tmp(*,67)=tmp(*,65)
;			tmp(*,66)=tmp(*,65)
;						  ; ��� 2x1
;			tmp(*,190)=tmp(*,191)
;			tmp(*,189)=tmp(*,191)
;			tmp(*,188)=tmp(*,191)
;			tmp(*,187)=tmp(*,184)
;			tmp(*,186)=tmp(*,184)
;			tmp(*,185)=tmp(*,184)
			if j eq 0 then hist=[hist,TARGET(k)]
				hist=[hist,'        '+files(j,k)]
				print, 'size', size(tmp)
				print, 'size', size(cube_obj(*,*,*,j,k))
				print, 'size', Nx
				print, Y_c
				cube_obj(*,*,*,j,k)=frames_WOLL2(tmp,YC=Y_c,DY=d_Y,H=H_frame,/bias)
					if J eq 0 then begin
						PA=sxpar(h,'PARANGLE')-sxpar(h,'ROTANGLE')+132.5
						sxaddpar,header,'NAME'+string(k+1,format='(I1)'),sxpar(h,'OBJECT')
						sxaddpar,header,'START'+string(k+1,format='(I1)'),sxpar(h,'START')
						sxaddpar,header,'NUMEXP'+string(k+1,format='(I1)'),Nexp(k)
						sxaddpar,header,'EXPTIME'+string(k+1,format='(I1)'),sxpar(h,'EXPTIME')
						sxaddpar,header,'RA'+string(k+1,format='(I1)'),sxpar(h,'RA')
						sxaddpar,header,'DEC'+string(k+1,format='(I1)'),sxpar(h,'DEC')
						sxaddpar,header,'A'+string(k+1,format='(I1)'),sxpar(h,'A')
						sxaddpar,header,'Z'+string(k+1,format='(I1)'),sxpar(h,'Z')
						sxaddpar,header,'PA'+string(k+1,format='(I1)'),PA
						sxaddpar,header,'CUBE'+string(k+1,format='(I1)'),cub(k)
					endif
	ENDFOR
		if k eq 0 then begin
			sxaddpar,header,'NAXIS',5
			sxaddpar,header,'DATE-OBS',sxpar(h,'DATE')
			sxaddpar,header,'PROG-ID',sxpar(h,'PROG-ID')
			sxaddpar,header,'AUTHOR',sxpar(h,'AUTHOR')
			sxaddpar,header,'OBSERVER',sxpar(h,'OBSERVER')
			sxaddpar,header,'DIR',rdir
			sxaddpar,header,'FOCUS',sxpar(h,'FOCUS')
			sxaddpar,header,'BINNING',sxpar(h,'BINNING')
			sxaddpar,header,'RATE',sxpar(h,'RATE')
			sxaddpar,header,'GAIN',sxpar(h,'GAIN')
			sxaddpar,header,'NODE',sxpar(h,'NODE')
			sxaddpar,header,'IMSCALE',sxpar(h,'IMSCALE')
			sxaddpar,header,'CAMFOCUS',sxpar(h,'CAMFOCUS')
			sxaddpar,header,'COLFOCUS',sxpar(h,'COLFOCUS')
			sxaddpar,header,'MODE',sxpar(h,'MODE')
			sxaddpar,header,'DISPERSE',sxpar(h,'DISPERSE')
			sxaddpar,header,'SLITWID',sxpar(h,'SLITWID')
			sxaddpar,header,'SLITMASK',sxpar(h,'SLITMASK')
			sxaddpar,header,'FILTERS',sxpar(h,'FILTERS')
			sxaddpar,header,'FILTPOS1',sxpar(h,'FILTPOS1')
			sxaddpar,header,'FILTPOS2',sxpar(h,'FILTPOS2')
		ENDIF
endfor
;������ �������� �������� � obj_i.fts
sxaddpar,header,'Z',Z
sxaddhist,hist,header

writefits,wdir+'obj_i.fts',cube_obj,header

;====================================================================
;'neon','flat','eta'
name=['neon','flat','eta']
for k=0,2 do begin
	hist='INITIAL FILES'
	if k eq 0 then  file=file_neo
	if k eq 1 then  file=file_fla
	if k eq 2 then  file=file_eta
	N=N_elements(file) & file=file(1:N-1)  & N=N-1
	;cube_out=fltarr(Nx-20,H_frame,Nray,N)
	;cube_out=fltarr(Nx-350-149,H_frame,Nray,N)
	cube_out=fltarr(Nx,H_frame,Nray,N)
		for j=0,N-1 do begin
			curr_cube=strmid(file(j),0,STRLEN(file(j))-6)
			print, file(j),STRLEN(file(j))-6
			tmp=FLOAT(READ_CUBE(wdir,curr_cube+'.zip',file(j),h))
			tmp=tmp(Xstart:Xfinish,*)
			;������������ ������ ������
;			tmp(*,68)=tmp(*,69)  ; ��� 2x4
;			tmp(*,67)=tmp(*,65)
;			tmp(*,66)=tmp(*,65)
;						  ; ��� 2x1
;			tmp(*,190)=tmp(*,191)
;			tmp(*,189)=tmp(*,191)
;			tmp(*,188)=tmp(*,191)
;			tmp(*,187)=tmp(*,184)
;			tmp(*,186)=tmp(*,184)
;			tmp(*,185)=tmp(*,184)

			tmp=shift(tmp,0,Y_shift)
			cube_out(*,*,*,j)=frames_WOLL2(tmp,YC=Y_c,DY=d_Y,H=H_frame,/bias)
			HIST=[HIST,file(j)]
		endfor
	sxaddhist,hist,h
	si=size(cube_out)
	if si(0) eq 4 then writefits,wdir+name(k)+'_i.fts',total(cube_out,4),h
		if si(0) eq 3 then writefits,wdir+name(k)+'_i.fts',cube_out,h

ENDFOR
tmp=FILE_SEARCH(wdir+'*.zip')
if tmp(0) ne '' then  FILE_DELETE,FILE_SEARCH(wdir+'*.zip')
END


LOGFILE=DIALOG_PICKFILE(/read,path='/data6/SCORPIO/sppol_pipeline_v2023.8/LOGS/',FILTER='*')
print,logfile
create_initial_data_WOLL2,LOGFILE,Y_shift=10
W_DIR=sxpar(read_table(LOGFILE),'w_dir')
LOADFILE,dir=w_dir,ysz=500,xsz=1130

ViewPol_2
end
