pro create_initial_data_WOLL2,LOGFILE


;анализ данных по кубам
name=['OBJ','star_0','star']
;create work directory.
wdir=def_wdir(LOGFILE)

if FILE_TEST(wdir) eq 0 then SPAWN,'mkdir '+wdir
;изничтожение старых  FITS-файлов
;tmp=FILE_SEARCH(wdir+'*.fts')
;if tmp(0) ne '' then  FILE_DELETE,FILE_SEARCH(wdir+'*.fts')
rdir=def_rdir(LOGFILE)
ext=def_ext(LOGFILE)
cube_filename=strarr(3)
file_obj=''  & file_nul=''  & file_sta=''
file_neo=''  & file_fla=''  & file_eta=''
for J=0,2 DO BEGIN
cube_filename(j)=sxpar(read_table(LOGFILE),name(j))
FILE_COPY,rdir+cube_filename(j)+'.zip',wdir,/OVERWRITE

cd,wdir,CURRENT=old_dir
spawn,'7z.exe x '+cube_filename[J]+ext,/hide
LIST=FILE_SEARCH(wdir+'*.fts')
Nlist=N_elements(list)
;селекция типов данных

;выделение режима наблюдений
;goto,cont

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
 	ENDFOR
 tmp=FILE_SEARCH(wdir+'*.fts')
if tmp(0) ne '' then  FILE_DELETE,FILE_SEARCH(wdir+'*.fts')
	ENDFOR

print,'obj : ',file_obj
print,'null: ',file_nul
print,'star: ',file_sta
print,'neon: ',file_neo
print,'flat: ',file_fla
print,'eta : ',file_eta
;определение пределов вырезания frames
curr_cube=strmid(file_eta(1),0,STRLEN(file_eta(1))-6)
H_frame=140
ima=READ_CUBE_WOLL2(wdir,curr_cube+'.zip',file_eta(1),h)
y_c=center_frames_WOLL2(ima)
if y_c(0) lt H_frame/2 then d_y=H_frame/2-Y_c(0)+5 ELSE d_y=0
;формирование куба данных объекта
HIST='INITIAL FILES'
cub=['object','unpolarized star','polarized star'] & Ncube=3
Nobj=N_elements(file_obj)-1  & file_obj=file_obj(1:Nobj)
Nnul=N_elements(file_nul)-1 & file_nul=file_nul(1:Nnul)
Nsta=N_elements(file_sta)-1 & file_sta=file_sta(1:Nsta)
ima=read_cube(wdir,cube_filename(0)+'.zip',file_obj(0),head)
Nx=sxpar(head,'NAXIS1')  & Ny=sxpar(head,'NAXIS2')

Ncube=3  & Nray=4
header=create_header_WOLL_2(Nx,Ny,Nray,Nobj,Ncube,0,1)
if Nnul gt Nobj then Nnul=Nobj
if Nsta gt Nobj then Nsta=Nobj
files=strarr(Nobj,Ncube)
files(0:Nobj-1,0)=file_obj
files(0:Nnul-1,1)=file_nul(0:Nnul-1)
files(0:Nsta-1,2)=file_sta(0:Nsta-1)
cube_obj=fltarr(Nx-20,H_frame,Nray,Nobj,Ncube)
Nexp=[Nobj,Nnul,Nsta]

TARGET=['************ TARGET: OBJECT ***********',$
        '** TARGET: UNPOLARIZED STANDARD STAR **',$
        '*** TARGET: POLARISED STANDARD STAR ***']
;
for k=0,Ncube-1 do begin
;модификация FITS-шапки
sxaddpar,header,'cube'+string(k+1,format='(I1)'),sxpar(read_table(LOGFILE),cub(k))
FOR J=0,Nexp(k)-1 DO BEGIN
tmp=FLOAT(READ_CUBE(wdir,cube_filename(k)+'.zip',files(J,k),h))
;маскирование плохой строки
tmp(*,47)=tmp(*,45) & tmp(*,46)=tmp(*,48)
if j eq 0 then hist=[hist,TARGET(k)]
hist=[hist,'        '+files(j,k)]
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
sxaddhist,hist,header
;for j=0,N_elements(header)-1 do print,header(j)
writefits,wdir+'obj_i.fts',cube_obj,header
;cube=fltarr(Nx,Ny,N_elements(file_fla)]
;формирование исходной записи FLAT
;for k=0,N_elements(file_fla)-1 do begin

name=['neon','flat','eta']
for k=0,2 do begin
hist='INITIAL FILES'
if k eq 0 then  file=file_neo
if k eq 1 then  file=file_fla
if k eq 2 then  file=file_eta
N=N_elements(file) & file=file(1:N-1)  & N=N-1
cube_out=fltarr(Nx-20,H_frame,Nray,N)

for j=0,N-1 do begin
curr_cube=strmid(file_fla(j),0,STRLEN(file_fla(j))-6)
tmp=FLOAT(READ_CUBE(wdir,curr_cube+'.zip',file(j),h))
;маскирование плохой строки
if k lt 2 then begin
tmp(*,47)=tmp(*,45)
tmp(*,46)=tmp(*,48)
endif
cube_out(*,*,*,j)=frames_WOLL2(tmp,YC=Y_c,DY=d_Y,H=H_frame,/bias)
HIST=[HIST,file(j)]
endfor
sxaddhist,hist,h
writefits,wdir+name(k)+'_i.fts',total(cube_out,4),h

ENDFOR
tmp=FILE_SEARCH(wdir+'*.zip')
if tmp(0) ne '' then  FILE_DELETE,FILE_SEARCH(wdir+'*.zip')
end
;LOGFILE='h:\red_data.pol\LOGS\3C390_110825.txt'
LOGFILE=DIALOG_PICKFILE(/read,path='h:\red_data.pol\LOGS\',FILTER='*.txt')
print,logfile
create_initial_data_WOLL2,LOGFILE
W_DIR=sxpar(read_table(LOGFILE),'w_dir')
LOADFILE,dir=w_dir,ysz=500,xsz=1130

ViewPol_2
end