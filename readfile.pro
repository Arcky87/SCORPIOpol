function READFILE,rdir,filename,ext,head,OVERSCAN=overscan,SCALE_X=scale_x,SCALE_Y=scale_y
;+
; NAME:
;	READFILE
; PURPOSE:
;	read files with SCORPIO-images. Will also read zip compressed FTS-files
; DESCRIPTION:
;
; CALLING SEQUENCE:
;	Result= READ_DFILE ( rdir, filename, ext, head)
;
; CATEGORY:
;	reduction MPFS-data
;
; INPUTS:
;	rdir - string name input directory for reading
;	filename - string name file withot extention
;	ext -	string name extention file for reading (value '.zip' or '.fts')
;		if extention e.q. '.zip' routine unzippiped files
; OUTPUTS:
;	Result = 2D ploat point array SCORPIO-image
;
; OPTIONAL OUTPUT:
;	head - string array contain FITS-header image
;
; OPTIONAL INPUT KEYWORDS:
;	no
;
; RESTRICTIONS:
;	no
;
; NOTES:
;	no
;
; PROCEDURES USED:
;	READ_FTS, SPAWN
;
; MODIFICATION HISTORY:
;       Written by Victor Afanasiev, Special Astrophysical Observatory RAS, Jul 1999
;-
if not(keyword_set(overscan)) then overscan=0
if not(keyword_set(scale)) then scale=1

head=0
if ext eq '.zip' or ext eq '.ZIP' then begin
file_zip=rdir+strmid(filename,0,strlen(filename)-2)+ext
file_fts=strupcase(filename)+'.FTS'
;if os_family() eq 'unix' then zip_pro='unzip ' else zip_pro='pkunzip '; for 32-byte XP
if os_family() eq 'unix' then zip_pro='unzip ' else zip_pro='7z.exe x '; for 64-byte XP
if os_family() eq 'unix' then $
SPAWN,zip_pro+file_zip+' '+file_fts else $
SPAWN,zip_pro+file_zip+' '+file_fts,/hide
result=READFITS(file_fts,head);,SILENT='silent')
if os_family() eq 'unix' then del_pro='rm -f ' else del_pro='del '
if os_family() eq 'unix' then $
SPAWN,del_pro+file_fts else $
SPAWN,del_pro+file_fts,/hide
endif
if ext eq '.fts' then begin
result=READFITS(rdir+filename+'.fts',head);,/silent)
endif
;удаление overscan
a=size(result)
result=result(overscan:a(1)-1,0:a(2)-overscan-1)
a=size(result)
result=congrid(result,a(1)*SCALE_X,a(2)*SCALE_Y)
a=size(result)
sxaddpar,head,'NAXIS1',a(1)
sxaddpar,head,'NAXIS2',a(2)
;scale=0.5
;overscan=20
;a=size(result)
;Nx=a(1)-overscan
;Ny=a(2)-overscan
;map=intarr(Nx*scale+overscan,Ny*scale+scale)
;map(overscan:Nx*scale+overscan-1,0:Ny*scale-1)=congrid(result(overscan:Nx+overscan-1,0:Ny-1),Nx*scale,Ny*scale)


;print,size(map)
return,result
end
ima=READFILE( 'h:\obs_data.pol\s111117\','s9340110','.zip',h,overscan=20,scale_x=0.5,scale_y=0.5)

end
