function READ_FILE,rdir,filename,ext,head
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
;	Result = 2D ploat point array MPFS-image
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
head=0
if ext eq '.zip' or ext eq '.ZIP' then begin
file_zip=rdir+strmid(filename,0,strlen(filename)-2)+ext
file_fts=strupcase(filename)+'.FTS'
;if os_family() eq 'unix' then zip_pro='unzip ' else zip_pro='pkunzip '; for 32-byte XP
if os_family() eq 'unix' then zip_pro='unzip ' else zip_pro='7z.exe x '; for 64-byte XP
if os_family() eq 'unix' then $
SPAWN,zip_pro+file_zip+' '+file_fts else $
SPAWN,zip_pro+file_zip+' '+file_fts,/hide
result=READFITS(file_fts,head,SILENT='silent')
if os_family() eq 'unix' then del_pro='rm -f ' else del_pro='del '
if os_family() eq 'unix' then $
SPAWN,del_pro+file_fts else $
SPAWN,del_pro+file_fts,/hide
endif
if ext eq '.fts' then begin
result=READFITS(rdir+filename+'.fts',head,/silent)
endif
return,result
end
ima=READ_FILE( 'f:\obs_data\s100716\','s7750323','.zip',h)
print,h
end
