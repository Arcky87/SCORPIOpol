function read_cube,wdir,cubename,filename,head
;+
; NAME:
;	READ_FILE
; PURPOSE:
;	read files in zipped  SCORPIO cubes.
; DESCRIPTION:
;
; CALLING SEQUENCE:
;	Result= READ_CUBE ( wdir, cubename, filename, head)
;
; CATEGORY:
;	reduction SCORPIO-data
;
; INPUTS:
;	rdir - string name input directory for reading
;	cubename - string name cube with extintion '.zip'
;	filename - string name file with extention '.fts'
; OUTPUTS:
;	Result = 2D ploat point image with FITS-header
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
;       Written by Victor Afanasiev, Special Astrophysical Observatory RAS, mar, 2014

;SPAWN,'C:\Users\username\7-Zip\7z.exe x '+wdir+' '+cubename+' '+filename,/hide
;SPAWN, wdir+'7z.exe x '+wdir+' '+cubename+' '+filename,/hide
SPAWN, '7zzs x '+wdir+cubename+' '+filename + ' -o'+wdir;,/hide

RES=READFITS(wdir+filename,head);,/silent)

file_delete,wdir+filename
return,RES
end
wdir='D:\SCORPIO\red_data.pol\NGC4151_160307\'
cubename='s138711.zip'
filename='s13871124.fts'
ima=read_cube(wdir,cubename,filename,head)
print,size(ima)
end
