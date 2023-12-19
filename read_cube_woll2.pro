function read_cube_WOLL2,wdir,cubename,filename,head
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

;SPAWN,'7z.exe x '+wdir+' '+cubename+' '+filename,/hide
SPAWN, '/home/elias/bin/7zzs x '+wdir+cubename+' '+filename + ' -o'+wdir;,/hide

RES=READFITS(wdir+filename,head,/silent)

file_delete,wdir+filename
return,RES
end
wdir='h:\red_data.pol\Arp102b_131103\'
cubename='S114201.zip'
filename='s11420108.fts'
ima=read_cube(wdir,cubename,filename,head)
print,size(ima), 'this?'
end
