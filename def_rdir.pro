function DEF_RDIR,FILELOG
;+
; NAME:
;	DEF_RDIR
; PURPOSE:
;	DEFINIRION DIRECTORY FOR READING DATA FROM LOG FILE
; DESCRIPTION: 
;	
; CALLING SEQUENCE:
;	Result=DEF_RDIR(FILELOG)
;
; CATEGORY:
;	reduction MPFS-data		
;
; INPUTS:
;	LOGFILE = file name of LOG observation (in format FITS-header)
;	
;
; OUTPUTS:
;	Result = string scalar name directory for reading 
;		
; OPTIONAL OUTPUT:
;	no
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
;	SXPAR
;
; MODIFICATION HISTORY:
;       Written by Victor Afanasiev, Special Astrophysical Observatory RAS, Jul 1999
;-
on_error,2
openr,UNIT,FILELOG,/GET_LUN
a=' '
readf,UNIT,a
TABLE=a
WHILE NOT EOF(UNIT) do begin
readf,UNIT,a
TABLE=[TABLE,a]
ENDWHILE
close,UNIT
FREE_LUN, UNIT
RESULT=sxpar(TABLE,'R_DIR')
RETURN,RESULT
END
