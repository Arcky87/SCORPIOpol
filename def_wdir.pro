function DEF_WDIR,FILELOG
;+
; NAME:
;	DEF_WDIR
; PURPOSE:
;	DEFINIRION DIRECTORY FOR WRITING DATA 
; DESCRIPTION: 
;	
; CALLING SEQUENCE:
;	Result=DEF_WDIR(FILELOG)
;
; CATEGORY:
;	reduction MPFS-data		
;
; INPUTS:
;	LOGFILE = file name of LOG observation (in format FITS-header)
;	
;
; OUTPUTS:
;	Result = string scalar name directory for writing
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
RESULT=sxpar(TABLE,'W_DIR')
RETURN,RESULT
END
