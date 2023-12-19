function DEF_ext,FILELOG
;+
; NAME:
;	DEF_NAME
; PURPOSE:
;	definition file extention for SCORPIO-data
; DESCRIPTION:
;
; CALLING SEQUENCE:
;	Result=DEF_EXT(FILELOG)
;
; CATEGORY:
;	reduction SCORPIO-data
;
; INPUTS:
;	LOGFILE = file name of LOG observation (in format FITS-header)
;
; OUTPUTS:
;	Result = string scalar
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
;       Written by Victor Afanasiev, Special Astrophysical Observatory RAS, Oct 2001
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
EXTENTION=SXPAR(TABLE,'EXT')
RETURN,EXTENTION
END
