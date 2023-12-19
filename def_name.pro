function DEF_NAME,FILELOG,TYPE_EXP,N_EXP
;+
; NAME:
;	DEF_NAME
; PURPOSE:
;	definition filename of different type exposure MPFS-data
; DESCRIPTION: 
;	
; CALLING SEQUENCE:
;	Result=DEF_NAME(FILELOG,TYPE_EXP,N_EXP)
;
; CATEGORY:
;	reduction MPFS-data		
;
; INPUTS:
;	LOGFILE = file name of LOG observation (in format FITS-header)
;	TYPE_EXP = string scalar type exposure (values: 'bias','obj','star','flat',
;		   'eta','star','test')
;
; OUTPUTS:
;	Result = string scalar filename for reading (without extention)
;		
; OPTIONAL OUTPUT:
;	N_exp - number exposures with every value TYPE_EXP
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
FILE=sxpar(TABLE,TYPE_EXP)
FILE=(STR_SEP(FILE,','))
NIGTH=SXPAR(TABLE,'NIGTH')
FILENAME=NIGTH+FILE
N_EXP=N_elements(FILENAME)
IF N_exp EQ 1 then begin
sw=byte(file)
if sw(0) eq 32 then N_exp=0
endif
RETURN,FILENAME
END
