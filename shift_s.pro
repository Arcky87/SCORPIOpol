function shift_s,vector,dx
;+
; NAME:        
;	 SHIFT_S
; PURPOSE:      Remap spectrum by linear interpolation
; CATEGORY:
; CALLING SEQUENCE:
;       In = shift_s ( Vector , dx )
; INPUTS:
;       Vector	= Image to be shifted
;	Dx	= Shift in X-direction
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       In
; COMMON BLOCKS:
;       None.
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; MODIFICATION HISTORY:
;     Written by Victor Afanasiev, Special Astrophisical observatory,  July, 1999.
;-

on_error,2
 ax=0
in=shift(vector,rfix(dx))
dx=dx-rfix(dx)
if dx lt 0 then ax=in -shift(in,-1)
if dx gt 0 then ax=shift(in,1)-in
in=IN+ax*DX
return,IN
END
