function  arctan,Q,U
;Teta=0
;if Q eq 0 then goto, put
Teta=(ATAN(U/Q)*180/!PI)
if Q gt 0 and U ge 0 then TETA_0=0
if Q gt 0 and U lt 0 then TETA_0=360
if Q lt 0 then TETA_0=180
;put:
if Q eq 0 and U gt 0 then TETA_0=90
if Q eq 0 and U lt 0 then TETA_0=270
return, double(TETA+TETA_0)/2.0
END