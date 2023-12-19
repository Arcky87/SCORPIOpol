;correct calculation arctangens
function  calc_atan,x,y
N=N_elements(x)
fi=fltarr(N)
for k=0,N-1 do begin
if x(k) eq 0 and y(k) gt 0 then fi(k)=90
if x(k) eq 0 and y(k) lt 0 then fi(k)=270
if x(k) ne 0 then begin
c=abs(y(k)/x(k))
fi(k)=atan(c)*180/!PI
if y(k) ge 0 and x(k) lt 0 then fi(k)=180-fi(k)
if y(k) lt 0 and x(k) lt 0 then fi(k)=fi(k)+180
if y(k) lt 0 and x(k) gt 0 then fi(k)=360-fi(k)
endif
endfor
return,fi
end
