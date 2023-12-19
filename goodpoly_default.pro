;+
; NAME:
;	goodpoly
; PURPOSE: (one line)
;	Robust fitting of a polynomial to data.
; DESCRIPTION:
;	This is a multi-pass fitting routine that fits a fixed order polynomial
;	to the input data.  After each pass, the scatter of the fit relative
;	to the fitted line is computed.  Each point is examined to see if it
;	falls beyond THRESH sigma from the line.  If is does, it is removed
;	from the data and the fit is tried again.  This will make two attempts
;	to remove bad data.
; CATEGORY:
; CALLING SEQUENCE:
;	coeff = goodpoly(x,y,order,thresh,yfit,newx,newy)
; INPUTS:
;	x      - Input dataset, independant values.
;	y      - Input dataset, dependant values.
;	order  - Order of the polynomial fit (linear = 1).
;	thresh - Sigma threshold for removing outliers.
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;	yfit   - Fitted values for y that match the input vector.
;	newx   - X values from input that were considered good.
;	newy   - Y values from input that were considered good.
;	Return value is the set of polynomial coefficients.
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;	Written 1991 Feb., Marc W. Buie, Lowell Observatory
;  93/11/12, MWB, Program fixed to return a computed y for all input x.
;-
function goodpoly,x,y,order,thresh,yfit,newx,newy

arlen=n_elements(x)
xx = x
yy = y

; Initial fit with all the data.
print,'Poly_fit with order ',order,format='(a,i2)'
;print, 'xx:   ' , xx, '   yy:   ', yy ; IY
coeff=poly_fit(xx,yy,order,yfit=yfit)
print, 'yfit:   ', yfit
chisq0=total((yy-yfit)^2)/(arlen-order)
flat = (yy-yfit)+(total(yfit)/arlen)
sigma = stdev(flat,mean)
print, 'mean:   ', mean  ; IY
print, 'sigma:   ', sigma ; IY
snr   = mean/sigma
print,'   Initial chisq=',chisq0,'   S/N=',snr,$
   format='(a,f6.1,a,f6.1)'

;Remove all points beyond threshold sigma
print, 'thresh =   ', thresh, '  sigma =  ', sigma ; IY
print, 'Flat:  ', flat ;IY
good=where( abs(flat-mean) lt thresh*sigma,goodnum)
nbad = arlen-goodnum
xx=xx(good)
yy=yy(good)
arlen=n_elements(xx)

; Second pass fit with bad points removed.
coeff=poly_fit(xx,yy,order,yfit=yfit)
chisq1=total((yy-yfit)^2)/(arlen-order)
;print, 'yy2:  ', yy, 'yfit2:  ', yfit
;print, 'total: ', total((yy-yfit)^2)
flat = (yy-yfit)+(total(yfit)/arlen)
sigma = stdev(flat,mean)
snr   = mean/sigma
;print,'   Second pass  =',chisq1,'   S/N=',snr,' :  ',nbad,' points removed',$
   format='(a,f6.1,a,f6.1,a,i3,a)'

;Remove all points beyond threshold sigma
good=where( abs(flat-mean) lt thresh*sigma,goodnum)
nbad = arlen-goodnum
xx=xx(good)
yy=yy(good)
arlen=n_elements(xx)

; Third pass fit with bad points removed.
coeff=poly_fit(xx,yy,order,yfit=yfit)
chisq2=total((yy-yfit)^2)/(arlen-order)
flat = (yy-yfit)+(total(yfit)/arlen)
sigma = stdev(flat,mean)
snr   = mean/sigma
;print,'   Third pass   =',chisq2,'   S/N=',snr,' :  ',nbad,' points removed',$
   format='(a,f6.1,a,f6.1,a,i3,a)'

newx=xx
newy=yy

yfit = poly(x,coeff)

return,coeff
end
