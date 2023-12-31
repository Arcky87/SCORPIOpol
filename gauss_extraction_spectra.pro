;spectra extraction Woll2 v2

function gauss_extraction_spectra,ima
a=size(ima) & Ndeg=1
Nx=a(1) & Ny=a(2)

spectra=fltarr(Nx)
y=findgen(Ny) & x=findgen(Nx)
 ;& wx=100

;knots of interpolation
Nk=fix(Nx/10)-5 & knots=fltarr(Nk)
for i=0,Nk-1 do knots(i)=x(10*i)+4

;gauss
H=fltarr(Nk) & W=fltarr(Nk)
for i=0, Nk-1 do begin
	slice=ima(knots(i),*)
	res=gaussfit(y(Ny/2-20:Ny/2+20),slice(Ny/2-20:Ny/2+20),a)
	H(i)=a(1) & W(i)=a(2)
endfor

;spline
y_tra=fltarr(Nx) & w_tra=fltarr(Nx)
H=interpol(H,knots,x) & W=interpol(W,knots,x)
		f=goodpoly(x,H,Ndeg,1)
		fit=0 & for i=0,Ndeg do  fit=fit+f(i)*x^i
		y_tra=fit
		f=goodpoly(x,W,Ndeg,1)
		fit=0 & for i=0,Ndeg do  fit=fit+f(i)*x^i
		w_tra=fit

;summ in 5 guassian width
spectra=fltarr(Nx)
for i=0, Nx-1 do begin
	slice=ima(x(i),(y_tra(i)-7*w_tra(i)):(y_tra(i)+7*w_tra(i)))
	spectra(i)=total(slice)
endfor

return, spectra

end

LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path='e:\LOGS\')
wdir=DEF_WDIR(LOGFILE)

cube_ima=readfits(wdir+'obj-sky.fts',h)
a=size(cube_ima) & ima=cube_ima(50:a(1)-1,*,*,*,*)

cube_spectra=fltarr(a(1)-50,a(3),a(4),a(5))

for o=0,2 do begin
	for e=0,15 do begin
		for p=0,3 do begin
		print, o, e, p
			cube_spectra(*,p,e,o)=gauss_extraction_spectra(ima(*,*,p,e,o))
		endfor
	endfor
endfor

cgdisplay, wid=0
!p.multi=[0,1,4]
	for e=0,4 do begin
			cgplot,cube_spectra(*,0,e,0), yrange=[-500,3000], title=string(e)
			cgplot,cube_spectra(*,1,e,0), yrange=[-500,3000]
			cgplot,cube_spectra(*,2,e,0), yrange=[-500,3000]
			cgplot,cube_spectra(*,3,e,0), yrange=[-500,3000]
			wait,1
	endfor

writefits,wdir+'spectra.fit',cube_spectra,h

end