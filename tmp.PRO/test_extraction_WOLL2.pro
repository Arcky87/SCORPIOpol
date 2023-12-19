;test_extraction WOLLASTON-2
wdir=read_table('h:\red_data.pol\wdir_old.txt')
Ndir=N_elements(wdir)

FOR D=18,Ndir-1 DO BEGIN
print,systime()
print,wdir(D)
wait,1
;goto,cont
obj=readfits(wdir(D)+'obj-sky.fts',h,/silent)
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2') & Npol=sxpar(h,'NAXIS3') & Nex=sxpar(h,'NAXIS4') & Ntarget=sxpar(h,'NAXIS5')
Nexp=intarr(Ntarget)
for k=0,2 do Nexp(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))

;print,Nexp
angle=['0','45','90','135']
spectra=fltarr(Nx,Npol,Nex,Ntarget)
avg_spectra=fltarr(Nx,Npol,Ntarget)
ys=Ny/2 & wy=10
y=findgen(Ny)
for T=0,Ntarget-1 do begin
 	for P=0,Npol-1 do begin
		for j=0,Nexp(T)-1 do begin
			for x=0,Nx-1 do begin
			yfit= MPFITPEAK(y(ys-wy:ys+wy), obj(x,ys-wy:ys+wy,P,j,T), A)
			spectra(x,P,j,T)=total(Yfit)

		endfor
;print,'target='+string(T+1,format='(I1)')+' angle='+angle(p)+' exp='+string(J+1,format='(I2)')
;wait,0.02
		;window,2
		;!P.multi=[0,1,1]
	;	plot,spectra(*,P,j,T),xst=1

	endfor
for x=0,Nx-1 do  avg_spectra(x,P,T)=median(spectra(x,P,0:Nexp(T)-1,T))*Nexp(T)
endfor

endfor
writefits,wdir(D)+'avg_spectra.fit',avg_spectra,h
print,systime()
ENDFOR
end
cont:
s=readfits(wdir+'avg_spectra.fit',h)
Nx=sxpar(h,'NAXIS1') & Npol=sxpar(h,'NAXIS2') & Ntarget=sxpar(h,'NAXIS3')

Nx=1800
s=s(0:1799,*,*)
wave=findgen(Nx)*sxpar(h,'CDELT1') +sxpar(h,'CRVAL1')
window,2,xsize=1000,ysize=800
!P.multi=[0,1,4]
for k=0,Npol-1 do plot,s(*,k,2),xst=1
Q=fltarr(Nx,Ntarget)   & U=Q
for k=0,Ntarget-1 do begin

Q(*,k)=(s(*,0,k)-s(*,1,k))/(s(*,0,k)+s(*,1,k))
U(*,k)=(s(*,2,k)-s(*,3,k))/(s(*,2,k)+s(*,3,k))
endfor
;нуль
Q_null=LOWESS(findgen(Nx),Q(*,1),Nx/2 ,2,2)
U_null=LOWESS(findgen(Nx),U(*,1),Nx/2 ,2,2)
window,0
!P.multi=[0,1,2]
plot,Q_null,xst=1
plot,U_null,xst=1
for k=0,Ntarget-1 do begin
Q(*,k)=Q(*,k)-Q_null & U(*,k)=U(*,k)-U_null
P=SQRT(Q^2+U^2)   &
FI=fltarr(Nx,3)

FLUX=total(s,2)
ENDFOR
T=0

box=20
avg_Q=LOWESS(findgen(Nx),Q(*,T),box,2,2)

avg_U=LOWESS(findgen(Nx),U(*,T),box,2,2)
avg_P=sqrt(avg_Q^2+avg_U^2)
avg_FI=calc_atan(avg_Q,avg_U)/2; & R=where(avg_FI gt 178) & avg_FI(R)=avg_FI(R)-180

Window,3,xsize=600,ysize=900
!P.multi=[0,1,5]
!P.charsize=2

plot,wave,FLUX(*,T),xst=1
plot,wave,avg_Q,xst=1
plot,wave,avg_U,xst=1
plot,wave,avg_P,xst=1
plot,wave,avg_FI,xst=1
R=where(wave gt 5000 and wave lt 6000)
robomean,avg_P(R),3,0.5,P,err_P
print,P,err_P

end


