;test identefication
wdir='h:\red_data.pol\mkn509_141021\'
neon=readfits(wdir+'neon.fts',h)
N=numlines(wdir+'etalon.old')
etalon=fltarr(5,N)
openr,1,wdir+'etalon.txt'
readf,1,etalon
close,1

Ns=4  & G=0.25
Ny=sxpar(h,'NAXIS2')  & wy=20 & Nx=sxpar(h,'NAXIS1') & x=findgen(Nx)
print,Nx
spectra=total(neon,2)
for j=0,Ns-1 do begin
R=where(spectra(*,j) lt 1) & spectra(R,j)=1 & spectra(*,j)=spectra(*,j)^G
fi_peak,x ,spectra(*,j),0,ipix,xpk,ypk,bkpk,ipk
fon=INTERPOL(bkpk,xpk,x)
spectra(*,j)=spectra(*,j)-fon & R=where(spectra(*,j) lt 0) & spectra(R,j)=0
endfor

window,2,xsize=1600,ysize=1000
!P.multi=[0,1,4]
for k=0,3 do begin
plot,spectra(*,k),xst=1
print,etalon(k,*)
for j=0,N-1 do oplot,[1,1]*etalon(k+1,j),[0,1]*1e3,linestyle=2
endfor
end
