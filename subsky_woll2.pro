pro subsky_WOLL2,dir,FLAT=flat,NSDEG=NSdeg,NSCUT=NScut,PLOT=plot,YSHIFT=yshift
if not(keyword_set(NSdeg)) then  NSdeg=2
if not(keyword_set(NScut)) then  NScut=3
; extraction spectra WOLL2

cube=readfits(dir+'obj_lin.fts',h,/silent)
bin=sxpar(h,'BINNING') & bin=str_sep(bin,'x') & bin=FLOAT(bin(1)) & print,'bin',bin
co=1.5
a=size(cube)
Nx=a(1) & Ny=a(2) & Npol=a(3) & Nmax=a(4)  & Ncube=a(5)
Nexp=fltarr(Ncube)
for k=0,Ncube-1 do Nexp(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))

map=fltarr(Nx,Ny)
;����������� �������� ���� � ���������������� ����������� �������
if keyword_SET(flat) then begin
avg_flat=readfits(dir+'avg_flat.fts',/silent)
avg_ratio=readfits(dir+'avg_ratio.fts',/silent)
for j=0,Ncube-1 do begin
for i=0,Nexp(j)-1 do begin
for k=0,Npol-1 do cube(*,*,k,i,j)=cube(*,*,k,i,j)/avg_flat(*,*,k)
for k=0,1      do cube(*,*,2*k,i,j)=cube(*,*,2*k,i,j)*avg_ratio(*,*,k)
endfor & endfor
endif
titl=['object','unpolarized star','polarized star']
;��������� ���� ������� ����

;writefits,dir+'cube-flat.fts',cube,h

p=0
if keyword_set(plot) then p=1
 for j=0,Ncube-1 do begin
	for i=0,Nexp(j)-1 do begin
y=findgen(Ny)  & w=30/2.0   ;!!! default w=10/2.0

a=[0,500,50,4]
if p eq 1 then Window,3,xsize=100+a(1)*2 ,ysize=(a(2)+1)*a(3)*co,title=dir+',    '+titl(j)+',   exposure'+string(i+1)
robomean,cube(*,*,*,i,j),3,0.5,avg_map,rms_map

for k=0,Npol-1  do begin
map=cube(*,*,k,i,j)
;����������� ������� ���������� ������� ����
V=total(map(Nx/2-w:Nx/2+w,*),1)
robomean,V,NScut,0.5,avg_V,rms_V
R=where(V lt avg_V+3*rms_V)   ;!!!! 3
if p eq 1 then tv,255-bytscl(congrid(map,a(1),a(2)*co),avg_map-rms_map,avg_map+5*rms_map),100,(a(2)+1)*k*co
if p eq 1 then begin
;plot,V,y,position=[0,(a(2)*bin/4.0+1)*k*co,100,(a(2)*bin/4.0+1)*(k+1)*co],/device,/noerase,color=2^16-1,xcharsize=1e-5
plot,V,y,position=[0,(a(2)+1)*k*co,100,(a(2)+1)*(k+1)*co],/device,/noerase,color=2^16-1,xcharsize=1e-5,XST=1,yst=1
oplot,V(R),y(R),psym=6,symsize=0.25,color=1e5
endif
for x=0,Nx-1 do begin
f=goodpoly(y(R),map(x,R),NSdeg,1,Yfit)
	if k+i+j eq 0 then begin
		if (x gt 1500 and x lt 1550) then print, 'coef', f
	endif
sky=0 & for s=0,NSdeg do sky=sky+f(s)*y^s
map(x,*)=map(x,*)-sky
cube(x,*,k,i,j)=map(x,*)
endfor




;if p eq 1 then tv,255-bytscl(congrid(map,a(1)/4,a(2)*co),-rms_map,5*rms_map),100+a(1)/4,(a(2)+1)*bin/4.0*k*co
if p eq 1 then tv,255-bytscl(congrid(map,a(1),a(2)*co),-rms_map,5*rms_map),100+a(1),(a(2)+1)*k*co ELSE $
print,'subsky ',titl(j)+',   exposure'+string(i+1),NSdeg

endfor
;if p eq 1 then wait,0.05

;print,'subsky ',titl(j)+',   exposure'+string(i+1),NSdeg
;wait,0.1
	endfor
endfor

;��������� �������� ����� � �������� � ������
if keyword_set(yshift) then begin
		for t=0,Ncube-1 do begin
	for e=0,Nexp(t)-1 do begin
for p=0,Npol-1 do begin
;ma=obj_cube(*,*,p,e,t)
cube(*,*,p,e,t)=corr_tra(cube(*,*,p,e,t),Ndeg=1,win=200)
cube(*,*,p,e,t)=corr_tra(cube(*,*,p,e,t),/plot,Ndeg=2,win=200)
wait,0.1
endfor
	endfor
		endfor
endif

if keyword_set(flat) then sxaddhist,'flat-field corrected',h
if keyword_set(Yshift) then sxaddhist,'shift along slit corrected',h
;
writefits,dir+'obj-sky.fts',cube,h

NSfin:

END

wdir='h:\red_data.pol\Sy1\Mkn79_151106\'
wdir='h:\red_data.pol\Sy1\NGC3516_170131\'
wdir='h:\red_data.pol\Sy1\IRAS03450_141020\'
wdir='h:\red_data.pol\Sy1\NGC5548_150325\'
neon=readfits(wdir+'neon_lin.fts')

a=size(neon)
neon=total(neon(*,a(2)*0.5-a(2)/8:a(2)*0.5+a(2)/8,*),2)
M=10 & w=3 & xcross=findgen(2*M+1)-M

for k=0,a(3)-1 do begin
ycross=cross_norm(neon(*,k),neon(*,0),M)
f=goodpoly(xcross(M-w:M+w),ycross(M-w:M+w),2,3,Fit)
dx=-f(1)/f(2)/2 & print,dx
endfor
subsky_WOLL2,wdir,/flat,NSdeg=1,NScut=1,/plot;,/Yshift

;endfor
END
