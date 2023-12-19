;test ext
function  extraction,ima,PLOT=plot,WTITLE=wtitle,YPOS=ypos,PSF=psf
if not(keyword_set(Wtitle)) then wtitle=''
;экстракция спектрoв из изображения
a=size(ima)
y=findgen(a(2))
spectra=fltarr(a(1),2)
wy=20
if total(ima) eq 0 then goto,fin
for i=0,1 do begin
for x=0,a(1)-1 do begin
if keyword_set(ypos) then ylim=[ypos(i)-wy,ypos(i)+wy] $
else ylim=[a(2)/2*i,a(2)/2*(i+1)-1]
if not(keyword_set(psf)) then begin
gau=gaussfit(y(ylim(0):ylim(1)),ima(x,ylim(0):ylim(1)),G,Nterms=4)
spectra(x,i)=total(gau-g(3))
endif
if keyword_set(psf) then begin
vector=fltarr(ylim(1)-ylim(0)+1)
kernel=fltarr(ylim(1)-ylim(0)+1)
vector(*)=ima(x,ylim(0):ylim(1))
kernel(*)=psf(x,ylim(0):ylim(1))
spectra(x,i)=max(convol(vector,kernel))/sqrt(total(kernel^2))*total(kernel)
endif
endfor
endfor
if keyword_set(plot) then begin
window,plot,title=wtitle
!P.multi=[0,1,2]
plot,spectra(*,0),xst=1
plot,spectra(*,1),xst=1
!P.multi=[0,1,2]
endif
fin:
return,spectra
end

LOGFILE='g:\spectraPOL.log\GRW+70_100813.txt'
LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path='h:\spectraPOL.log')
wdir=sxpar(read_table(LOGFILE),'w_dir')
cube=readfits(wdir+'obj-sky.fts',h)
a=size(cube)
spectra=fltarr(a(1),2,a(3),a(4))
ji=1 & j=3

spectra_gau=fltarr(a(1))  & spectra_opt=spectra_gau

y=[100,170]
;psf=total(total(cube(a(1)/2-a(1)/8:a(1)/2+a(1)/8,1)-1,*,0),1),2)
;psf=fltarr(a(1),y(1)-y(0))
;for x=0,a(1)-1 do  psf(x,*)=total(cube(x,y(0):y(1)-1,i,3),1)
;psf=total(cube(*,*,*,j),3)
;psf=psf(y(0):y(1)-1)
;for x=0,a(1)-1 do begin
;psf=cube(x,y(0):y(1)-1,i,0)
;vector=fltarr(y(1)-y(0)); & psf=fltarr(y(1)-y(0))
;psf=cube(x,y(0):y(1)-1,i,0)
;vector(*)=cube( x,y(0):y(1)-1,i,j)
;kernel=fltarr(y(1)-y(0))
;kernel(*)=total(psf,1)
;kernel=kernel/max(kernel)
 ;plot,vector,xst=1,psym=10
;gau=gaussfit(findgen(y(1)-y(0)),vector,G)
;oplot,gau
;spectra_gau(x)=total(gau-G(3))
;spectra_opt(x)=max(convol(vector,kernel))/sqrt(total(kernel^2))*total(kernel)
;print,total(gau-G(3)),max(convol(vector,vector))/sqrt(total(vector^2)),max(convol(vector,vector))/sqrt(total(vector^2))/total(gau)
;plot,convol(vector,vector),xst=1
;endfor
;Window,2
;!P.multi=[0,1,2]
;plot,spectra_gau/LOWESS(findgen(a(1)),spectra_gau,a(1)/8,3),xst=1,yrange=[0.6,1.1]
;plot,spectra_opt/LOWESS(findgen(a(1)),spectra_opt,a(1)/8,3),xst=1,yrange=[0.6,1.1]
;END

for j=0,a(4)-1 do begin
for i=0,a(3)-1 do begin
spectra(*,*,i,j)=extraction(cube(*,*,i,j),plot=2,wtitle='i='+string(i,format='(I2)')+'  j='+string(j,format='(I2)'),ypos=[80,240]) ;ELSE $
;endif else $
;spectra(*,*,i,j)=extraction(cube(*,*,i,j),plot=2,wtitle='i='+string(i,format='(I2)')+'  j='+string(j,format='(I2)'),ypos=[135,300],psf=total(cube(*,*,*,0),3))
;,ypos=[137,301]),ypos=[30,193]

endfor
endfor

;sxaddpar,h,'NAME1','star near 0716+71'
writefits,wdir+'obj_spectra.fit',spectra,h
end
