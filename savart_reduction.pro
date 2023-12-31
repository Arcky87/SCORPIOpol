;Savart_reduction
dir='e:\Phaethon_171218\' & cube='S15990'
;dir='e:\Phaethon_171219\' & cube='S16000'
dir='e:\Phaethon_171226\' & cube='S16020'
file='obj' & scale=1e6

;goto,cont
N=numlines(dir+file+'.txt')
print,N
tab=intarr(10,N)
openr,1,dir+file+'.txt'
readf,1,tab
close,1

ray=['o','e']
ang=['0 deg','45 deg']
flux=fltarr(2,2,N)  & err_flux=flux
w=16 ; ���������� ����
map=fltarr(2*w,2*w,2)
S=7

dN=5
;
openw,1,dir+file+'_flux.txt'
		for j=0,N-1 do begin
	for k=0,1 do begin
file_in=dir+cube+string(tab(dN*k,j),format='(I3)')+'.fts'
ima=readfits(file_in,/silent)/scale
	x_o=tab(1+dN*k,j) & y_o=tab(2+dN*k,j)
	x_e=tab(3+dN*k,j) & y_e=tab(4+dN*k,j)
	map(*,*,0)=ima(x_o-w:x_o+w-1,y_o-w:y_o+w-1) ; ������������ ���
	yfit=MPFIT2DPEAK(map(*,*,0), A, /TILT,/moffat )
	x_o=x_o+dx & y_o=y_o+dy
	map(*,*,0)=ima(x_o-w:x_o+w-1,y_o-w:y_o+w-1) ; ������������ ��� �������� ������
	map(*,*,1)=ima(x_e-w:x_e+w-1,y_e-w:y_e+w-1) ; �������������� ���
	yfit=MPFIT2DPEAK(map(*,*,1), A, /TILT,/moffat )
	x_e=x_e+dx & y_e=y_e+dy
	map(*,*,1)=ima(x_e-w:x_e+w-1,y_e-w:y_e+w-1) ; �������������� ��� �������� ������
for i=0,1 do begin
window,0,xsize=2*w*S,ysize=2*w*S*3,title='exp='+string(j,format='(I1)')+' angle='+ang(k)+' ray='+ray(i)
robomean,map(*,*,i),3,0.5,avg_map,rms_map
tv,255-bytscl(congrid(map(*,*,i),2*w*S,2*w*S),avg_map-3*rms_map,avg_map+5*rms_map)
fit=MPFIT2DPEAK(map(*,*,i), A, /TILT,/moffat )
tv,255-bytscl(congrid(fit,2*w*S,2*w*S),avg_map-3*rms_map,avg_map+5*rms_map),0,2*w*S
tv,255-bytscl(congrid(map(*,*,i)-fit+A(0),2*w*S,2*w*S),avg_map-3*rms_map,avg_map+5*rms_map),0,4*w*S
flux(i,k,j)=total(fit)-A(0) & err_flux(i,k,j)=SQRT(total(((map(*,*,i)-fit)^2)))
;print,flux(i,k,j),err_flux(i,k,j)
endfor
	endfor
printf,1,flux(0,0,j),err_flux(0,0,j),flux(1,0,j),err_flux(1,0,j),flux(0,1,j),err_flux(0,1,j),flux(1,1,j),err_flux(1,1,j),format='(8F10.5)'
		endfor
close,1
cont:
N=numlines(dir+file+'_flux.txt')
tab=fltarr(8,N)
openr,1,dir+file+'_flux.txt'
readf,1,tab
close,1
flux=reform(tab,2,4,N )
print,flux(0,3,*)
window,0,xsize=600,ysize=900
!P.multi=[0,1,3]
!P.charsize=2
Q=(flux(0,0,*)-flux(0,1,*))/(flux(0,0,*)+flux(0,1,*))
plot,Q,xst=1,psym=6,xrange=[-1,N]
U=(flux(0,2,*)-flux(0,3,*))/(flux(0,2,*)+flux(0,3,*))
plot,U,xst=1,psym=6,xrange=[-1,N]
P=SQRT(U^2+Q^2)*SQRT(2)
plot,P,xst=1,psym=6,xrange=[-1,N]
robomean,P(2:N-1),3,0.5,avg_P,rms_P
;print,avg_P,rms_P
end
