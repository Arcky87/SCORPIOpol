function  angle_calculation,Q,U
Teta=ATAN(U/Q)*180/!PI
if Q gt 0 and U ge 0 then TETA_0=0
if Q gt 0 and U lt 0 then TETA_0=360
if Q lt 0 then TETA_0=180
if Q eq 0 and U gt 0 then TETA_0=90
if Q eq 0 and U lt 0 then TETA_0=270
return, (TETA+TETA_0)/2
END

;calibration angle

path='/home/elias/SCORPIO/sppol_pipeline_v2023.8/'
dir='Sy1'
dir='AGN'
dir='tinatin'
FILTER='R'
value_F=readfits(path+filter+'.fts',H)
 wave_F=FINDGEN(sxpar(h,'NAXIS1'))*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')

log_dir=path+'LOGS/'

;выбор LOG-файлов
LOGFILE=FILE_SEARCH(log_dir+'*.txt')
N=N_elements(LOGFILE) & ANALYZER=intarr(N)
for k=0,N-1 do begin
LOG=read_table(LOGFILE(k))
ANALYZER(k)=FIX(strmid(sxpar(LOG,'ANALYZER'),10,1))
endfor
R=where(ANALYZER EQ 2, M) & LOGFILE=LOGFILE(R)
wdir=strarr(M)
for k=0,M-1 do begin
tmp=str_sep(sxpar(read_table(LOGFILE(k)),'w_dir'),'/')
wdir(k)=path+dir+'/'+tmp(N_elements(tmp)-2)+'/'

endfor
openw,1,path+dir+'_'+filter+'.txt'
FOR L=0,M-1 DO BEGIN
if FILE_TEST(wdir(L)+'flat_lin.fts') EQ 0 THEN GOTO,FIN
if slope_flat(wdir(L)) gt 0 then type=' new' ELSE type=' old'

cube=readfits(wdir(L)+'obj_lin.fts',h,/silent)
Nx=sxpar(h,'Naxis1')& Ny=sxpar(h,'Naxis2')& Npol=4 & Nc=3
Nexp=intarr(Nc) & Name=strarr(Nc) & PA=fltarr(Nc)
for k=0,Nc-1 do begin
num=string(k+1,format='(I1)')
Nexp(k)=sxpar(h,'NUMEXP'+num)
Name(k)=sxpar(h,'NAME'+num)
PA(k)=sxpar(h,'PA'+num)
endfor
;формирование кривой пропускания фильтра
wave=findgen(Nx)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
V= INTERPOL(value_F,wave_F,wave,/SPLINE)
V=smooth(V,10)
R=where(V lt 0,ind) & if ind gt 0 then V(R)=0
;*********
plt=0
;*********
if plt eq 1 then window,0,xsize=400,ysize=250
!P.multi=[0,1,1]
if plt eq 1 then plot,wave,V,xst=1
cube_V=fltarr(Ny,Npol,Nexp(0),Nc)
			for j=0,Nexp(0)-1 do begin
		for ky=0,Ny-1 do begin
	for kp=0,Npol-1 do begin
for kc=0,Nc-1 do begin
cube_V(ky,kp,j,kc)=TOTAL(cube(*,ky,kp,j,kc)*V(*))
endfor
	endfor
		endfor
			endfor
cube_pol=fltarr(Npol,Nexp(0),Nc)
y=findgen(Ny) & V=y

T=0
!P.multi=[0,1,4]
		for kc=0,Nc-1 do begin
	for j=0,Nexp(kc)-1 do begin
if plt eq 1 then window,1,xsize=400,ysize=800,title='cube'+string(kc+1,format='(I3)')+'  exp'+string(j+1,format='(I3)')
for k=0,Npol-1 do begin
;вычитание фона
V(*)=cube_V(*,k,j,kc)
f=goodpoly([y(0:20),y(Ny-20:Ny-1)],[V(0:20),V(Ny-20:Ny-1)],2,3)
fon=f(0)+f(1)*y+f(2)*y^2
V=V-fon
yfit=mpfitpeak(y,V,A,/moffat)
cube_pol(k,j,kc)=total(yfit)
if plt eq 1 then plot,y,V,xst=1,title=string(cube_pol(k,j,kc))
if plt eq 1 then oplot,y,yfit,color=3e5
endfor
wait,T
	endfor
		wait,T
			endfor
; определение чувствительности поляризационных каналов
kQ=fltarr(Nexp(1)) & kU=kQ

for k=0,Nexp(1)-1 do begin
kQ(k)=cube_pol(0,k,1)/cube_pol(1,k,1)
kU(k)=cube_pol(2,k,1)/cube_pol(3,k,1)
;print,kQ(k),kU(k)
endfor
robomean,kQ,3,0.5,avg_kQ,rms_kQ
robomean,kU,3,0.5,avg_kU,rms_kU

printf,1,wdir(L),' ',type,'  ',Name(2),PA(2)
;if  avg_kQ  eq 0 OR avg_kU then goto,fin

;avg_kQ=1 & avg_kU=1
;нуль-пункт поляризации
kc=1
Q=fltarr(Nexp(kc)) & U=Q & P=Q & TETA=P
for k=0,Nexp(kc)-1 do begin
Q(k)=(cube_pol(0,k,kc)-cube_pol(1,k,kc)*avg_kQ)/(cube_pol(0,k,kc)+cube_pol(1,k,kc)*avg_kQ)
U(k)=(cube_pol(2,k,kc)-cube_pol(3,k,kc)*avg_kU)/(cube_pol(2,k,kc)+cube_pol(3,k,kc)*avg_kU)
endfor

;print,avg_U_0*100,rms_U_0*100
;вычисление параметров поляризации  стандарта поляризации
kc=2
Q=fltarr(Nexp(kc)) & U=Q & P=Q & TETA=P
for k=0,Nexp(kc)-1 do begin
Q(k)=(cube_pol(0,k,kc)-cube_pol(1,k,kc)*avg_kQ)/(cube_pol(0,k,kc)+cube_pol(1,k,kc)*avg_kQ)
U(k)=(cube_pol(2,k,kc)-cube_pol(3,k,kc)*avg_kU)/(cube_pol(2,k,kc)+cube_pol(3,k,kc)*avg_kU)
Q(k)=Q(k)-avg_Q_0 & U(k)=U(k)-avg_U_0

P(k)=SQRT(Q(k)^2+U(k)^2)
TETA(k)= angle_calculation(Q(k),U(k))
;print,Q(k)*100,U(k)*100,P(k)*100,TETA(k)
endfor
robomean,Q,3,0.5,avg_Q,rms_Q
robomean,U,3,0.5,avg_U,rms_U
printf,1,'kQ=',avg_kQ,rms_kQ,'kU=',avg_kU,rms_kU,format='(A5,2F7.4,A6,2F7.4)'
printf,1,' Q=',avg_Q*100,rms_Q*100,'U=',avg_U*100,rms_U*100,format='(A5,2F7.2,A6,2F7.2)'

robomean,P*100,3,0.5,avg_P,rms_P
robomean,TETA,3,0.5,avg_TETA,rms_TETA
;printf,1,'P=',avg_P,rms_P,'TETA=',avg_TETA,rms_TETA,format='(A5,2F7.2,A6,2F8.1)'
fin:
print,L
wait,2
ENDFOR
close,1
end

