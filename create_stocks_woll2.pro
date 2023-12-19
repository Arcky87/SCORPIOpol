;test WOLLASTON2

pro create_stocks_WOLL2,wdir,CUT=cut,PLOT=plot
if not(keyword_set(cut)) then cut=100
;проверка шкалы длин волн
neon=readfits(wdir+'neon_lin.fts',h)
Nx=sxpar(h,'NAXIS1')
Ny=sxpar(h,'NAXIS2')
yc=Ny/2  & wy=20
vector=total(neon(*,yc-wy:yc+wy,*),2)
M=10 & xcross=findgen(2*M+1)-M
for k=0,3 do begin
ycross=cross_norm(vector(*,k),vector(*,0),M)
gau=gaussfit(xcross,ycross,G)
print,G(1),g(2)*2.345/sqrt(2)
endfor
;калибровка нулевого стандакта
s=readfits(wdir+'avg_spectra.fit',h)
spectra=readfits(wdir+'spectra.fit',hs)

Nx=sxpar(h,'NAXIS1') & Npol=sxpar(h,'NAXIS2') & Ncube=sxpar(h,'NAXIS3')
name=strarr(Ncube) & PA=fltarr(Ncube) & Nexp=intarr(Ncube)
for k=0,2 do begin
Nexp(k)=sxpar(hs,'NUMEXP'+string(k+1,format='(I1)'))
Name(k)=sxpar(hs,'NAME'+string(k+1,format='(I1)'))
PA(k)=sxpar(hs,'PA'+string(k+1,format='(I1)'))
endfor

wx=20
x=findgen(Nx) & wave=x*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')

beg=cut
cut=Nx-cut
;& deg=5
R_01=s(*,0,1)/s(*,1,1)
R_23=s(*,2,1)/s(*,3,1)
;сглаживание
fit=LOWESS(x(beg:cut),R_01(beg:cut),Nx/2 ,3,2)
R01=INTERPOL( fit,x(beg:cut),x)
fit=LOWESS(x(beg:cut),R_23(beg:cut),Nx/2 ,3,2)
R23=INTERPOL( fit,x(beg:cut),x)
if keyword_set(plot) then begin
window,0,xsize=600,ysize=427
!P.multi=[0,1,2]
plot,wave,R_01,xst=1,charsize=1,yrange=[0.5,1.5],$
	ytitle='R01',xtitle='Wavelength, A',yst=1,title='unpolarized star '+name(1)
xyouts,wave(Nx/2),1.3,'R01=flux(0)/flux(90)',align=0.5
oplot,wave,R01,color=3e5,thick=2
plot,wave,R_23,xst=1,charsize=1,yrange=[0.5,1.5],$
	ytitle='R23',xtitle='Wavelength, A',yst=1
xyouts,wave(Nx/2),1.3,'R23=flux(45)/flux(135)',align=0.5
oplot,wave,R23,color=3e5,thick=2
endif
stoks=fltarr(Nx,4,3) ;!!!
cube_stoks=fltarr(Nx,4,Nexp(0),3) ;!!!


Q=fltarr(Nx,3) & U=Q

for k=0,Ncube-1 do begin
stoks(*,0,k)=s(*,0,k)+R01*s(*,1,k)+s(*,2,k)+R23*s(*,3,k)
stoks(*,1,k)=(s(*,0,k)-R01*s(*,1,k))/(s(*,0,k)+R01*s(*,1,k))
stoks(*,2,k)=(s(*,2,k)-R23*s(*,3,k))/(s(*,2,k)+R23*s(*,3,k))
for j=0,Nexp(k)-1 do begin
cube_stoks(*,0,j,k)= spectra(*,0,j,k)+R01*spectra(*,1,j,k)+spectra(*,2,j,k)+R23*spectra(*,3,j,k)
cube_stoks(*,1,j,k)=(spectra(*,0,j,k)-R01*spectra(*,1,j,k))/(spectra(*,0,j,k)+R01*spectra(*,1,j,k))
cube_stoks(*,2,j,k)=(spectra(*,2,j,k)-R23*spectra(*,3,j,k))/(spectra(*,2,j,k)+R23*spectra(*,3,j,k))
endfor
endfor

Q_null=Smooth(stoks(*,1,1),wx,/edge_truncate)
U_null=SMooth(stoks(*,2,1),wx,/edge_truncate)
;***************************************
;запись результата
writefits,wdir+'avg_stoks.fit',stoks,h
writefits,wdir+'stoks.fit',cube_stoks,hs
;***************************************

Q_null=0
U_null=0
for k=0,Ncube-1 do begin
stoks(*,1,k)=stoks(*,1,k)-Q_null
stoks(*,2,k)=stoks(*,2,k)-U_null
for j=0,Nexp(k)-1 do begin
cube_stoks(*,1,j,k)=cube_stoks(*,1,j,k)-Q_null
cube_stoks(*,2,j,k)=cube_stoks(*,2,j,k)-U_null
endfor
endfor
;
window,1,xsize=600,ysize=427, title='unpolarized star after depol cor'
T=1
robomean,stoks(Nx/2-Nx/4:Nx/2+Nx/4,1,T)*100,2,0.5,avg_val,rms_val
plot,stoks(*,1,T)*100,yrange=[-5,5],xst=1,charsize=1,$
	title='Q='+string(avg_val,format='(F7.2)')+$
	'!9 +!3'+string(rms_val,format='(F5.2)')

oplot,[0,Nx],[0,0],color=3e5,thick=2
robomean,stoks(Nx/2-Nx/4:Nx/2+Nx/4,2,T)*100,2,0.5,avg_val,rms_val
plot,stoks(*,2,T)*100,yrange=[-5,5],xst=1,charsize=1,$
	title='U='+string(avg_val,format='(F7.2)')+$
	'!9 +!3'+string(rms_val,format='(F5.2)')

oplot,[0,Nx],[0,0],color=3e5,thick=2

;;***********
FLUX=fltarr(Nx,3) & avg_Q=fltarr(NX)  & avg_U=avg_Q
FLUX(*,*)=stoks(*,0,*)
T=0

;box=5
;avg_Q=LOWESS(findgen(Nx),Q(*,T),box,2,2)

avg_Q(*)=stoks(*,1,T)-Q_null
;avg_U=LOWESS(findgen(Nx),U(*,T),box,2,2)
avg_U(*)=stoks(*,2,T)-U_null
avg_P=sqrt(avg_Q^2+avg_U^2)
avg_FI=fltarr(Nx)
for x=0,Nx-1 do  avg_FI(x)=calc_atan(avg_Q(x),avg_U(x))/2; & R=where(avg_FI gt 178) & avg_FI(R)=avg_FI(R)-180

Window,2,xsize=600,ysize=900
!P.background=2^16-1 & !P.color=0
!P.multi=[0,1,5]
!P.charsize=2
R=where(wave gt 5000 and wave lt 6000)
RR=where(wave gt 4500 and wave lt 7500,IND) & PRINT,IND
RR=findgen(Nx)
PRINT,AVG_p
robomean,avg_P(R)*100,3,0.5,P,err_P
robomean,avg_FI(R) ,3,0.5,F,err_F
plot,wave(RR),FLUX(RR,T),xst=1,$
title=name(T)+'   PAslit='+string(PA(T),format='(F7.1)')
plot,wave(RR),avg_Q(RR)*100,xst=1,yrange=[-5,5],yst=1
plot,wave(RR),avg_U(RR)*100,xst=1,yrange=[-5,5],yst=1
plot,wave(RR),avg_P(RR)*100,xst=1,yrange=[0,10],yst=1,$
	title='P(V)='+string(P,format='(F5.2)')+'!9 +!3'+string(err_P,format='(F5.2)')
plot,wave(RR),avg_FI(RR),xst=1,title='!9P!3(V)='+string(F,format='(F7.1)')+'!9 +!3'+string(err_F,format='(F4.1)'),$
subtitle=wdir
;
;

FLUX=fltarr(Nx,3) & avg_Q=fltarr(NX)  & avg_U=avg_Q
FLUX(*,*)=stoks(*,0,*)
T=2

;box=5
;avg_Q=LOWESS(findgen(Nx),Q(*,T),box,2,2)

avg_Q(*)=stoks(*,1,T)-Q_null
;avg_U=LOWESS(findgen(Nx),U(*,T),box,2,2)
avg_U(*)=stoks(*,2,T)-U_null
avg_P=sqrt(avg_Q^2+avg_U^2)
avg_FI=fltarr(Nx)
;for x=0,Nx-1 do  avg_FI(x)=calc_atan(avg_Q(x),avg_U(x))/2; & R=where(avg_FI gt 178) & avg_FI(R)=avg_FI(R)-180
for x=0,Nx-1 do  avg_FI(x)=arctan(avg_Q(x),avg_U(x))

Window,3,xsize=600,ysize=900
!P.background=2^16-1 & !P.color=0
!P.multi=[0,1,5]
!P.charsize=2
R=where(wave gt 5000 and wave lt 6000)
RR=where(wave gt 4500 and wave lt 7500,IND) & PRINT,IND
RR=findgen(Nx)
PRINT,AVG_p
robomean,avg_P(R)*100,3,0.5,P,err_P
robomean,avg_FI(R) ,3,0.5,F,err_F
plot,wave(RR),FLUX(RR,T),xst=1,$
title=name(T)+'   PAslit='+string(PA(T),format='(F7.1)')
plot,wave(RR),avg_Q(RR)*100,xst=1,yrange=[-5,5],yst=1
plot,wave(RR),avg_U(RR)*100,xst=1,yrange=[-5,5],yst=1
plot,wave(RR),avg_P(RR)*100,xst=1,yrange=[0,10],yst=1,$
	title='P(V)='+string(P,format='(F5.2)')+'!9 +!3'+string(err_P,format='(F5.2)')
plot,wave(RR),avg_FI(RR),xst=1,title='!9P!3(V)='+string(F,format='(F7.1)')+'!9 +!3'+string(err_F,format='(F4.1)'),$
subtitle=wdir

END


path='d:\'
LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path=path+'LOGS')
;LOGFILE=path+'LOGS\'+'E1821+643_140324.txt'
tmp=str_sep(LOGFILE,'\')
path=tmp(0)
for j=1,N_elements(tmp)-3 do begin
path=path+'\'+tmp(j)
endfor
path=path+'\'
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
wdir=path+wdir(N_elements(wdir)-2)+'\'
print,wdir
cut=450
create_stocks_WOLL2,wdir,/plot, CUT=cut
end