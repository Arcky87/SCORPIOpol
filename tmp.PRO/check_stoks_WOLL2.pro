;check_stoks_WOLL2


log_dir='h:\red_data.pol\'
LOGFILE=DIALOG_PICKFILE(/read,path=log_dir+'LOGS\',FILTER='*.txt')

;LOGFILE=DIALOG_PICKFILE(/read,path='h:\red_data.pol\LOGS\',FILTER='*.txt')
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
	wdir=log_dir+wdir(N_elements(wdir)-2)+'\'
;print,wdir
;wdir=log_dir+'Mkn1501_141120\'
;�������� ����� ���� ����
neon=readfits(wdir+'neon_lin.fts',h)
Nx=sxpar(h,'NAXIS1')
Ny=sxpar(h,'NAXIS2')
beg=200
cut=1800
yc=Ny/2  & wy=20
vector=total(neon(*,yc-wy:yc+wy,*),2)
M=10 & xcross=findgen(2*M+1)-M
for k=0,3 do begin
ycross=cross_norm(vector(*,k),vector(*,0),M)
gau=gaussfit(xcross,ycross,G)
print,G(1),g(2)*2.345/sqrt(2)
endfor

;s=readfits(wdir+'avg_spectra.fit',h)
spectra=readfits(wdir+'spectra.fit',hs)
Nexp=intarr(Ntarget)
for k=0,2 do Nexp(k)=sxpar(hs,'NUMEXP'+string(k+1,format='(I1)'))

Nx=sxpar(h,'NAXIS1') & Npol=sxpar(h,'NAXIS2') & Ntarget=sxpar(h,'NAXIS3')
name=strarr(3) & PA=fltarr(3)
for j=0,2 do begin
name(j)=sxpar(h,'NAME'+string(j+1,format='(I1)'))
PA(j)=sxpar(h,'PA'+string(j+1,format='(I1)'))

endfor

window,1,xsize=600,ysize=427
!P.multi=[0,1,2]
x=findgen(Nx)

;cut=4200
;& deg=5
R01=s(*,0,1)/s(*,1,1)
R23=s(*,2,1)/s(*,3,1)
plot,R01,xst=1,charsize=1,color=1e5,yrange=[0,3]
fit=LOWESS(x(beg:cut),R01(beg:cut),Nx/2 ,3,1)
R01=INTERPOL( fit,x(beg:cut),x)
oplot,x,R01
plot,R23,xst=1,charsize=1,color=1e5,yrange=[0,3]
fit=LOWESS(x(beg:cut),R23(beg:cut),Nx/2 ,3,1)
R23=INTERPOL( fit,x(beg:cut),x)
oplot,x,R23

;
stoks=fltarr(Nx,4,3)
cube_stoks=fltarr(Nx,4,Nexp(0),3)
;Nx=2291
;s=s(200:Nx-200,*,*)
a=size(s) & Nx=a(1) & x=findgen(Nx)
wave=findgen(Nx)*sxpar(h,'CDELT1') +sxpar(h,'CRVAL1')
Q=fltarr(Nx,3) & U=Q
;R01=1 & R23=1
for k=0,Ntarget-1 do begin
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
Q_null=0 & U_null=0

for k=0,Ntarget-1 do begin
stoks(*,1,k)=stoks(*,1,k)-Q_null
stoks(*,2,k)=stoks(*,2,k)-U_null
for j=0,Nexp(k)-1 do begin
cube_stoks(*,1,j,k)=cube_stoks(*,1,j,k)-Q_null
cube_stoks(*,2,j,k)=cube_stoks(*,2,j,k)-U_null
endfor
endfor



writefits,wdir+'avg_stoks.fit',stoks,h
writefits,wdir+'stoks.fit',cube_stoks,hs
wx=50


window,0,xsize=600,ysize=427
T=1
robomean,stoks(Nx/2-Nx/4:Nx/2+Nx/4,1,T)*100,2,0.5,avg_val,rms_val
plot,stoks(*,1,T)*100,yrange=[-5,5],xst=1,charsize=1,$
	title='Q='+string(avg_val,format='(F7.2)')+$
	'!9 +!3'+string(rms_val,format='(F5.2)')

oplot,[0,Nx],[0,0],linestyle=2
robomean,stoks(Nx/2-Nx/4:Nx/2+Nx/4,2,T)*100,2,0.5,avg_val,rms_val
plot,stoks(*,2,T)*100,yrange=[-5,5],xst=1,charsize=1,$
	title='U='+string(avg_val,format='(F7.2)')+$
	'!9 +!3'+string(rms_val,format='(F5.2)')

oplot,[0,Nx],[0,0],linestyle=2
;***********
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
print,P,err_P


T=0
ISM=[-0.35,-0.14]/100.
ISM=[0,0]
avg_Q(*)=stoks(*,1,T)-Q_null-ISM(0)
;avg_U=LOWESS(findgen(Nx),U(*,T),box,2,2)
avg_U(*)=stoks(*,2,T)-U_null-ISM(1)
avg_P=sqrt(avg_Q^2+avg_U^2)
avg_FI=fltarr(Nx)
for x=0,Nx-1 do  avg_FI(x)=calc_atan(avg_Q(x),avg_U(x))/2; & R=where(avg_FI gt 178) & avg_FI(R)=avg_FI(R)-180

Window,3,xsize=600,ysize=900,xpos=600,ypos=0
!P.background=2^16-1 & !P.color=0
!P.multi=[0,1,5]
!P.charsize=2
R=where(wave gt 5000 and wave lt 6000)
RR=where(wave gt 4500 and wave lt 8000)
robomean,avg_P(R)*100,1,0.5,P,err_P
robomean,avg_FI(R) ,1,0.5,F,err_F
plot,wave(RR),FLUX(RR,T),xst=1,title=name(T)+'   PAslit='+string(PA(T),format='(F7.1)')
plot,wave(RR),avg_Q(RR)*100,xst=1,yrange=[-5,5],yst=1
plot,wave(RR),avg_U(RR)*100,xst=1,yrange=[-5,5],yst=1
plot,wave(RR),avg_P(RR)*100,xst=1,yrange=[0,10],yst=1,$
	title='P(V)='+string(P,format='(F5.2)')+'!9 +!3'+string(err_P,format='(F5.2)')
plot,wave(RR),avg_FI(RR),xst=1,title='!9P!3(V)='+string(F,format='(F7.1)')+'!9 +!3'+string(err_F,format='(F4.1)'),$
subtitle=wdir
print,P,err_P
END
window,4
!P.multi=[0,1,1]
;avg_vector=robust_estimation(avg_P,wx=10,tresh=1 )
avg_P=median(avg_P,10)
plot,wave(RR),avg_P(RR),xst=1
oplot,wave(RR),LOWESS(wave(RR),avg_P(RR),50,3,3),color=1e5,thick=2

end