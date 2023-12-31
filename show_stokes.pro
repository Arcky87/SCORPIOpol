
;
;LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path='h:\red_data.pol\LOGS\')
LOGFILE='h:\red_data.pol\LOGS\Mrk1502_141122.txt'
!P.charsize=1
wdir=sxpar(read_table(LOGFILE),'w_dir')
stoks=readfits(wdir+'stoks.fit',h)
 flux=readfits(wdir+'flux.fit')
;����������� ����-������
Nx=sxpar(h,'NAXIS1') & x=findgen(Nx)
Nexp=intarr(3) & Ncube=sxpar(h,'NAXIS4')
for k=0,2 do Nexp(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))

Q_null=LOWESS(x,total(stoks(*,1,0:Nexp(1)-1,1),3)/Nexp(1),Nx/2,2,2)
Q_null=LOWESS(x,Q_null,Nx/2,2,2)
U_null=LOWESS(x,total(stoks(*,2,0:Nexp(1)-1,1),3)/Nexp(1),Nx/2,2,2)
U_null=LOWESS(x,U_null,Nx/2,2,2)
for i=0,Ncube-1 do begin
for j=0,Nexp(i)-1 do begin
stoks(*,1,j,i)=stoks(*,1,j,i)-Q_null
stoks(*,2,j,i)=stoks(*,2,j,i)-U_null
endfor & endfor

if FILE_TEST(wdir) eq 0 then sent=readfits(wdir+'sent.fts') else sent=1
;sent=readfits(wdir+'sent.fts'

lambda=findgen(Nx)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
date=sxpar(h,'DATE-OBS')
Texp=sxpar(h,'EXPTIME'+string(j+1,format='(I1)'))
juldate,str_sep(date,'-'),JD

print,date,JD
Nx=sxpar(h,'NAXIS1')
PA=sxpar(h,'PA'+string(j+1,format='(I1)'))
print,PA
name_obj=sxpar(h,'Name'+string(j+1,format='(I1)'))
print,'PA=',PA
PA=PA*!PI/180

;**************************************************************************
wx=2 &  xrng=[4500,7500] & tresh=1	& j=1  & xrng=[4500,7500]  & dP=5 & err=0
Wave_c=5500  & width_wave=500

;**************************************-***********************************
print,Nx
Npos=FIX(float(Nx)/FLOAT(wx/2))-2
xpos=findgen(Npos)*(wx/2)+wx/2
wave=lambda(xpos)
R=where(wave gt xrng(0) and wave lt xrng(1))
Rc=where(wave gt wave_c-width_wave/2 and wave lt wave_c+width_wave/2)
; ���������� ��������� ������ ���������� ������ � ����� �� ���� �����������
F=fltarr(Npos,Ncube)
avg_I=F & avg_Q=F & avg_U=F & avg_P=F & avg_FI=F & & rms_I=F & rms_Q=F  & rms_U=F  &  rms_P=F & rms_FI=F
j=0
for k=0,Npos-1 do begin
;average intensity
robomean,stoks(xpos(k)-wx/2:xpos(k) +wx/2,0,*,j),tresh,0.5,avg_val,rms_val
avg_I(k)=avg_val & rms_I(k)=rms_val
robomean,100*stoks(xpos(k)-wx/2:xpos(k) +wx/2,1,*,j),tresh,0.5,avg_val,rms_val
avg_Q(k)=avg_val  & rms_Q(k)=rms_val
robomean,100*stoks(xpos(k)-wx/2:xpos(k) +wx/2,2,*,j),tresh,0.5,avg_val,rms_val
avg_U(k)=avg_val & rms_U(k)=rms_val
endfor

avg_P=sqrt(avg_Q^2+avg_U^2) & rms_P=sqrt(rms_Q^2+rms_U^2)/sqrt(2)
avg_FI=calc_atan(avg_Q,avg_U)/2  & rms_FI=rms_P*28.2

; ���������� ��������� ������ � ������� ������
print,RC
robomean,avg_Q(Rc),TRESH,0.5,mean_Q,err_Q
print,mean_Q,err_Q
robomean,avg_U(Rc),TRESH,0.5,mean_U,err_U
print,mean_U,err_U
robomean,avg_P(Rc),TRESH,0.5,mean_P,err_P
print,mean_P,err_P
robomean,avg_FI(Rc),TRESH,0.5,mean_FI,err_FI
print,mean_FI,err_FI
;**********************************
dH=0.02
window,2,xsize=800,ysize=1000
!P.multi=[0,1,1]
;intensity
order=FIX(ALOG10(max(Flux(Nx*0.2:Nx-1))))
print,order
plot,lambda,flux/10.^order,xrange=xrng,xst=1,thick=1,ytickinterval=1,$
	position=[0.13,.80+dH,.99,.95+dH],/norm,xcharsize=1e-5,$
    	;title=STRCOMPRESS(name(j)+', '+date+', Texp='+$
		 ; string(exptime(j),format='(I4)')+' s, PA!Dslit!N= '+string(PA(j),format='(F6.1)')+' deg'),$
	ytitle='FLux, 10!U'+string(order,format='(I1)')+'!N, ADU/px'



;,yrange=[-1,1]*0.05,yst=1,thick=2
;plot,wave,avg_I*avg_P,xrange=xrng,xst=1;,yrange=[-1,1]*0.05,yst=1
;Q-Stokes
robomean,avg_Q(R),3,0.5,avg_val,rms_val
plot,wave,avg_Q,xrange=xrng,xst=1,yrange=avg_val+[-1,1]*dP ,yst=1,psym=10,thick=err+1	,$
	position=[0.13,.50+dH,.99,.65+dH],/norm,/noerase,xcharsize=1e-5,ytitle='Q-Stoks, %'
if err eq 1 then   oploterr,wave,avg_Q,rms_Q

xyouts,(xrng(0)+xrng(1))/2,avg_val+dp*0.75,'Q('+string(Wave_c,format='(I4)')+') ='+$
	string(mean_Q,format='(F6.2)')+'!9 +!3'+string(err_Q,format='(F5.2)')+' %',align=0.5
goto,fin
;U-stokes
robomean,avg_U(R),3,0.5,avg_val,rms_val
plot,wave,avg_U,xrange=xrng,xst=1,yrange=avg_val+[-1,1]*dP ,yst=1,psym=10,thick=err+1
if err eq 1 then   oploterr,wave,avg_U,rms_U,psym=10
xyouts,(xrng(0)+xrng(1))/2,avg_val+dp*0.75,'U('+string(Wave_c,format='(I4)')+') ='+$
	string(mean_U,format='(F6.2)')+'!9 +!3'+string(err_U,format='(F5.2)')+' %',align=0.5

; P  degree of linear polarization

plot,wave,avg_P,xrange=xrng,xst=1,yrange=[0,1]*dP  ,yst=1,psym=10,thick=err+1
if err eq 1 then   oploterr,wave,avg_P,rms_P,psym=10
xyouts,(xrng(0)+xrng(1))/2, dp*0.875,'P('+string(Wave_c,format='(I4)')+') ='+$
	string(mean_P,format='(F6.2)')+'!9 +!3'+string(err_P,format='(F5.2)')+' %',align=0.5

; angle polarization plane
robomean,avg_FI(R),3,0.5,avg_val,rms_val
plot,wave,avg_FI,xrange=xrng,xst=1,yrange=avg_val+[-1,1]*90 ,yst=1,psym=10,thick=err+1
if err eq 1 then   oploterr,wave,avg_FI,rms_FI,psym=10
xyouts,(xrng(0)+xrng(1))/2,avg_val+90*0.75,'FI('+string(Wave_c,format='(I4)')+') ='+$
	string(mean_FI,format='(F6.2)')+'!9 +!3'+string(err_FI,format='(F5.2)')+' DEG',align=0.5




fin:



END
set_plot,'PS'
device,file=wdir+'res.ps',/portrait,xsize=19,ysize=28,xoffset=1,yoffset=1
!P.multi=[0,1,6]
!P.charsize=2
w=10
tresh=2

lim=[-1.5,1.5]
xrng=[4500,lambda(Nx-1)]
xrng=[4500,7280]
flux=stoks(*,0,j)/sent
Rc=where(ABS(lambda-6000) lt sxpar(h,'CDELT1'))
;print,lambda(Rc(0))
flux=flux/flux(Rc(0))
flux_P=flux*sqrt(stoks(*,1,j)^2+stoks(*,2,j)^2)*100

plot,lambda,flux,xst=1,xrange=xrng,title=name_obj+' '+date+', JD'+$
	string(JD,format='(I5)')+', Texp='+string(Texp,format='(I4)')+' s, PA_slit='+$
    string(sxpar(h,'PA'+string(j+1,format='(I1)')),format='(F6.1)')+' deg',xtitle='Wavelength, A',ytitle='Flux'
wind=LOWESS(findgen(Nx),stoks(*,0,j),Nx/2,1)
wind=wind/max(wind)
R=where(wind le 0.5) & wind(R)=1
wind=LOWESS(findgen(Nx),wind,Nx/2,2)
wind=1/wind
;������ ������� �� ����
writefits,wdir+'spectra.fts',flux,h
;������������ ����������� ��������

Npos=(Nx-wind(0)*w/2 -wind(Nx-1)*w/2)/w-1

xpos=findgen(Npos)*w+wind(0)*w/2
print,xpos(0),xpos(Npos-1),Nx-1-w
wave=lambda(xpos)

wave_c=6000
w_wave=500
wave_min=wave_c-w_wave  & wave_max=wave_c+w_wave

R_min=where(ABS(wave-wave_min) lt w)
R_max=where(ABS(wave-wave_max) lt w)
REG=[R_min(0),R_max(0)]
print,'REG',REG
ampl=wind(xpos)
Q=fltarr(Npos) & U=fltarr(Npos) & rms_U=U  & rms_Q=Q
f_P=fltarr(Npos) & rms_f_P=f_P
for k=0,Npos-1 do begin
;Intensity
robomean,flux_P(xpos(k)-w/2*wind(xpos(k)):xpos(k)+w/2*wind(xpos(k))),tresh,0.5,avg_val,rms_val
f_P(k)=avg_val
rms_f_P(k)=rms_val

robomean,stoks(xpos(k)-w/2*wind(xpos(k)):xpos(k)+w/2*wind(xpos(k)),1,j),tresh,0.5,avg_val,rms_val
Q(k)=avg_val
rms_Q(k)=rms_val
robomean,stoks(xpos(k)-w/2*wind(xpos(k)):xpos(k)+w/2*wind(xpos(k)),2,j),tresh,0.5,avg_val,rms_val
U(k)=avg_val
rms_U(k)=rms_val
endfor

;������������� ���������� ��������������� �������
f_P_tmp=f_P
;f_P_tmp(0:29)=f_P_tmp(30)

for k=0,3 do begin
cont=LOWESS(findgen(Npos),f_P_tmp,Npos/2,1)
robomean,f_P_tmp-cont,3,0.5,mean,rms
R=where(f_P_tmp-cont gt rms,ind)
if ind gt 1 then f_P_tmp(R)=cont(R)
endfor
;

;������������� ������� ����������
mul=2.35482
Res=MULTIGAUS ( wave,f_P-cont,[6900], FWHM=[300],YFIT=YFIT);,/plot)
broad=gaussian(wave,[res.max,res.center,res.FWHM/mul])
R=where(ABS(wave-res.center) lt w)
ratio=res.max/cont(R(0))
;print,SQRT(sigma.center)
V_broad=(res.center/6563-1)*3E5

Rb=where(wave gt 6600 and wave lt 7200)
gau=gaussfit(wave(Rb),f_P(Rb),G,Nterms=4,sigma=s)
;print,g(1),g(2)*mul
;print,s(1),s(2)*mul
;���������� ��������� �������� ��� 3�390.3
z=1.056
narrow=[4861,4959,5007,6563]
;find lines peaks
wl=3
dL=lambda(1)-lambda(0)
for j=0,3 do begin
R=where(ABS(lambda-narrow(j)*z) lt dL)
pos=R(0)
parab=goodpoly(lambda(pos-wL:pos+wL),flux(pos-wL:pos+wL),2,3,Yfit)
pos=-parab(1)/parab(2)/2
narrow(j)=(pos/narrow(j)-1)*3E5
endfor
robomean,narrow,3,0.4,avg_Vsys,rms_Vsys
dV=V_broad-avg_Vsys
dVrms=s(1)/6563*3E5
dVrms=201
;print,dV,dVrms
FWHM=res.FWHM/6563*3E5  & FWHMrms=s(2)*mul/6563*3E5
;print,FWHM,FWHMrms


;���������
dV=dV;-300
plot,wave,f_P ,xst=1,xrange=xrng ,psym=10,ytitle='Flux x P(%)',xtitle='Wavelength, Angstrom';,$
	;title='(V-Vsys)='+string(dV,format='(I5)')+'!9+!3'+string(dVrms,format='(I3)')+$
	;' km/s, FWHM='+string(FWHM,format='(I5)')+' km/s, k='+$
	;string(ratio,format='(F5.2)'),yrange=[0,10]
	;oplot,w;e,cont
;	oplot,wave,(cont+broad)
	;oplot,wave(Rb)*gau,thick=2
PA=PA;-180
;PA=-PA

;u=-u
Q_new=Q*cos(2*PA)+U*sin(2*PA)
U_new=Q*sin(2*PA)-U*cos(2*PA)
;Q_new=Q
;U_new=U
;
;���������� ������� ��������
wc='('+string(wave_c,format='(I4)')+') ='
robomean,Q_new(REG(0):REG(1))*100,tresh,0.5,avg_Q ,rms_Q_new
robomean,U_new(REG(0):REG(1))*100,tresh,0.5,avg_U,rms_U_new

plot,wave,U_new*100,xst=1,xrange=xrng,yrange=lim+avg_U,psym=10,yst=1 ,$
	title='Q'+wc+string(avg_U,format='(F7.2)')+' !9+!3'+string(rms_U_new,format='(F5.2)')+' %',$
	xtitle='Wavelength, Angstrom',ytitle='Q-Stoks, %'
plot,wave,Q_new*100,xst=1,xrange=xrng,yrange=lim+avg_Q,psym=10,yst=1,$
	title='U'+wc+string(avg_Q,format='(F7.2)')+' !9+!3'+string(rms_Q_new,format='(F5.2)')+' %',$
	xtitle='Wavelength, Angstrom',ytitle='U-Stoks, %'
P=sqrt(U_new^2+Q_new^2)
robomean,P(REG(0):REG(1))*100,tresh,0.5,avg_P,rms_P
plot,wave,P*100,xst=1,xrange=xrng,yrange=lim+avg_P,psym=10,yst=1,$
	title='P'+wc+string(avg_P,format='(F7.2)')+' !9+!3'+string(rms_P,format='(F5.2)')+' %',$
	xtitle='Wavelength, Angstrom',ytitle='P, %'

FI=-atan(U_new/Q_new)*180/!PI/2+180
FI=calc_atan(u_new,q_new)/2;+180
;FI=-PA+FI+5.3;+180
robomean,FI(REG(0):REG(1)),tresh,0.5,avg_FI,rms_FI

R=where(FI lt 100,ind) & if ind gt 1 then FI(R)=FI(R)+180
plot,wave,FI,xst=1,xrange=xrng,psym=10,yrange=[-30,40]+avg_FI,$
	title='PA'+wc+string(avg_FI,format='(F7.1)')+' !9+!3'+string(rms_FI,format='(F5.1)')+' deg',$
	xtitle='Wavelength, Angstrom',ytitle='PA, deg'

xyouts,0,0,LOGFILE,/norm,charsize=0.75
device,/close
sxdelpar, h, ['NAME2','CUBE2','START2','EXPTIME2','RA2','DEC2','A2','Z2','PA2',$
              'NAME3','CUBE3','START3','EXPTIME3','RA3','DEC3','A3','Z3','PA3']
writefits,wdir+'flux.fit',flux,h
sxaddpar,h,'CRVAL1',wave(0)
sxaddpar,h,'CDELT1',wave(1)-wave(0)
out_stokes=fltarr(N_elements(wave),5)
out_stokes(*,0)=f_P
out_stokes(*,1)=U_new*100
out_stokes(*,2)=Q_new*100
out_stokes(*,3)=P*100
out_stokes(*,4)=FI

writefits,wdir+'stokes.fit',out_stokes,h
Window,2,xsize=400,ysize=750
!P.multi=[0,1,5]
tmp=readfits(wdir+'stokes.fit')

for j=0,4 do plot,tmp(*,j),xst=1
;set_plot,'WIN'
END
set_plot,'PS'
device,file=wdir+'H_alpha.ps',xsize=14,ysize=22.5,xoffset=4,yoffset=4
Lc=6685.

plot,lambda,Flux*10,xst=1,xrange=[6200,7200],position=[0.15,0.15,0.95,0.35],/norm,$
	 XTICKINTERVAL=500,yrange=[0.6,5.9]*11,yst=1,YTICKINTERVAL=20,Yminor=2,$
	 ytitle='Flux',xtitle='Wavelength (A)'
oplot,[1,1]*Lc,[0,1000],linestyle=2
plot,wave,P*100,xst=1,xrange=[6200,7200],position=[0.15,0.35,0.95,0.55],/norm,$
 XTICKINTERVAL=500,psym=10,xcharsize=1e-5,yrange=[0.1,1.9],yst=1,$
	 ytitle='Polarization(%)',YTICKINTERVAL=0.5
oplot,[1,1]*Lc,[0,1000],linestyle=2
plot,wave,f_P/10,xst=1,xrange=[6200,7200],position=[0.15,0.55,0.95,0.75],/norm,$
 XTICKINTERVAL=500,psym=10,xcharsize=1e-5,Yminor=5,yrange=[0.05,0.39],yst=1,$
	 ytitle='Polarized Flux',YTICKINTERVAL=0.1
oplot,[1,1]*Lc,[0,1000],linestyle=2
plot,wave,FI,xst=1,xrange=[6200,7200],position=[0.15,0.75,0.95,0.95],/norm,$
 XTICKINTERVAL=500,psym=10,xcharsize=1e-5,yrange=[10,210],yst=1,$
 YTICKINTERVAL=50,yminor=5, ytitle='!7h!3(degrees)',title=date+', JD'+	string(JD,format='(I5)')
oplot,[1,1]*Lc,[0,1000],linestyle=2
device,/close
set_plot,'WIN'
END
