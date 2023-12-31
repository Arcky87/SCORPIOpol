
pro kepler_motion,VEL,FI,PLOT=plot,PATH=path


END
wdir='h:\red_data.pol\AGN\PG0844+349_141121\'
cube_stoks=readfits(wdir+'stoks.fit',h)
Nx=sxpar(h,'NAXIS1') & Npol=sxpar(h,'NAXIS2') & Nexp=fltarr(3)
for k=0,2 do Nexp(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))  & print,Nexp
z=sxpar(h,'Z')  & name=sxpar(h,'NAME1')  & PA=sxpar(h,'PA1')
wave=findgen(Nx)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
xrng=[wave(0),wave(Nx-1)]
;xrng=[6500,7500]
xrng=[4500,7500]


Fobj=fltarr(Nx)  & Qobj=Fobj  & Uobj=Fobj &  Fstar=Fobj & Qstar=Fobj & Ustar=Fobj
Fnorm=fltarr(Nexp(0)) & Qnorm=Fnorm & Unorm=Fnorm
;����������� �������� ������������
for k=0,Nexp(0)-1 do  begin
Fnorm(k)=total(cube_stoks(Nx/2-Nx/4:Nx/2+Nx/4,0,k,0))/(Nx/2+1)
Qnorm(k)=total(cube_stoks(Nx/2-Nx/4:Nx/2+Nx/4,1,k,0))/(Nx/2+1)
Unorm(k)=total(cube_stoks(Nx/2-Nx/4:Nx/2+Nx/4,2,k,0))/(Nx/2+1)
endfor
Fnorm=Fnorm/median(Fnorm)
Qnorm=Qnorm/median(Qnorm)
Unorm=Unorm/median(Unorm)
;����������� �������� ������������ � ��������� ����������� ������
for x=0,Nx-1 do begin
Fobj(x)=median(cube_stoks(x,0,*,0)/Fnorm)
Qobj(x)=median(cube_stoks(x,1,*,0)/Qnorm)
Uobj(x)=median(cube_stoks(x,2,*,0)/Unorm)
Fstar(x)=median(cube_stoks(x,0,0:Nexp(1)-1,1))
Qstar(x)=median(cube_stoks(x,1,0:Nexp(1)-1,1))
Ustar(x)=median(cube_stoks(x,2,0:Nexp(1)-1,1))
endfor
;����������� ���������
atm_abs=total(readfits(wdir+'atm_abs.fit'),2)/4 & atm_abs=shift_s(atm_abs,1.5)
x=findgen(Nx)
;����������� ���������������� �����������
;�� ��������� ������� �����������
Qbias=LOWESS(x,Qstar,Nx/2,3,3) & Qstar=Qstar-Qbias
Ubias=LOWESS(x,Ustar,Nx/2,3,3) & Ustar=Ustar-Ubias
Qobj=Qobj-Qbias  & Uobj=Uobj-Ubias
;����������� ������������ ���������������
print,'zero polarization star ',sxpar(h,'NAME2')
tab_dir='d:\standards\data\'
star_name=['BD+28','BD+33','G191B']
 tab_name=['bd28d4211','bd33d4642','fg191b2b']+'.dat'
        R=where(STRUPCASE(strmid(sxpar(h,'NAME2'),0,5)) eq star_name,ind_sent)
CASE STRING(ind_sent,FORMAT='(I1)') OF
'0':	begin
		;�������� ������� ������ ������������ ����������������
		sent=1
		ind=FILE_TEST(wdir+'sent.fts')
		if ind eq 1 then sent=readfits(wdir+'sent.fts') ELSE $
		print,'no table spectrophotometric standard!!
		end
'1': 	begin
		star_obs=Fstar

		;window,3
		;!P.multi=[0,1,3]
		;plot,wave,star_obs,xst=1
		star_obs=create_continuum(star_obs)
		;oplot,wave,star_obs,color=1e5
		star_tab=read_st(tab_dir+tab_name(R),wave,/print)
		;plot,wave,star_tab,xst=1
		star_tab=create_continuum(star_tab)
		;oplot,wave,star_tab,thick=2
		sent=star_obs/star_tab & sent=sent/sent(Nx/2)
		;plot,wave,sent,xst=1;,yrange=[0,10]
		sent=create_continuum(sent)
		;oplot,wave,sent
		writefits,wdir+'sent.fts',sent
		hs=headfits(wdir+'sent.fts')
		sxaddpar,hs,'STAR',sxpar(h,'NAME2')
		sxaddpar,hs,'CRVAL1',sxpar(h,'CRVAL1')
		sxaddpar,hs,'CDELT1',sxpar(h,'CDELT1')
		modfits,wdir+'sent.fts',0,hs
		end
ENDCASE
Fobj=Fobj/atm_abs/sent
;Fcnt=create_continuum(Fobj,SIGN=1)

;��������� ������� ������� � �������
cut=0.75
ff=goodpoly(x(0:Nx*cut),Qobj(0:Nx*cut),1,3,fit)
fit=ff(0)+ff(1)*x & Qobj=Qobj-LOWESS(x,Qobj,Nx/4,3,3)+fit
ff=goodpoly(x(0:Nx*cut),Uobj(0:Nx*cut),1,3,fit)
fit=ff(0)+ff(1)*x & Uobj=Uobj-LOWESS(x,Uobj,Nx/4,3,3)+fit
;���������� ���������� (Q,U) �������������� �������� � �������� �����
PA_null=317.3
TETA=2*(PA-PA_null)*!PI/180
Q=Qobj*cos(TETA)-Uobj*sin(TETA)
U=Qobj*sin(TETA)+Uobj*cos(TETA)
Qobj=Q & Uobj=U
S=10 & Qobj=median(Qobj,5) & Uobj=median(Uobj,5)
window,2,xsize=500,ysize=1000
!P.multi=[0,1,3]
plot,wave,Fobj,xst=1,xrange=xrng
;oplot,wave,Fcnt
oplot,[1,1]*6562.8*(1+z),[0,1e8],linestyle=2,psym=10
plot,wave,Qobj,xst=1,xrange=xrng,yrange=[-1,1]*0.02,psym=10
robomean,Qobj-LOWESS(x,Qobj,20,3,3),3,0.5,val,rms  & print,val,rms
plot,wave,Uobj,xst=1,xrange=xrng,yrange=[-1,1]*0.02,psym=10

bobo=1
end
