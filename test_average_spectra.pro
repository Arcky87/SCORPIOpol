;test avg
wdir='h:\red_data.pol\Sy1\Mkn79_151207\'
stoks=readfits(wdir+'avg_stoks.fit',h)
sent=readfits(wdir+'sent.fts')
a=size(spectra)
Nx=a(1) & Np=a(2) & Nt=a(3)

window,0,xsize=700,ysize=1100
!P.multi=[0,1,5]
!P.charsize=2
t=2
I=stoks(*,0,t)/sent
Q=stoks(*,1,t)
U=stoks(*,2,t)
P=SQRT(Q^2+U^2)
FI=calc_atan(Q,U)/2
plot,I,xst=1
plot,Q*100,xst=1,yrange=[-2,2]*4
plot,U*100,xst=1,yrange=[-2,2]*4
plot,P*100,xst=1,yrange=[0,2]*4
plot,FI,xst=1;,yrange=[100,120]
;;вычисление параметров Стокса
;
;Q=(spectra(*,1,t)-spectra(*,0,t))/(spectra(*,1,t)+spectra(*,0,t))
;U=(spectra(*,3,t)-spectra(*,2,t))/(spectra(*,3,t)+spectra(*,2,t))
;P=sqrt(Q^2+U^2)
;FI=calc_atan(Q,U)/2
;window,2,xsize=900,ysize=1100
;plot,Q,xst=1
;plot,U,xst=1
;plot,P,xst=1
;plot,FI,xst=1
end