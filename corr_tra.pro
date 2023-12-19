;исправление остаточной кривизны траекторий
function corr_tra,ima,WIN=win,PLOT=plot,NDEG=Ndeg
if not(keyword_set(win)) then win=20
if not(keyword_set(Ndeg)) then Ndeg=3

wx=win/2  & wy=50
a=size(ima) & Nx=a(1) & Ny=a(2)
x=findgen(Nx) & y=findgen(Ny)
;формирование узлов
Npos=Nx/win
xpos=findgen(Npos)*win+wx   & ypos=fltarr(Npos)
for k=0,Npos-1 do begin
Vy=total(ima(xpos(k)-wx:xpos(k)+wx,*),1)
Vy=median(Vy,10)  ;чистка частиц
max_value=max(Vy,Nmax) & w=3

;Nmax=Ny/2
f=goodpoly(y(Nmax-10:Nmax+10),Vy(Nmax-10:Nmax+10),2,3)
ypos(k)=-f(1)/f(2)/2
print,xpos(k),Nmax,ypos(k)
endfor

;аппроксимация кривизны нитки спектра
f=goodpoly(xpos,ypos,Ndeg,3,Yfit)
dy=0 & for j=0,Ndeg do dy=dy+f(j)*x^j & dy=dy-Ny/2
WINDOW,0,xsize=1000,ysize=610
wy=Ny/4

region=ima(*,Ny/2-wy:Ny/2+wy) & robomean,region(Nx/2-20:Nx/2+20,*),3,0.5,avg_reg,rms_reg
region(*,wy)=fltarr(Nx)
tv,bytscl(congrid(region,950,200),avg_reg-3*rms_reg,avg_reg+10*rms_reg),40,400
plot,xpos,ypos-Ny/2,psym=6,xrange=[0,Nx],xst=1,yrange=[-wy,+wy],yst=1,$
	position=[40,200,990,400],/device,/noerase
oplot,x,dy & oplot,[0,Nx],[0,0],linestyle=2
;исправление кривизны
for kx=0,Nx-1 do ima(kx,*)=shift_s(ima(kx,*),-dy(kx))
region=ima(*,Ny/2-wy:Ny/2+wy) & robomean,region(Nx/2-20:Nx/2+10,*),3,0.5,avg_reg,rms_reg
region(*,wy)=fltarr(Nx)
tv,bytscl(congrid(region,950,200),avg_reg-3*rms_reg,avg_reg+5*rms_reg),40,0

;
	if keyword_set(plot) then begin

	endif
return,ima
end




wdir='h:\red_data.pol\Sy1\ngc3516_170131\'
wdir='h:\red_data.pol\Sy1\IRAS03450_141020\'
wdir='h:\red_data.pol\Sy1\NGC5548_150325\'
obj_cube=readfits(wdir+'obj-sky.fts',h)
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2') & Np=sxpar(h,'NAXIS3')
Nt=sxpar(h,'NAXIS5') & Nexp=intarr(Nt) & gain=sxpar(h,'GAIN')
for k=0,Nt-1 do Nexp(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))
ima=fltarr(Nx,Ny)
e=2 & p=2 & t=0

;ima(*,*)=obj_cube(*,*,p,e,t)
		for t=0,Nt-1 do begin
	for e=0,Nexp(t)-1 do begin
for p=0,Np-1 do begin
;ima=obj_cube(*,*,p,e,t)
obj_cube(*,*,p,e,t)=corr_tra(obj_cube(*,*,p,e,t),Ndeg=2,win=200);,Ndeg=2)
obj_cube(*,*,p,e,t)=corr_tra(obj_cube(*,*,p,e,t),/plot,Ndeg=2,win=200);,Ndeg=2)
;;bj_cube(*,*,p,e,t)=out
wait,0.25
endfor
	endfor
		endfor
writefits,wdir+'obj-sky.fts',obj_cube,h
end