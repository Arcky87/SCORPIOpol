function extraction_spectra,ima,BROAD=broad,POS=pos,WIDTH=width,PLOT=plot
a=size(ima)
Ndeg=1
if not(keyword_set(broad)) then broad=0
if not(keyword_set(pos)) then pos=[a(2)/2]
if not(keyword_set(width)) then width=20

N_narrow=N_elements(pos)
N_comp=N_narrow+1
spectra=fltarr(a(1),N_comp)
y=findgen(a(2))
x=findgen(a(1))
Nx=a(1) & Ny=a(2) & wx=100
print,Nx,Ny
Vy=total(ima(Nx/2-wx:Nx/2+wx,*),1)

window,2
plot,Vy,xst=1
res=MULTIGAUS ( y,Vy ,[1]*Ny/2, FWHM=[5],/plot,YFIT=YFIT)
;oplot,Yfit
;goto,fin
;remove broad component
	if keyword_set(broad) then begin
FOR kx=0,a(1)-1 DO BEGIN
	;убирание компонент
	slice=ima(kx,*)
	for i=0,N_comp-2 do begin
		gau=gaussfit(y(pos(i)-width/2:pos(i)+width/2),slice(pos(i)-width/2:pos(i)+width/2),G)
		gau=gauss_profile(y(pos(i)-width/2:pos(i)+width/2),G(0:2))
		slice(pos(i)-width/2:pos(i)+width/2)=slice(pos(i)-width/2:pos(i)+width/2)-gau
	endfor
	;выделение широкой компоненты
	gau_broad=gaussfit(y,slice,G)
	ima(kx,*)=ima(kx,*)-gau_broad
	spectra(kx,N_comp-1)=total(gau_broad)
ENDFOR
	endif




;декомпозиция узкой компоненты
	;построение траекторий компонент
Npos=50 &  w=a(1)/Npos
xpos=findgen(Npos)*w+w/2
ypos=fltarr(Npos,N_narrow)
wpos=fltarr(Npos,N_narrow)
	for j=0,Npos-1 do begin
	vector=total(ima(xpos(j)-w/2:xpos(j)+w/2,*),1)
	res=MULTIGAUS ( findgen(a(2)),vector ,POS, FWHM=fltarr(N_narrow)+5); [,/FIXpos][,/FIXfwhm], [,/PLOT=PLOT])
;                          [,/SILENT] [ABSORP=ABSORP] [,/LIM_FWHM=], [,/DOUBL_LEN]
;                          [,/DOUBL_RATIO=] [,/FIXRATIO=],YFIT=YFIT,sigma=sigma)
ypos(j,*)=res.center
wpos(j,*)=res.FWHM
	endfor
;аппроксимация траекторий спектров узких компонентов их ширин
start_pos=10
	y_tra=fltarr(a(1),N_narrow)
	w_tra=fltarr(a(1),N_narrow)
	for j=0,N_narrow-1 do begin
		f=goodpoly(xpos(start_pos:Npos-1),ypos(start_pos:Npos-1,j),Ndeg,1)
		fit=0 & for i=0,Ndeg do  fit=fit+f(i)*x^i
		y_tra(*,j)=fit
		f=goodpoly(xpos(start_pos:Npos-1),wpos(start_pos:Npos-1,j),Ndeg,1)
		fit=0 & for i=0,Ndeg do  fit=fit+f(i)*x^i
		w_tra(*,j)=fit
	endfor
if keyword_set(plot) then begin
	window,2,xsize=600,ysize= 400
	!P.multi=[0,1,N_narrow]
	plot,[0,a(1)],[0,a(2)],/nodata,xst=1,yst=1
	for k=0,N_narrow-1 do begin
		oplot,xpos,ypos(*,k),psym=6,symsize=0.25
		oplot,x,y_tra(*,k)
			endfor
	plot,[0,a(1)],[0,10],/nodata,xst=1,yst=1
	for k=0,N_narrow-1 do begin
		oplot,xpos,wpos(*,k),psym=6,symsize=0.25
		oplot,x,w_tra(*,k)
			endfor
endif
	;аппрoксимация компонент гауссианами

;for kx=0,a(1)-1 do begin
;res=MULTIGAUS ( y,ima(kx,*) ,y_tra(kx,*), FWHM=w_tra(kx,*),/FIXpos,/FIXfwhm)
;spectra(kx,0:N_narrow-1)=res.flux(*)
;endfor
;if keyword_set(plot) then begin
;window,3,xsize=600,ysize= 400
;!P.multi=[0,1,N_comp]
;titl=strarr(N_comp)+'narrow'
;titl(N_comp-1)='broad'
;for j=0,N_comp-1 do plot,spectra(*,j),xst=1,title='component    '+titl(j)+string(J)
;!P.multi=[0,1,1]
;endif
;wait,1
fin:
return,spectra
end

;LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path='h:\red_data.pol\URAN\LOGS\')
wdir=DEF_WDIR(LOGFILE)

cube=readfits(wdir+'obj-sky.fts',h)
a=size(cube) & ima=total(cube(*,0:a(2)-1,*,0,0),3)
;spectra=extraction_spectra(ima,pos=[80],width=20,/plot)
spectra=extraction_spectra(ima,pos=[40],width=40,/plot)

cgdisplay, wid=10
cgplot, spectra

end