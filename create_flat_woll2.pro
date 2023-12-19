pro create_flat_WOLL2,dir,PLOT=plot

print,'create flat'

flat=readfits(dir+'flat_lin.fts',h,/silent)
grating='grating '+sxpar(h,'DISPERSE')
a=size(flat) & Nx=a(1)  & Ny=a(2)  & Npol=a(3)

avg_flat=flat   		;create array - average flat
map=fltarr(Nx,Ny) 		; ?
norm_X=fltarr(Nx,Npol) 	;X normalization
x=findgen(Nx) & y=findgen(Ny)
ratio=fltarr(Nx,Ny,2)

;;+++++++++++++++++++++++++++++++
;N=Npol
;normal=fltarr(N)
;
;for k=0,N-1 do begin
;	robomean,flat_cube(Nx/6.0:5.0*Nx/6.0,Ny/6.0:5.0*Ny/6.0,k),3,0.5,avgg,rms
;	normal(k)=avgg
;endfor
;
;normal=normal/total(normal)*N
;
;ima=fltarr(Nx,Ny)
;for kx=0,Nx-1 do begin
;	for ky=0,Ny-1 do begin
;		ima(kx,ky)=median(flat_cube(kx,ky,*)/normal)
;	endfor
;endfor
;
;slices=fltarr(Nx)
;for k=0, Nx-1 do begin
;	slices(k)=total(ima(k,*))/Ny
;	ima(k,*)=ima(k,*)/slices(k)
;endfor
;;+++++++++++++++++++++++++++++++

ratio(*,*,0)=flat(*,*,1)/flat(*,*,0) & ratio(*,*,1)=flat(*,*,3)/flat(*,*,2)

;формирование нормировки плоского поля
;нормировка по X
w=5
;norm_X=total(total(flat,3),2)/Npol/Ny
cgdisplay, wid=10
!p.multi=[0,1,4]
for k=0,Npol-1 do begin
	norm_X(*,k)=total(flat(*,Ny/2-w:Ny/2+w,k),2)/(2*w+1)
	cgplot,norm_X(*,k),color='blue'
	;norm_X(*,k)=lowess(indgen(Nx),norm_X(*,k),Nx/80,2,2)  ;!!! 60
	cgoplot,norm_X(*,k)
	for j=0,Ny-1 do begin
		flat(*,j,k)=flat(*,j,k)/norm_X(*,k)
	endfor
endfor
;
;;;нормировка по Y
; wx=200  & wy=10  & norm_Y=fltarr(Npol)
;for k=0,Npol-1 do norm_Y(k)=total(flat(Nx/2-wx:Nx/2+wx,Ny/2-wy:Ny/2+wy,k))/(2*wx+1)/(2*wy+1)
;;print,norm_Y
rat='ratio '+['  90/0','135/45']
;;;
;		for k=0,Npol-1 do begin
;	flat(*,*,k)=flat(*,*,k)/norm_Y(k)
;		ENDFOR
;;;
;;;сглаживание вдоль щели
;	avg_flat=flat
;			for k=0,Npol-1 do begin
;	for kx=0,Nx-1 do begin
;ff=goodpoly(y,avg_flat(kx,*,k),2,3,Yfit)
;avg_flat(kx,*,k)=Yfit
;	endfor
;	endfor
;
;;сглаживание вдоль дисперсии
;for k=0,Npol-1 do begin
;	for ky=0,Ny-1 do avg_flat(*,ky,k)=smooth(avg_flat(*,ky,k),Ny/2);,/edge_truncate)
;	;for ky=0,Ny-1 do avg_flat(*,ky,k)=lowess(indgen(Nx),avg_flat(*,ky,k),Nx/2,1,1);smooth(avg_flat(*,ky,k),Ny/4,/edge_truncate)
;endfor

ratio(*,*,0)=avg_flat(*,*,1)/avg_flat(*,*,0) & ratio(*,*,1)=avg_flat(*,*,3)/avg_flat(*,*,2) ; added for tests IY
;ratio(*,*,0)=flat(*,*,1)/flat(*,*,0) & ratio(*,*,1)=flat(*,*,3)/flat(*,*,2)


	if keyword_set(plot) then begin
cgdisplay,wid=0,xsize=650,ysize=900,title=dir+'  '+grating
ang='angle='+['  0',' 90',' 45','135']+' deg'
!p.multi=[0,1,6]
;cgplot,[0,Nx],[0,6],/nodata,/norm,xst=1,yst=1
		for k=0,Npol-1 do begin
			map(*,*)=avg_flat(*,*,k)
			;map(*,*)=avg_flat(*,*,k)
			;cgimage,255-bytscl(congrid(map,650,150),0.5,1.5),0,150*k
			;cgplot,[0,Nx],[0,6],/nodata,/norm,xst=1,yst=1
			cgimage,(congrid(map,650,150)), minvalue=0.2, XRange=[0,Nx], YRange=[0.5,1.5]
			cgoplot,findgen(Nx),map(*,Ny/2),color='orange'
			cgoplot,[0,Nx],[1,1]*0.5+ k,linestyle=2;,color=3e5
			xyouts,240,1.3,ang(k),align=0.5,charsize=2,CHARTHICK=2
		endfor
		for j=0,1 do begin
			;cgplot,[0,Nx],[0,6],/nodata,/norm,xst=1,yst=1
			cgimage,(congrid(ratio(*,*,j),650,150,1)), minvalue=0.2, XRange=[0,Nx], YRange=[0,2]
			cgoplot,findgen(Nx),ratio(*,Ny/2,j), color='orange'
			;oplot,[0,Nx],[1,1]*0.5+j+4,linestyle=2,color=3e5
			xyouts,240,1.7,rat(j),align=0.5,charsize=2,CHARTHICK=2
		endfor
	endif
writefits,dir+'tmp_flat.fts',flat,h
writefits,dir+'avg_flat.fts',avg_flat,h
writefits,dir+'avg_flat.fts',flat,h
writefits,dir+'avg_ratio.fts',ratio,h
end

log_dir='/data6/SCORPIO/sppol_pipeline_v2023.8/'
LOGFILE=DIALOG_PICKFILE(/read,path=log_dir+'LOGS/',FILTER='*.txt')
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'/')
wdir=log_dir+wdir(N_elements(wdir)-2)+'/'
print,wdir
create_flat_WOLL2,wdir,/plot
end
