;
function correction_image,ima,C,Yc=yc,Hslit=hslit,PLOT=plot,TITLE=title,TRIM_EDGE=trim_edge
;исправление искажение спектра сравнения
a=size(ima) & Nx=a(1) & Ny=a(2)
ima_corr=WARP_TRI(c(*,0),c(*,1),c(*,2),c(*,3),ima)
;ima_corr=ima_corr(*,Yc-Hslit/2:Yc+Hslit/2-1) ; commented for checking purposes. Return in deploy
if not keyword_set(trim_edge) or (N_ELEMENTS(trim_edge) ne 2) then trim_edge = [0, 0]
ima_corr=ima_corr(trim_edge[0]:Nx-trim_edge[1]-1,Yc-Hslit/2:Yc+Hslit/2-1)
if keyword_set(plot) then begin
;Window,2,xsize=Nx/2,ysize=Ny+Ny,TITLE=title ;Hslit
cgDisplay, Nx/2, Ny+Ny, Title=title
!p.multi = [0, 1, 2] 
co=50
map=congrid(ima,Nx/2,Ny) & robomean,map,3,0.5,avg_map,rms_map
cgImage, 255-bytscl(map,avg_map-rms_map*co/3,avg_map+rms_map*co), /SCALE, /KEEP_ASPECT, MARGIN=0.1,TITLE='Original Image'
map=congrid(ima_corr,Nx/2,Ny) ;Hslit
cgImage, 255-bytscl(map,avg_map-rms_map*co/3,avg_map+rms_map*co), /SCALE, /KEEP_ASPECT, MARGIN=0.1,TITLE='Corrected Image'
!p.multi = 0
;co=50
;tv,255-bytscl(map,avg_map-rms_map*co/3,avg_map+rms_map*co)
;map=congrid(ima_corr,Nx/2,Ny) ;Hslit
;tv,255-bytscl(map,avg_map-rms_map*co/3,avg_map+rms_map*co),0,Ny+1
endif
return,ima_corr
end
