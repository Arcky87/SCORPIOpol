;
function correction_image,ima,C,Yc=yc,Hslit=hslit,PLOT=plot,TITLE=title
;исправление искажение спектра сравнения
a=size(ima) & Nx=a(1) & Ny=a(2)
ima_corr=WARP_TRI(c(*,0),c(*,1),c(*,2),c(*,3),ima)
ima_corr=ima_corr(*,Yc-Hslit/2:Yc+Hslit/2-1)
if keyword_set(plot) then begin
Window,2,xsize=Nx/2,ysize=Ny+Hslit,TITLE=title
map=congrid(ima,Nx/2,Ny) & robomean,map,3,0.5,avg_map,rms_map
co=50
tv,255-bytscl(map,avg_map-rms_map*co/3,avg_map+rms_map*co)
map=congrid(ima_corr,Nx/2,Hslit)
tv,255-bytscl(map,avg_map-rms_map*co/3,avg_map+rms_map*co),0,Ny+1
endif
return,ima_corr
end
