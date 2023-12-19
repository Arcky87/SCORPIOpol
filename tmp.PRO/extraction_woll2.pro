;test extraction_WOLL2
function extraction_WOLL2,image,BROAD=broad,POS=pos,PLOT=plot,GAUSS=gauss
a=size(image)
Nx=a(1)  & Ny=a(2)
window,2
x=Nx*0.7
contrast=0.8
y=findgen(Ny)  & wy=50
narrow=fltarr(Nx)
broad=fltarr(Nx)
for x=0,Nx-1 do begin
Vy=image(x,*)
;определение  ширины узкой компоненты
;narrow=gaussfit(y,Vy,GN,Nterms=4)
Vn= MPFITPEAK(y, Vy, AN,/moffat)
;pos=GN(1)
Vb=MPFITPEAK(y,smooth(Vy-Vn*contrast,10),AB,/gauss)
;print,Gb
plot,Vy,xst=1,psym=10
oplot,Vb,color=1e5
oplot,[0,Ny-1],[0,0],color=1E5
Vn=MPFITPEAk(y,Vy-Vb,AN,/moffat)
 broad(x)=total(Vb)
narrow(x)=total(Vn)
oplot,narrow,color=1e5
oplot,narrow+broad
oplot,Vy-narrow-broad,psym=10
endfor
window,3
!P.multi=[0,1,2]
plot,narrow,xst=1
plot,broad,xst=1

;Result=MULTIGAUS ( y, Vy, [Gn(1),Gb(1)], FWHM=[Gn(2),Gb(2)]*2.345,/plot,/fixFWHM)
; [,/FIXpos][,/FIXfwhm], [,/PLOT=PLOT])
res=fltarr(Nx)
return,res
end
dir='h:\red_data.pol\NGC3227_140325\'
dir='h:\red_data.pol\NGC5548_140325\'
dir='h:\red_data.pol\NGC4051_140324\'
dir='h:\red_data.pol\E1841+643_140324\'
cube=readfits(dir+'obj-sky.fts',h)
res=extraction_WOLL2(cube(*,*,0,0),/plot)

end