DIR='h:\red_data.pol\AGN\NGC7469_151210\'
cube=readfits(dir+'obj-sky.fts',h)
print,size(cube)
Nx=2201  & Ny=320 & Npol=4 & Nexp=20
map=fltarr(Nx,Ny)
window,2,xsize=Nx/2,ysize=Ny
map(*,*)=cube(*,*,0,0,0)
w=100
tv,bytscl(congrid(map,Nx/2,Ny),0,max(map)/10)
vector_tot=total(cube(Nx/2-w:Nx/2+w,*,*,*,0),4)/Nexp
vector_tot=total(vector_tot,3)/4
vector_tot=total(vector_tot,1)/(2*w+1)
window,3
!P.multi=[0,1,1]
plot,vector_tot,xst=1
end