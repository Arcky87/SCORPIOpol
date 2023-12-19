;COMPARE  woll_2
FLAT_old=READFITS('H:\red_data.pol\AGN\Mkn382_131104\flat_lin.fts',h)
FLAT_new=READFITS('H:\red_data.pol\AGN\Mkn382_160408\flat_lin.fts',h)
print,size(flat_new)

flat=flat_old

Nx=1851 & Ny=80
window,2,xsize=200+Nx/2,ysize=4*Ny*2,title='WOLLASTON-2 old'
tv,255-bytarr(200+Nx/2,8*Ny)
for k=0,3 do begin



tv,255-bytscl(congrid(flat(*,*,k),Nx/2,2*Ny,1),0,3e5),200,2*Ny*k
vector=total(flat(Nx/2-Nx/4:Nx/2+Nx/4,*,k),1)/(Nx/2-1)
plot,vector,findgen(Ny),xrange=[0,1.5e5],xst=1,xcharsize=1e-5,$
position=[0,Ny*2*k,200,Ny*2*(k+1)],/device,/noerase,color=1
endfor
end
