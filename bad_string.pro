;correction _bad string
dir='h:\red_data.pol\AGN\E1841+643_131105\'
suff='flat'
ima=readfits(dir+'obj_i.fts',h)
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'Naxis2')
Nexp=sxpar(h,'NAXIS4')
x=findgen(Nx) & y=findgen(Ny)
kx=1600  & ex=0 & cu=0
plot,y,ima(kx,*,0,ex,cu),xst=1,xrange=[60,80],psym=10
badpos=70
oplot,[1,1]*badpos,[-1e5,1e7],linestyle=2
for cu=0,2 do begin
for ex=0,Nexp-1 do begin
for kx=1,Nx-2 do begin
ima(kx,badpos,3,ex,cu)=(ima(kx-1,badpos,3,ex,cu)+ima(kx+1,badpos,3,ex,cu))/2
endfor & endfor & endfor
writefits,dir+'obj_i.fts',ima,h
end