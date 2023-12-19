function  resize_ima,ima,xbin=xbin,ybin=ybin,overscan=overscan
if not(keyword_set(overscan)) then overscan=0
a=size(ima) &
Nx=a(1)-overscan  & Ny=a(2)-overscan
out=fltarr(overscan+Nx/xbin,overscan+Ny/ybin)
tmp= congrid(ima(overscan:a(1)-1,overscan:a(2)-1),Nx/xbin,Ny/ybin)
out(overscan:Nx/xbin-1+overscan,overscan:Ny/ybin-1+overscan)=tmp
return,out
end
;dir='h:\obs_data.pol\s160406\'
;files=FILE_SEARCH(DIR+'*.fts')
;N=N_elements(files)
;for k=0,N-1 do begin
;ima=readfits(files(k),h)
;ima=resize_ima(ima,xbin=2,ybin=2,overscan=20)
;writefits,files(k),FIX(ima),h
;endfor
dir='h:\obs_data.pol\s161123\'
files=dir+'s147108'+string(findgen(17)+26,format='(I2.2)')+'.fts'
files=dir+'s147105'+string(findgen(9)+1,format='(I2.2)')+'.fts'
;files=dir+'s147113'+string(findgen(7)+1,format='(I2.2)')+'.fts'
print,files

N=N_elements(files)
for k=0,N-1 do begin
ima=readfits(files(k),h)
ima=resize_ima(ima,xbin=1,ybin=2,overscan=20)
sxaddpar,h,'BINNING','2x2'
writefits,files(k),FIX(ima),h
endfor
print,files
end