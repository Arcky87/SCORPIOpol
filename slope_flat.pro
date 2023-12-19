;проверка типа анализатора по плоскому полю
function slope_flat,wdir
flat=readfits(wdir+'flat_lin.fts',/silent)
a=size(flat)  & Nx=a(1)  & Ny=a(2) & Npol=a(3)  & wx=100
flat=total(flat(Nx/2-wx:Nx/2+wx,*,*),1)
slope=fltarr(Npol)
for k=0,Npol-1 do begin
f=goodpoly(findgen(Ny),flat(*,k),1,3)
slope(k)=f(1)
endfor
slope=FIX((slope(0)+slope(1)-slope(2)-slope(3))/total(ABS(slope)))
return,slope
end
log_dir='h:\red_data.pol\Sy1\'
LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path=log_dir+'LOGS\')

wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
wdir=log_dir+wdir(N_elements(wdir)-2)+'\'
print,wdir,slope_flat(wdir)
end