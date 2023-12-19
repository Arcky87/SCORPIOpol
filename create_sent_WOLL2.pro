;create_sent
function create_continuum,spectra,REP=rep
if not(keyword_set(rep)) then rep=4
N=N_elements(spectra)
for k=0,rep-1 do begin
cont=LOWESS(findgen(N),spectra,N/4,3)
robomean,spectra-cont,3,0.5,avg_cont,rms_cont
R=where(spectra-cont lt -rms_cont, ind)
if ind gt 0 then spectra(R)=cont(R)
endfor
return,spectra
end

dir='h:\red_data.pol\AGN\'
dir='h:\red_data.pol\Sy1\'
;dir='h:\red_data.pol\TINATIN2\'
path=dir+'LOGS\'
LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path=path)
   	name_out=str_sep(FILE_BASENAME(LOGFILE),'.')
   		name_out=name_out(0)
			wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
				wdir=dir+wdir(N_elements(wdir)-2)+'\'

spectra=readfits(wdir+'spectra.fit',h)
Nx=sxpar(h,'NAXIS1')
wave_0=sxpar(h,'CRVAL1')  & d_wave=sxpar(h,'CDELT1')
wave=findgen(Nx)*d_wave+wave_0
;определение имени таблицы
star_name=sxpar(h,'NAME2')
print,star_name
input=''

read,input,prompt='number of table?(BD+33d2642=0,BD+28d4211=1,g191b2b=2  '
table='d:\standards\data\'+['fbd33d2642','fbd28d4211','fg191b2b']+'.dat'
table=table(fix(input))

window,3,xsize=600,ysize=800,title='wdir  '+wdir+'  star  '+star_name
!P.multi=[0,1,3]
x=findgen(Nx)
star_obs=total(spectra(*,*,*,1),3)
star_obs=total(star_obs,2)
star_obs=median(star_obs,3)
plot,wave,star_obs,xst=1,charsize=2,title='star_obs',yrange=[0,1.1]*max(star_obs),yst=1
star_obs=create_continuum(star_obs)
oplot,wave,star_obs,thick=2
star_tab=read_st(table,wave,/print)
plot,wave,star_tab,xst=1,charsize=2,title='star_tab'
star_tab=create_continuum(star_tab)
oplot,wave,star_tab,thick=2
sent=star_obs/star_tab & sent=sent/sent(Nx/2)
plot,wave,sent,xst=1,charsize=2;,title='relative sensitivity',xtitle='Wavelength, A'
sent_smooth=LOWESS(findgen(Nx),sent,Nx/4,2,2)
oplot,wave,sent_smooth,thick=2

writefits,wdir+'sent.fts',sent_smooth
sent=readfits(wdir+'sent.fts',h)
sxaddpar,h,'CRVAL1',wave_0
sxaddpar,h,'CDELT1',d_wave
sxaddpar,h,'STAR',star_name
writefits,wdir+'sent.fts',sent,h

end