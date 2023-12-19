;input redshift in header
dir='h:\red_data.pol\'
;LOGFILE=pickfile(/read,path=dir+'LOGS\',filter='*.txt')
logfile=DIR+'LOGS\'+'E1821+643_140324.txt'
   	name_out=str_sep(FILE_BASENAME(LOGFILE),'.')
   		name_out=name_out(0)
			wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
				wdir=dir+wdir(N_elements(wdir)-2)+'\'
print,name_out
print,logfile
h=read_table(logfile)
print,sxPAR(H,'Rdir')


END
input=''
read,input,prompt='input value redshift  '
z=float(input)
;������ �������� �������� � LOGFILE
h=read_table(LOGFILE) & N=N_elements(h)
h=[h(0:3),'z       ='+string(z,format='(F34.6)')+' / REDSHIFT' ,h(N-6:N-1)]
openw,1,LOGFILE
for j=0,N_elements(h)-1 do printf,1,h(j)
close,1

file=['obj_i.fts','obj.fts','obj_lin.fts','obj-sky.fts','spectra.fit','stoks.fit']
for k=0,5 do begin
if FILE_TEST(wdir+file(k)) eq 1 then begin
h=headfits(wdir+file(k))
sxaddpar,h,'Z',z,after='NAME1',' /REDSHIFT'
modfits,wdir+file(k),0,h
;print,h
endif ELSE print,'MISSING FILE '+file(k)+'!!'
endfor
end