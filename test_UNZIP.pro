;test UNIZP+FTS


pro test_UNZIP
work_dir='h:\WOLLASTON-2.lib\'
zip=DIALOG_PICKFILE(/read,path='h:\obs_data.pol\',FILTER='*.zip')
;print,zip
zip_name=FILE_BASENAME(zip)
zip_dir=FILEPATH(zip)

SPAWN,'7z.exe l '+zip+' > '+work_dir+'zip',/LOG_OUTPUT
zip_tab=STRCOMPRESS(read_table(work_dir+'zip'))
N=N_elements(zip_tab)
for k=0,N-1 do begin
tmp=str_sep(zip_tab(k),' ')
M=N_elements(tmp)
zip_tab(k)=tmp(M-1)
endfor

zip_tab=zip_tab(12:N-3)  & N=N_elements(zip_tab)
for k=0,N-1 do begin

openw,1,work_dir+'list_zip\'+zip_tab(k)
printf,1,zip_tab(k)
close,1
endfor
fts=DIALOG_PICKFILE(/read,path=work_dir +'list_zip',FILTER='*.fts')
fts=FILE_BASENAME(fts)
print,'7z.exe x '+zip+' '+ fts
SPAWN,'7z.exe x '+zip+' '+ fts +' > info.txt'
FILE_DELETE,FILE_SEARCH(work_dir+'list_zip\'+'*.fts')

;file_copy,FTS,work_dir+'fts_zip\'+fts
;file_delete,FTS
view_desktop,work_dir+fts
file_delete,FTS
end
test_UNZIP
end