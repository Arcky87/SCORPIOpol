;view file
pro viewfile,dir
file=DIALOG_PICKFILE(PATH=dir,/read,filter='*.fts')

view_desktop,file
end