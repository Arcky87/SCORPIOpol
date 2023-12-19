;prepare date
pro prepare_data
FILELOG=DIALOG_PICKFILE(/read,path='d:\red.data\LOGS\',FILTER='*.txt')
create_initial_data_WOLL2,FILELOG

W_DIR=sxpar(read_table(FILELOG),'w_dir')
LOADFILE,dir=w_dir,ysz=500,xsz=1130
ViewPol_2
end
prepare_data
end