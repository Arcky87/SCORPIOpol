pro  LOADFILE,REP=rep,DIR=dir,YSZ=ysz,XSZ=xsz
COMMON IMAblock,file,h,frame,ima,WDIR,Nx,Ny,Nz,Nt,sc,bin,plane,target,d_lambda,lambda_0,over,x_size,y_size,range,level,avg_ima,rms_ima
COMMON WINblock,win_num,swin,planeW,draw,tra,ind_tra,planeT
if keyword_set(dir) then WDIR=dir
if keyword_set(xsz) then x_size=xsz
if keyword_set(ysz) then y_size=ysz
file=DIALOG_pickfile(PATH=wdir,/read,filter='*.fts')
frame=readfits(file,h)

;if os_family() eq 'unix' then sep='/' else sep='\'
;file=str_sep(file,sep) & file=file(N_elements(file)-1) & file=str_sep(file,'.') & file=file(0)
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2')
Nz=sxpar(h,'NAXIS3')
Nt=sxpar(h,'NAXIS4')
print,Nz,Nt
d_lambda=sxpar(h,'CDELT1')
lambda_0=sxpar(h,'CRVAL1')
sc=float(y_size)/float(Ny)
Naxis=sxpar(h,'NAXIS')
if Naxis eq 3 then begin
Nz=sxpar(h,'NAXIS3')
bin=1

if bin eq 0 then begin
bin=1
goto, cont1
endif
;bin=fix(st;r_sep(bin,'x')) & bin=float(bin(0))/bin(1)
cont1:
XY=float(Nx)/float(Ny)
sc=float(y_size)/float(Ny)
;x_size=y_size*XY*float(bin)
if keyword_set(rep) then begin
WIDGET_CONTROL, draw, GET_VALUE=win_num,DRAW_XSIZE=x_size
WIDGET_CONTROL, draw, GET_VALUE=win_num,DRAW_YSIZE=y_size
;WIDGET_CONTROL, draw, GET_VALUE=win_num,SET_DRAW_VIEW=[(x_size-y_size)/2,0]
print,x_size,y_size
endif
ima=congrid(frame,x_size,y_size,Nz)
endif
if Naxis eq 2 then begin
Nz=0 & ima=congrid(frame,x_size,y_size)
endif
if Naxis eq 4 then begin
tmp=fltarr(Nx,Ny,Nz)
ima=congrid(frame(*,*,*,target),x_size,y_size,Nz)
endif
plane=0
rms_ima=stdev(ima(*,*,plane),avg_ima)
;чтение файла траектории

if Nz ne 0 then begin
file_tra=STRMID(file,0,strlen(file)-strlen(FILE_BASENAME(FILE)))+'tra.fit'
if FILE_TEST(file_tra) eq 0 then tra=0 ELSE tra=readfits(file_tra)
endif
END