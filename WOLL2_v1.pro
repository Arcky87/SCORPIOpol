;SCORPIO_POL  for WOLL2 version    jul 2015
function MOD_value,value
value=value+1
value=value-2*(value/2)
return,value
end
;**********************************************************************
PRO WOLL2_RUN
  COMMON SET_widget,obj_pos_W,slitpos_W,yrange_W,overscan_W,beg_wave_W,end_wave_W,d_wave_W,BUTTONS
  COMMON SET_step,step,strob,x_slit,NScut,atm_lines,atm_window ,Ypos,Ywid,star,ampl,shift_abs
  COMMON SET_key,key1,key2,key3,key_1,key_2,key_3
  COMMON set_COL,background,color
  COMMON LOG,dLOG,text_log,log_dir,LOGFILE,idl_dir,gview,wdir
  COMMON LIN,beg_lambda,end_lambda,d_lambda,N_disp
  COMMON IMATYPE,ima1,ima3,ima5,ind_binX,ind_binY,ind_Neta
  COMMON PAR,par2,par5,par6_1,par6_2,par11_1,par11_2,par11_3
  COMMON SET_value,val
  COMMON SET_event,ev
  COMMON obj_pos,objpos,overscan,overscan_val,bias,bias_val,width,yrange_val
  COMMON GEOM,x_beg,w_x,N_tra,Deg_Tra,D_y,wg_x,wg_y,TRESH,H_slit,Y_shi
  COMMON EXTR,pos,w,comp,Ncomp,Rx,Ry,NSdeg
   COMMON STOKS,shift_PA,out_end,out_beg,W_integr,P_ISM,PA_ISM
;step#1 create traectory

 if step(1) eq 1 then begin

     traectory_WOLL2,wdir,WX=w_x,NDEG=Deg_Tra,Xbeg=x_beg,plot=key1(1)
 endif

;step#2 create geometry model

  if step(2) eq 1 then begin

  ;построение геометрической модели искажений

cube=readfits(wdir+'neon_i.fts',h)
bin=sxpar(h,'BINNING') & bin=str_sep(bin,'x') & bin=FLOAT(bin(1))
tra=readfits(wdir+'tra.fit')
;Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2')
geometry_coeff=fltarr(1000,4,4)
angle=' for   polarization angle = '+['0 deg','90 deg','45 deg','135 deg']
print,'D_y=',D_y
for k=0,3 do begin
geometry_WOLL2,cube(*,*,k),tra(*,*,k),X0,Y0,X1,Y1,plot=key1(2),scale=4,tresh=1,TITLE=angle(k),Dy=D_y ;default sc=4, thresh=1, no eps keyword IY
wait,1
Nc=N_elements(X0)
coords=[X0,Y0,X1,Y1]

C=reform(coords,Nc,4)
geometry_coeff(0:Nc-1,*,k)=C

endfor
writefits,wdir+'geometry.fit',geometry_coeff
print,'Model geometry is created!'
  endif

;step#3 correction geometry
 if step(3 ) eq 1 then begin
    ;коррекция геометрических искажений
angle=['0 deg','90 deg','45 deg','135 deg']
name_target=['  obj','  star_0','  star']
;****************************
scale=0.357 ; масштаб вдоль щели
;****************************
;Hslit =80 ;высота щели в рх *
;****************************
xpk0=[-25.75,0,25.75] ; координаты точечного теста в arcsec
cube_eta=readfits(wdir+'eta_i.fts',h_eta)
cube_neon=readfits(wdir+'neon_i.fts',h_neon) & print,'cube neon:  ',size(cube_neon)
cube_flat=readfits(wdir+'flat_i.fts',h_flat)
 cube_obj=readfits(wdir+'obj_i.fts',h_obj)
 Nx=sxpar(h_obj,'NAXIS1')
 Ny=sxpar(h_obj,'NAXIS2')
 Npol=sxpar(h_obj,'NAXIS3')
 Nexp=sxpar(h_obj,'NAXIS4')
 Ncub=sxpar(h_obj,'NAXIS5')
 eta_corr=fltarr(Nx,H_slit,Npol)
neon_corr=fltarr(Nx,H_slit,Npol)
flat_corr=fltarr(Nx,H_slit,Npol)
 obj_corr=fltarr(Nx,H_slit,Npol,Nexp,Ncub)
geometry_coeff=readfits(wdir+'geometry.fit')

FOR k=0,Npol-1 DO BEGIN
R=where(geometry_coeff(*,0,k) gt 0,ind)
C=fltarr(ind,4) & C(*,*)=geometry_coeff(0:ind-1,*,k)
;C=readfits(dir+'geometry'+string(k+1,format='(I1)')+'.fit')
;определение центра щели
eta=cube_eta(*,*,k)
eta_new=WARP_TRI(c(*,0),c(*,1),c(*,2),c(*,3),eta)
wx=10  & y=findgen(Ny) & wy=3
Vy=total(eta_new(Nx/2-Nx/8:Nx/2+Nx/8,*),1)
d_y=100
;Vy(0:d_y)=0  & Vy(Ny-d_y:Ny-1)=0
fi_peak,y,Vy,max(Vy)/5,ipix,xpk,ypk,bkpk,ipk
PRINT,'IPK',IPK
window,0

plot,y,Vy,xst=1
for j=0,ipk-1  do begin
f=goodpoly(y(xpk(j)-wy:xpk(j)+wy),Vy(xpk(j)-wy:xpk(j)+wy),2,3,fit)
xpk(j)=-f(1)/f(2)/2
endfor
xpk0=xpk0/scale
xpk=xpk
Yc=xpk(1)
print,'Yc=',Yc
PRINT,'h_SLIT', H_SLIT
neon_corr(*,*,k)=correction_image(cube_neon(*,*,k),C,Yc=xpk(1),Hslit=H_slit,plot=key1(3),TITLE='neon '+angle(k)) & WAIT,key1(3)
flat_corr(*,*,k)=correction_image(cube_flat(*,*,k),C,Yc=xpk(1),Hslit=h_slit,plot=key1(3),TITLE='flat '+angle(k)) & WAIT,key1(3)
 eta_corr(*,*,k)=correction_image(cube_eta(*,*,k),C,Yc=xpk(1),Hslit=h_slit,plot=key1(3),TITLE='eta '+angle(k))  & WAIT,key1(3)
 ;exposure
 for j=0,Ncub-1 do begin
 for i=0,Nexp-1 do begin
 if total(cube_obj(*,*,k,i,j)) ne 0 then begin
; if k eq 0 then cube_obj(*,*,k,i,j)=remove_bad_row(cube_obj(*,*,k,i,j),[66,67])
 obj_corr(*,*,k,i,j)=correction_image(cube_obj(*,*,k,i,j),C,Yc=xpk(1),Hslit=H_slit,plot=key1(3),TITLE='obj '+angle(k)+' exp='+string(i+1,format='(I2)')+'  cube='+string(j+1,format='(I2)') )
   WAIT,key1(3)
print,angle(k),name_target(j),'   exposure',i
    endif
 endfor
 endfor

ENDFOR
writefits,wdir+'eta.fts',eta_corr,h_eta
writefits,wdir+'neon.fts',neon_corr,h_neon
writefits,wdir+'flat.fts',flat_corr,h_flat
writefits,wdir+'obj.fts',obj_corr,h_obj
print,'Correction is made successfully!'

endif

;step#4 dispersion curve creation
    if step(4 ) eq 1 then begin


  if FILE_TEST(wdir+'etalon.txt') eq 1 then  goto,disp_1 ELSE $
  Res= DIALOG_MESSAGE( 'No file etalon.txt, you want to create it ?', /QUESTION)
  if res eq 'No' then goto, disp_2 ELSE  create_etalon_WOLL2,wdir     ;comment: before 'wdir' was 'LOGFILE'
  disp_1:

  neon=readfits(wdir+'neon.fts',h)
  grating=sxpar(h,'DISPERSE')
  Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'NAXIS2')
  Npol=sxpar(h,'NAXIS3') & Ndeg=3
  ;построение дисперсионной кривой
  ;чтение таблицы линий
  Ntab=numlines(wdir+'etalon.txt')
  tab=fltarr(5,Ntab)
  openr,1,wdir+'etalon.txt'
  readf,1,tab
  close,1
  d_wave=tab(0,0) & min_wave=tab(1,0) & max_wave=tab(2,0)
  tab=tab(*,1:Ntab-1) & Ntab=Ntab-1
  ident_table=fltarr(2,Ntab)
  ident_table(0,*)=tab(0,*)
  Disp=dblarr(Ny,N_disp+1,Npol)
  FOR k=0,Npol-1 DO BEGIN
  ident_table(1,*)=tab(k+1,*)
  if key1(4) eq 0 then DISP(*,*,k)=dispersion_WOLL2(neon(*,*,k),ident_table,N_DEG=N_disp)  ELSE $
  DISP(*,*,k)=dispersion_WOLL2(neon(*,*,k),ident_table,N_DEG=N_disp,plot=k+1)
  ENDFOR

  ;вычисление пределов линеаризации
  wave=fltarr(Nx,Npol)
  for k=0,Npol-1 do begin
  D=total(disp,1)/Ny
  tmp=0 & for j=0,Ndeg do tmp=tmp+D(j,k)*findgen(Nx)^j
  wave(*,k)=tmp
  endfor
  beg_lambda=FIX(max(wave(0,*))/10+1)*10
  end_lambda=FIX(min(wave(Nx-1,*))/10)*10

  robomean,wave-shift(wave,1),3,0.5,avg_dwave
  ;d_lambda=FIX(avg_dwave*2)*0.5
  d_lambda=avg_dwave

  N_lin=FIX(end_lambda-beg_lambda)/d_lambda+1
  mkhdr,h,disp
    sxaddpar,h,'NAXIS1',Ny,' LENGTH OF SLIT, px'
  sxaddpar,h,'NAXIS2',Ndeg+1,' DEGREE OF POLYNOM + 1'
  sxaddpar,h,'NAXIS3',Npol,'POLARIZATION PLANE (0,90,45 and 135 degree) '
  sxaddpar,h,'lambda_0',beg_lambda,' INITIAL WAVELENGTH, A'
   print, 'lambda_0',beg_lambda,' INITIAL WAVELENGTH, A'
  sxaddpar,h,'d_lambda',d_lambda,'  Dispersion, A/px'
   print, 'd_lambda',d_lambda,'  Dispersion, A/px'
  sxaddpar,h,'N_lin',N_lin,' Length linearised spectra'
  sxaddpar,h,'grating',grating,' Name of grating'
  writefits,wdir+'dispersion.fit',disp,h
  ;
    WIDGET_CONTROL,beg_wave_W,SET_VALUE=string(beg_lambda,format='(I4)')
  WIDGET_CONTROL,end_wave_W,SET_VALUE=string(end_lambda,format='(I4)')
  WIDGET_CONTROL,d_wave_W,SET_VALUE=string(d_lambda,format='(F3.1)')
  disp_2:
    endif
;step#5 linerisation
  if step(5) eq 1 then begin
     DISP=readfits(wdir+'dispersion.fit',hd)
  print,'DISP:  ', size(DISP)
  stop
     Ny=sxpar(hd,'NAXIS1') & Npol=sxpar(hd,'NAXIS3')
  Nlin=FIX((end_lambda-beg_lambda)/d_lambda)+1
  PARAM_LIN=[beg_lambda,d_lambda,Nlin]
  print,'start linerization ',PARAM_lin, format='(A20,I6,F5.2,I6)'
  ;goto,fin_lin
   ;линеаризация  калибровочных спектров
  cube_lin=fltarr(Nlin,Ny,Npol)
  type=['eta','flat','neon']

  for k=0,2 do begin
  print,'linerization  '+type(k)

     cube_ini=readfits(wdir+type(k)+'.fts',h,/silent)
  for j=0,Npol-1 do begin

  ima=linerisation_WOLL2(cube_ini(*,*,j),DISP(*,*,j),PARAM=[beg_lambda,d_lambda,Nlin])

  robomean,congrid(ima,Nlin,Ny),3,0.5,avg_ima,rms_ima

  cube_lin(*,*,j)=ima
  endfor
   sxaddpar,h,'Naxis1',Nlin
   sxaddpar,h,'CRVAL1',beg_lambda,' INITIAL WAVELENGTH, A'
   sxaddpar,h,'CDELT1',d_lambda,AFTER='CRVAL1',' DISPERSION, A/px'
   writefits,wdir+type(k)+'_lin.fts',cube_lin,h
  ENDFOR
  ;линеаризация  спектра объекта
  print,'linerization spectra object'
  cube_obj=readfits(wdir+'obj.fts',h,/silent)
  Npol=4
  Nexp=sxpar(h,'NAXIS4')
  print,Nexp
  Ncub=3
  cube_lin=fltarr(Nlin,Ny,Npol,Nexp,Ncub)

  for i=0,Ncub-1 do begin

  for j=0,Nexp-1 do begin
  if total(cube_obj(*,*,*,j,i)) ne 0 then begin
  ;Window,2,xsize=Nlin/2,ysize=(Ny+1)*4,title=dir+'   cube '+string(i,format='(I2)')+'  exp '+string(j,format='(I3)')
  for k=0,Npol-1 do begin
  ima=linerisation_WOLL2(cube_obj(*,*,k,j,i),DISP(*,*,k),PARAM=[beg_lambda,d_lambda,Nlin])

  cube_lin(*,*,k,j,i)=ima
  endfor
  endif
  endfor
    endfor
  sxaddpar,h,'Naxis1',Nlin
  sxaddpar,h,'CRVAL1',beg_lambda,' INITIAL WAVELENGTH, A'
  sxaddpar,h,'CDELT1',d_lambda,AFTER='CRVAL1',' DISPERSION, A/px'
  writefits,wdir+'obj_lin.fts',cube_lin,h
  fin_lin:
  endif

;step#6 N.S. substraction
    if step(6) eq 1 then begin

    create_flat_WOLL2,wdir,plot=key1(6)
    print,'sky_substraction'
    print,'NSdeg',NSdeg,'NScut',NScut
    subsky_WOLL2,wdir,plot=key1(6),/flat,/Yshift;,NSDEG=NSdeg
     endif
;step#7 spectra extraction
    if step(7) eq 1 then begin

    ;проверка наличия файла
    lines_atm=FLOAT(str_sep(atm_lines,','))
    print,'atm_lines',lines_atm
    if  file_test(wdir+'obj-sky.fts') eq 1 then begin

    if key2(7) eq 1 then atm_absorbtion_WOLL2,LOGFILE,LINE=lines_atm,region=atm_window,star=star;,plot=key1(7)
  ;step#7 spectra extraction
    if step(7) eq 1 then begin
    print,'strob',strob
    type=[4,0,1]
    cube=readfits(wdir+'obj-sky.fts',h)
    print,'Rx=',Rx,'  Ry=',Ry
  print, 'Ywid=',Ywid,' Ypos=',Ypos
  ;Ywid=40
     extraction_WOLL2,wdir,plot=key1(7),WY=Ywid,Yc=Ypos,ATM_ABS=key2(7),AMPL=ampl,SHIFT_X=shift_abs ;,YC=yc,WY=wy
  stoks_WOLL2,wdir
    endif

  endif
      endif



;step#8 create Stoks vectors
  if step(8) eq 1 then begin
        ;create_stocks_WOLL2,wdir
     ;avg_stoks_analyz
  endif
;step#9 determination redshhift
    if step(9) eq 1 then begin
    z_measure,LOGFILE,idl_dir
    endif
  goto,fin
    err_log:
    print,'LOGFILE is empty!!!
    fin:
END
;************************************************************************************
PRO WOLL2_Event, Event
  COMMON SET_widget,obj_pos_W,slitpos_W,yrange_W,overscan_W,beg_wave_W,end_wave_W,d_wave_W,BUTTONS
  COMMON SET_step,step,strob,x_slit,NScut,atm_lines,atm_window ,Ypos,Ywid,star,ampl,shift_abs
 COMMON SET_key,key1,key2,key3,key_1,key_2,key_3
  COMMON set_COL,background,color
  COMMON LOG,dLOG,text_log,log_dir,LOGFILE,idl_dir,gview,wdir
  COMMON IMATYPE,ima1,ima3,ima5,ind_binX,ind_binY,ind_Neta
  COMMON PAR,par2,par5,par6_1,par6_2,par11_1,par11_2,par11_3
 COMMON LIN,beg_lambda,end_lambda,d_lambda,N_disp
  COMMON SET_value,val
  COMMON SET_event,ev
  COMMON obj_pos,objpos,overscan,overscan_val,bias,bias_val,width,yrange_val
  COMMON GEOM,x_beg,w_x,N_tra,Deg_Tra,D_y,wg_x,wg_y,TRESH,H_slit,Y_shi
  COMMON EXTR,pos,w,comp,Ncomp,Rx,Ry,NSdeg
   COMMON STOKS,shift_PA,out_end,out_beg,W_integr,P_ISM,PA_ISM
  ;!P.background=background  & !P.color=color
 ;set INIT value keywords
 key_1=['no','no','plot','plot','plot','plot','no','erg','PS',  'plot','plot']
 key_2=['no','no',  'no',  'no',  'no',  'no','no','mJy','no','manual', 'mJy']
 key_3=['no','no',  'no',  'no',  'no',  'no','no','no', 'no',    'no',  'no']
 R=where(key1 eq 0) & key_1(R)='no' & R=where(key2 eq 0) & key_2(R)='no'

  WIDGET_CONTROL,Event.Id,GET_UVALUE=Ev

  CASE Ev OF
  'extr'      :  BEGIN
       WIDGET_CONTROL,event.id,GET_VALUE=strob
       ;if strob eq 1 then for k=15,20  do WIDGET_CONTROL,BUTTONS(k),SENSITIVE=1 ELSE $
       ;	for k=15,20  do WIDGET_CONTROL,BUTTONS(k),SENSITIVE=0

       END
  ;reading & edition LOG-file
  'set_PATH':   BEGIN
      log_dir=DIALOG_PICKFILE(PATH=log_dir,/DIRECTORY)
      WIDGET_CONTROL, dLOG, SET_VALUE=log_dir
      end

  'LOG_input':  BEGIN
          ;filename=PICKFILE(/read,PATH=log_dir+'LOGS\',filter='*.txt')
          WIDGET_CONTROL,BUTTONS(17),sensitive=0
          CD,log_dir+'LOGS/',CURRENT=scorp_dir
       filename=pickfile(/read,filter='*.txt')
       CD,scorp_dir
    LOGFILE=filename
       print,logfile
      name_out=FILE_BASENAME(filename)
      WIDGET_CONTROL, text_LOG, SET_VALUE=name_out
    print,sxpar(read_table(LOGFILE),'W_DIR','/')
    wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'/')
    wdir=log_dir+wdir(N_elements(wdir)-2)+'/'
    print,wdir

      ;goto,cont1
       for j=0,8 do WIDGET_CONTROL,BUTTONS(J),SENSITIVE=1
       ;
      ;wdir=DEF_WDIR(LOGFILE)
    ;print,FILE_TEST(wdir+'log.fit')
      if   FILE_TEST(wdir+ 'eta_i.fts') eq 1 and FILE_TEST(wdir+'neon_i.fts')$
      and  FILE_TEST(wdir+'flat_i.fts') eq 1 and FILE_TEST(wdir+'obj_i.fts')$
       then WIDGET_CONTROL,BUTTONS(3),SENSITIVE=0 $
        ELSE WIDGET_CONTROL,BUTTONS(3),SENSITIVE=1
    if   FILE_TEST(wdir+ 'obj-sky.fts') eq 1 then begin
    Ypos=sxpar(headfits(wdir+ 'obj-sky.fts'),'NAXIS2')/2
    Ywid=20
     WIDGET_CONTROL,BUTTONS(17),sensitive=1
    WIDGET_CONTROL,BUTTONS(18),SET_VALUE=string(Ypos,format='(I2)')
    WIDGET_CONTROL,BUTTONS(19),SET_VALUE=string(Ywid,format='(I2)')
    endif
    ;if FILE_TEST(wdir+'log.fit') eq 1 then WIDGET_CONTROL,BUTTONS(1),SENSITIVE=0 $
              ; ELSE WIDGET_CONTROL,BUTTONS(1),SENSITIVE=1

        ;x_slit=string(def_slit(logfile),format='(I4)')
        ;WIDGET_CONTROL,slitpos_W,SET_VALUE=x_slit
        ;obj_name=def_name(LOGFILE,'obj',Nobj)
        ;tmp_obj=READFILE(def_rdir(LOGFILE),obj_name(0),def_ext(LOGFILE),h)
      ;		Ny=sxpar(h,'NAXIS2') & yrange_val='0,'+string(Ny-1,format='(I4)')
        ;WIDGET_CONTRO,yrange_W,SET_VALUE=yrange_val
      ;openw,1,'wdir.txt'
      ;printf,1,def_wdir(LOGFILE)
      ;close,1
     if FILE_TEST(def_wdir(LOGFILE)+'dispersion.fit') eq  1 then begin
     hd=headfits(def_wdir(LOGFILE)+'dispersion.fit')
     beg_lambda=sxpar(hd,'lambda_0') & d_lambda=sxpar(hd,'d_lambda')
     end_lambda=(sxpar(hd,'N_lin')-1)*d_lambda+beg_lambda
     WIDGET_CONTROL,beg_wave_W,SET_VALUE=string(beg_lambda,format='(I5)')
     WIDGET_CONTROL,end_wave_W,SET_VALUE=string(end_lambda,format='(I5)')
     WIDGET_CONTROL,d_wave_W,SET_VALUE=string(d_lambda,format='(F3.1)')
     endif
    cont1:
       end
  'LOG_edit':begin & editor,LOGFILE & end


  'Yshift': begin
         WIDGET_CONTROL,event.id,GET_VALUE=T
         Y_shi=FLOAT(T(0))
         end
  'del_FTS': BEGIN
     tmp=FILE_SEARCH(wdir+'*.fts')
     if tmp(0) ne '' then  FILE_DELETE,FILE_SEARCH(wdir+'*.fts')
     WIDGET_CONTROL,BUTTONS(2),SENSITIVE=0
     WIDGET_CONTROL,BUTTONS(3),SENSITIVE=1
     end
'datacube':BEGIN

   create_initial_data_WOLL2,LOGFILE,Y_shift=Y_shi
   WIDGET_CONTROL,BUTTONS(3),SENSITIVE=0
   WIDGET_CONTROL,BUTTONS(2),SENSITIVE=1
   END
'binX':  begin
   ind_binX=event.index
   print,ind_binX
   end
'binY': ind_binY=event.index
'Neta': ind_Neta=event.index
  'view':   begin
     LOADFILE,dir=wdir,ysz=500,xsz=1130
   ViewPol_2
     end
  ;set step reduction, keywords & parameters

; create traectory	#1
  val(1):   begin & step(1)=mod_value(step(1))
     if step(1) eq 0 then Wdelete,22
     if step(1) eq 0 then Wdelete,2
     end
        'k0':   begin & key2(1)=mod_value(key2(1)) & END    ;plot
 ;'k1':   begin &	key1(1)=mod_value(key1(1)) & END    ;correction geometry
'Xbeg':  begin & WIDGET_CONTROL,event.id,GET_VALUE=T & x_beg=FLOAT(T(0)) & end
'wx':   begin & WIDGET_CONTROL,event.id,GET_VALUE=T & w_x=FLOAT(T(0)) & END
'Ntra':  begin & WIDGET_CONTROL,event.id,GET_VALUE=T & N_tra=FLOAT(T(0)) & END
'DegTra':  begin & WIDGET_CONTROL,event.id,GET_VALUE=T & Deg_Tra=FLOAT(T(0)) & END

;create geometry model  #2
'Dy':   begin & WIDGET_CONTROL,event.id,GET_VALUE=T & D_y=FLOAT(T(0)) & end
'wgx':   begin & WIDGET_CONTROL,event.id,GET_VALUE=T & wg_x=FLOAT(T(0)) & END
'wgy':    begin & WIDGET_CONTROL,event.id,GET_VALUE=T & wg_y=FLOAT(T(0)) & END
'trsh':  begin & WIDGET_CONTROL,event.id,GET_VALUE=T & TRESH=FLOAT(T(0)) & END

  val(2):  begin
     step(2)=mod_value(step(2))
     if step(2) eq 0 then begin
     Wdelete,2
     Wdelete,3
     endif
     END
   'k2':   begin & key1(2)=mod_value(key1(2)) & END
   'k_2':   begin & key2(2)=mod_value(key2(2)) & END
    'p2':   BEGIN WIDGET_CONTROL,event.id,GET_VALUE=par2 & end

;correction geometry	#3

  val(3):  begin & step(3)=mod_value(step(3)) & END
    'k3':   begin & key1(3)=mod_value(key1(3)) & END
   'k_3':   begin & key2(3)=mod_value(key2(3)) & END
 'Hslit':  begin & WIDGET_CONTROL,event.id,GET_VALUE=T & H_slit=FLOAT(T(0)) & END


; create 2D-dispersion curve	#4
  val(4):  begin
     step(4)=mod_value(step(4))
     if step(4) eq 0 then begin
     for j=0,4 do Wdelete,j
     endif
     END
'DEG_DISP': begin & WIDGET_CONTROL,event.id,GET_VALUE=T & N_disp=FLOAT(T(0)) & END
    'k4':   begin & key1(4)=mod_value(key1(4)) & END
   'k_4':   begin & key2(4)=mod_value(key2(4)) & END

;   linerization spectra
  val(5):  begin & step(5)=mod_value(step(5))
     if step(5) eq 1 then begin
     if file_test(wdir+'dispersion.fit') eq 0 then begin
     R=Dialog_message(/error,'CANNOT FIND dispersion curve')
     RETURN
     endif
     hd=headfits(wdir+'dispersion.fit')
     beg_lambda=sxpar(hd,'lambda_0')
   end_lambda=(sxpar(hd,'N_lin')-1)*sxpar(hd,'d_lambda')+beg_lambda

  WIDGET_CONTROL,beg_wave_W,SET_VALUE=string(beg_lambda,format='(I4)')
  WIDGET_CONTROL,end_wave_W,SET_VALUE=string(end_lambda,format='(I5)')
  WIDGET_CONTROL,d_wave_W,SET_VALUE=string(d_lambda,format='(F3.1)')
     endif
     END ;#5
    'k5':   begin & key1(5)=mod_value(key1(5)) & END
   'k5m':   begin & key2(5)=mod_value(key2(5)) & END
'beg_wave': begin & WIDGET_CONTROL,event.id,GET_VALUE=t & beg_lambda=float(t(0)) & END
'end_wave': begin & WIDGET_CONTROL,event.id,GET_VALUE=t & end_lambda=float(t(0)) & END
'd_wave': begin & WIDGET_CONTROL,event.id,GET_VALUE=t & d_lambda=float(t(0)) & END




;   sky substraction
  val(6):  begin & step(6)=mod_value(step(6)) & END ;#6
    'k6':   begin & key1(6)=mod_value(key1(6)) & END

'NSdeg':begin & WIDGET_CONTROL,event.id,GET_VALUE=t & NSdeg=float(t(0)) & END
'NS_cut': begin & WIDGET_CONTROL,event.id,GET_VALUE=mm & NScut=mm & END


;extraction spectra


;#7



  val(7):  begin & step(7)=mod_value(step(7)) & END
    'k7':   begin & key1(7)=mod_value(key1(7)) & END
   'k_7': begin
      key2(7)=mod_value(key2(7))
      print,key2(7)
      if key2(7) eq 0 then WIDGET_CONTROL,BUTTONS(16),SENSITIVE=0 ELSE WIDGET_CONTROL,BUTTONS(16),SENSITIVE=1
      if key2(7) eq 0 then WIDGET_CONTROL,BUTTONS(20),SENSITIVE=0 ELSE WIDGET_CONTROL,BUTTONS(20),SENSITIVE=1
      END

'strob_width': begin & WIDGET_CONTROL,event.id,GET_VALUE=pp & width=pp & END
'Ypos': begin
  WIDGET_CONTROL,event.id,GET_VALUE=EE
  Ypos=FIX(EE(0)) & print,'Ypos=',Ypos
  END
'Ywid':  begin
  WIDGET_CONTROL,event.id,GET_VALUE=EE
  Ywid=FIX(EE(0))
  END
'atm_line':begin
  WIDGET_CONTROL,event.id,GET_VALUE=EE
  atm_lines=EE
  END
'atm_win': begin
  WIDGET_CONTROL,event.id,GET_VALUE=EE
  NSwindow=FIX(EE(0))
  END
'star': begin
  WIDGET_CONTROL,event.id,GET_VALUE=EE
  star=FIX(EE(0))
  END
'ampl': begin
  WIDGET_CONTROL,event.id,GET_VALUE=EE
  ampl=EE(0) & print,ampl
  END
'shift':begin
  WIDGET_CONTROL,event.id,GET_VALUE=EE
  shift_abs=EE(0) & print,shift_abs
  END
;create stoks


  val(8):  begin & step(8)=mod_value(step(8)) & END ;#8
    'k8':   begin & key1(8)=mod_value(key1(8)) & END
   'k_8':   begin & key2(8)=mod_value(key2(8)) & END
 'shift_PA': begin & WIDGET_CONTROL,event.id,GET_VALUE=aa & shift_PA=aa & END
 'out_end': begin & WIDGET_CONTROL,event.id,GET_VALUE=aa & out_end=aa & END
 'out_beg': begin & WIDGET_CONTROL,event.id,GET_VALUE=aa & out_beg=aa & END
 'integr': begin & WIDGET_CONTROL,event.id,GET_VALUE=aa & W_integr=aa & END
 'P_ISM': begin & WIDGET_CONTROL,event.id,GET_VALUE=aa & P_ISM=aa & END
 'PA_ISM': begin & WIDGET_CONTROL,event.id,GET_VALUE=aa & PA_ISM=aa & END

  val(9):  begin & step(9)=mod_value(step(9)) & END ;#9
    'k9':   begin & key1(9)=mod_value(key1(9)) & END

  'RUN' :  BEGIN & WOLL2_RUN & END

  'HELP':BEGIN & xdisplayfile,'mpfs_help.txt' & END
  'XLOA':BEGIN xloadct & END
  'EXIT': BEGIN

          WIDGET_CONTROL,/DESTROY,event.top
          END

  ENDCASE


END ; MAIN_EVENT

PRO WOLL2 ; MAIN WIDGETS-INTERFACE
  COMMON SET_widget,obj_pos_W,slitpos_W,yrange_W,overscan_W,beg_wave_W,end_wave_W,d_wave_W,BUTTONS
  COMMON SET_step,step,strob,x_slit,NScut,atm_lines,atm_window ,Ypos,Ywid,star,ampl,shift_abs
   COMMON SET_key,key1,key2,key3,key_1,key_2,key_3
  COMMON set_COL,background,color
  COMMON LOG,dLOG,text_log,log_dir,LOGFILE,idl_dir,gview,wdir
  COMMON IMATYPE,ima1,ima3,ima5,ind_binX,ind_binY,ind_Neta
  COMMON PAR,par2,par5,par6_1,par6_2,par11_1,par11_2,par11_3
  COMMON SET_VALUE,val
  COMMON SET_event,ev
   COMMON LIN,beg_lambda,end_lambda,d_lambda,N_disp
  COMMON obj_pos,objpos,overscan,overscan_val,bias,bias_val,width,yrange_val
  COMMON GEOM,x_beg,w_x,N_tra,Deg_Tra,D_y,wg_x,wg_y,TRESH,H_slit,Y_shi
  COMMON EXTR,pos,w,comp,Ncomp,Rx,Ry,NSdeg
   COMMON STOKS,shift_PA,out_end,out_beg,W_integr,P_ISM,PA_ISM
; set init parameters
;******************************************************
   log_dir='/hdd/Glagol/2019/WOLL-2/20191216/'
;*****************************************************
background=2^24-1 & color=0
Y_shi=0  & NSdeg=2 & NScut=3 & atm_lines='6890,7260,7620' & atm_window=300
strob=0 &  width=7
Ncomp=1 & comp=intarr(3) & pos=fltarr(3) & w=fltarr(3)
LOGFILE='empty.txt'
ind_binX=0
ind_binY=0
ind_Neta=1
Rx=1 & Ry=20
  ima1=0 & ima3=0 & ima5=0
  LOGFILE='empty.txt'
  val=['1','2','3','4','5','6','7','8','9','10','11']
  par2='3' & par5='3' & par6_1='3' & par6_2='3'& par7_1='4863.18'
  overscan_val=20 & objpos=0 & x_slit=0
  sw1=[0,1,0]
  step_title=['create traectory           ',$
           'create geometry model      ',$
           'correction geometry        ',$
           'create 2D-dispersion curve ',$
          'linearisation spectra',$
              'nigth sky substraction                                ',$
       'extraction spectra  ',$
       'create Stoks vectors    ',$
       'correction interstellar polarization                                ',$
       'determination redshift                 ',$
       'flux calibration                                   ']

    x1=21.5
 MAINBASE = WIDGET_BASE(col=1,MAP=1,xsize=x1*54,TITLE='Data reduction for Double Wollaston analyzer mode SCORPIO-2 Spectrograph')
base=WIDGET_BASE(mainbase,row=1,frame=1)

;SET directory LOGFILE,open LOGFILE, edition LOGFILE
base_LOG=WIDGET_BASE(mainbase,row=1,frame=1)
label_LOG=WIDGET_BUTTON(base_log,value='PATH',uvalue='set_PATH');,xsize=x1*4);,/align_center)

dLOG=WIDGET_TEXT(BASE_log,value=log_dir,XSIZE=25,/align_center)


BUT_LOG=WIDGET_BUTTON(BASE_LOG,VALUE='open LOG-file',uvalue='LOG_input')

text_LOG=WIDGET_TEXT(base_log, XSIZE=18,/align_center)

edit_LOG=WIDGET_BUTTON(BASE_LOG,VALUE='edit LOG-file',uvalue='LOG_edit',xsize=x1*3.7)

LABEL=widget_label(base_LOG,value='  Shift Y',/align_center)
Yshift=WIDGET_TEXT(BASE_LOG,VALUE='0',uvalue='Yshift',/edit,/all,/align_center,xsize=3)

BUT_DEL_FTS=WIDGET_BUTTON(BASE_LOG,VALUE='delete FTS',uvalue='del_FTS',xsize=x1*3.3)
BUT_CREATE_datacube=WIDGET_BUTTON(BASE_LOG,VALUE='create DATACUBE',uvalue='datacube',xsize=x1*5.5)
BUT_BINX=WIDGET_DROPLIST(BASE_LOG,VALUE=['1','2','4'],uvalue='binX',title='binX')
BUT_BINY=WIDGET_DROPLIST(BASE_LOG,VALUE=['1','2','4'],uvalue='binY',title='binY')
BUT_Neta=WIDGET_DROPLIST(BASE_LOG,VALUE=['3','5','13'],uvalue='Neta',title='Neta')
BUT_VIEW=WIDGET_BUTTON(BASE_LOG,VALUE='view frame',uvalue='view',xsize=x1*3.5)


;common title
BASE_0=WIDGET_BASE(mainbase,row=1,frame=1)

tmp=WIDGET_LABEL(BASE_0,value='STEP OF REDUCTION',xsize=x1*16,/align_center)
tmp=WIDGET_LABEL(BASE_0,value='KEYWORDS',xsize=x1*12,/align_center)
tmp=WIDGET_LABEL(BASE_0,value='PARAMETERS',xsize=x1*12,/align_center)

     ;working directory
ypd=4
for k=1,9 do begin
   BASE_1=WIDGET_BASE(mainbase,row=1,frame=1)
scale=25
if k eq 7 then scale=10
if k eq 7 then ypd=30 else ypd=4
if k eq 8 then scale=10
;if k eq 1 then scale=15
bgroup = cw_bgroup(base_1,step_title(k-1),/nonEXCLUSIVE,uvalue=val(k),xsize=x1*scale,set_value=step(k),/col,ypad=ypd)

if k eq 1 then begin
yrange_val='0,1024'
;bgroup = cw_bgroup(base_1, 'correction geometry',/nonEXCLUSIVE,uvalue='k1',set_value=key1(k))
bgroup = cw_bgroup(base_1,'plot',/nonEXCLUSIVE,uvalue='k0',set_value=key2(k))
x_beg=100 & w_x=40 & N_tra=12 & Deg_Tra=3
 LABEL=widget_label(base_1,value='  Start X ',/align_center)
 Xbeg=WIDGET_TEXT(BASE_1,VALUE=string(x_beg),uvalue='Xbeg',/edit,/all,/align_center,xsize=6)
 LABEL=widget_label(base_1,value='  Window at X ',/align_center)
 Wx=WIDGET_TEXT(BASE_1,VALUE=string(w_x),uvalue='wx',xsize=6,/edit,/all,/align_center)
 LABEL=widget_label(base_1,value='  Number of traectory ',/align_center)
 Ntra=WIDGET_TEXT(BASE_1,VALUE=string(N_tra),uvalue='Ntra',/edit,/all,/align_center,xsize=5)
 LABEL=widget_label(base_1,value='  Degree ',/align_center)
 DegTra=WIDGET_TEXT(BASE_1,VALUE='3',uvalue='DegTra',xsize=5,/edit,/all,/align_center)
  endif
if k eq 2 then begin
  bgroup = cw_bgroup(base_1, 'plot',/nonEXCLUSIVE,uvalue='k2',set_value=key1(k))
D_y=15 & wg_x=30 & wg_y=20 & TRESH=1
 LABEL=widget_label(base_1,value=' Expand at Y ',/align_center)
 dY=WIDGET_TEXT(BASE_1,VALUE='40',uvalue='Dy',/edit,/all,/align_center,xsize=6)
 LABEL=widget_label(base_1,value=' Window at X ',/align_center)
 Wgx=WIDGET_TEXT(BASE_1,VALUE='30',uvalue='wgx',xsize=6,/edit,/all,/align_center)
 LABEL=widget_label(base_1,value=' Strob width Y ',/align_center)
 Wgy=WIDGET_TEXT(BASE_1,VALUE='20',uvalue='wgy',/edit,/all,/align_center,xsize=6)
 LABEL=widget_label(base_1,value= ' Tresh, rms ',/align_center)
 TRSH=WIDGET_TEXT(BASE_1,VALUE='1',uvalue='trsh',xsize=5,/edit,/all,/align_center)





endif

if k eq 3 then begin
H_slit=80
 bgroup = cw_bgroup(base_1, 'plot',/nonEXCLUSIVE,uvalue='k3',set_value=key1(k))
  LABEL=widget_label(base_1,value= 'Length of slit ',/align_center)
 Hslit=WIDGET_TEXT(BASE_1,VALUE=string(H_slit),uvalue='Hslit',xsize=6,/edit,/all,/align_center)

endif

if k eq 4 then begin
N_disp=3
 bgroup = cw_bgroup(base_1, 'plot result approximation',    /nonEXCLUSIVE,uvalue='k4',set_value=key1(k))
 LABEL=widget_label(base_1,value= '   Degree of polynomial approximation  dispersion',/align_center)
 DEG_DISP=WIDGET_TEXT(BASE_1,VALUE=string(N_disp),uvalue='DEG_DISP',xsize=6,/edit,/all,/align_center)
;	bgroup = cw_bgroup(base_1, 'plot result approximation',    /nonEXCLUSIVE,uvalue='k_4',set_value=key2(k))

endif

if k eq 5 then begin

spectral_range_lab=widget_label(base_1,value='spectra range, AA',/align_center)
beg_wave_W=WIDGET_TEXT(BASE_1,VALUE=string(0,format='(I5)'),uvalue='beg_wave',xsize=5,/edit,/all,/align_center)
end_wave_W=WIDGET_TEXT(BASE_1,VALUE=string(0,format='(I5)'),uvalue='end_wave',xsize=6,/edit,/all,/align_center)
disperse_lab=widget_label(base_1,value='dispersion, A/px',/align_center)
d_wave_W=WIDGET_TEXT(BASE_1,VALUE=string(0,format='(I4)'),uvalue='d_wave',xsize=5,/edit,/all,/align_center)

endif
if k eq 6 then begin
 bgroup = cw_bgroup(base_1,'plot',/nonEXCLUSIVE,uvalue='k6',set_value=key1(k))


NS_cut_LAB=widget_label(base_1,value='N.S. cutoff level (rms)  ',/align_center)
    NS_cut=WIDGET_TEXT(BASE_1,VALUE=string(NScut,format='(I3)'),uvalue='NS_cut',xsize=4,/edit,/all,/align_center)
    NS_deg_lab=widget_label(base_1,value='  degree polynomial approximation N.S,  ',/align_center)
 NS_deg=WIDGET_TEXT(BASE_1,VALUE=string(NSdeg,format='(I4)'),uvalue='NSdeg',xsize=5,/edit,/all,/align_center)

endif
if k eq 7 then begin
 bgroup = CW_BGROUP(base_1,['strob','gauss','moffat'],/EXCLUSIVE,uvalue='extr',set_value=strob,/col)
 reg=widget_base(base_1,/col,/frame,ypad=1)
 lab=widget_label(reg,value='REGION OF EXTRACTION',ysize=20,yoffset=0 )
   Y_pos=CW_FIELD(reg,value='',uvalue='Ypos',title='  position (px)  ',xsize=4,/ALL_EVENTS,/string )
 Y_width=CW_FIELD(reg,value='',uvalue='Ywid',title='      width (px)  ',xsize=4,/ALL_EVENTS,/string )
   ; strob_width_LAB=widget_label(base_1,value='width',/align_center,yoffset=-20 )

   ; strob_width=WIDGET_TEXT(BASE_1,VALUE=string(width,format='(I3)'),uvalue='strob_width',xsize=4,/edit,/all,/align_center)
   tmp_LAB=widget_label(base_1,value=' ',xsize=x1*1.3,/align_center)
   base_2=WIDGET_BASE(base_1,/col,ysize=1 )
;	empty=widget_label(base_1,value='  ')
 bgroup = cw_bgroup(base_2, 'plot',   /nonEXCLUSIVE,uvalue='k7',set_value=key1(k))
 bgroup=cw_bgroup(base_2, 'correction atmosphere',   /nonEXCLUSIVE,uvalue='k_7',set_value=key2(k) )
 base_3=WIDGET_BASE(base_1,/col,ysize=1,ypad=1)
 lines_atm=CW_FIELD(base_3,value=atm_lines,uvalue='atm_line',title='absorbtion lines: ',xsize=15,/ALL_EVENTS,/STRING )
 window_atm=CW_FIELD(base_3,value=string(atm_window,format='(I3)'),uvalue='atm_win',title='            window: ',xsize=5,/ALL_EVENTS,/STRING )
 base_4=WIDGET_BASE(base_1,/col,ysize=1,ypad=1)
 star=1 & ampl=1 & shift_abs=0
 star_number=CW_FIELD(base_4,value=string(star,format='(I3)'),uvalue='star',title='   target: ',xsize=3,/ALL_EVENTS,/STRING )
       atm_ampl=CW_FIELD(base_4,value=string(ampl,format='(I3)'),uvalue='ampl',title='amplifier: ',xsize=3,/ALL_EVENTS,/FLOAT )
      atm_shift=CW_FIELD(base_4,value=string(shift_abs,format='(I3)'),uvalue='shift',title='shift_abs: ',xsize=3,/ALL_EVENTS,/FLOAT )
 ;tmp=widget_label(base_2,value=' ',ysize=7)

endif

if k eq 8 then begin
 bgroup = cw_bgroup(base_1,'correction depolarization',/nonEXCLUSIVE,uvalue='k8',set_value=key1(k),ypad=ypd)
 shiftPA_LAB=widget_label(base_1,value='shift PA',/align_center)
    shiftPA_W=WIDGET_TEXT(BASE_1,VALUE=string(x_slit),uvalue='shift_PA',xsize=5,/edit,/all,/align_center)
 empty=widget_label(base_1,value='deg   ')
 bgroup = cw_bgroup(base_1,'plot',/nonEXCLUSIVE,uvalue='k_8',set_value=key2(k),ypad=ypd)
 ;bgroup = cw_bgroup(base_1, 'mJy ', /nonEXCLUSIVE,  uvalue='k8',set_value=key1(k-1))
 ;bgroup = cw_bgroup(base_1, 'erg ', /nonEXCLUSIVE,  uvalue='k_8',set_value=key2(k-1))

output_range_lab=widget_label(base_1,value='output range, AA',/align_center)
out_beg_wave_W=WIDGET_TEXT(BASE_1,VALUE=string(0,format='(I4)'),uvalue='out_beg',xsize=5,/edit,/all,/align_center)
out_end_wave_W=WIDGET_TEXT(BASE_1,VALUE=string(0,format='(I4)'),uvalue='out_end',xsize=5,/edit,/all,/align_center)
width_lab=widget_label(base_1,value='integration window,px',/align_center)
width_integr=WIDGET_TEXT(BASE_1,VALUE=string(10,format='(I4)'),uvalue='integr',xsize=5,/edit,/all,/align_center)

endif

if k eq 9 then begin
tmp=WIDGET_LABEL(BASE_1,value='P(ISM),%',xsize=x1*4,/align_center)
  P_ISM_W=WIDGET_TEXT(BASE_1,VALUE=string(0,format='(I4)'),uvalue='P_ISM',xsize=5,/edit,/all,/align_center)
  tmp=WIDGET_LABEL(BASE_1,value='PA(ISM),deg',xsize=x1*4,/align_center)
  PA_ISM_W=WIDGET_TEXT(BASE_1,VALUE=string(0,format='(I4)'),uvalue='P_ISM',xsize=5,/edit,/all,/align_center)

endif

endfor

base2=WIDGET_BASE(mainbase,row=1,frame=1)
BUT_RUN=WIDGET_BUTTON(BASE2,VALUE='Run program',uvalue='RUN',xsize=x1*18)
tmp=WIDGET_LABEL(BASE2,value='configuration',xsize=x1*9,/align_center)
BUT_DEFAULT=WIDGET_BUTTON(BASE2,VALUE='default',uvalue='dflt')
BUT_LOAD=WIDGET_BUTTON(BASE2,VALUE='load',uvalue='load')
BUT_SAVE=WIDGET_BUTTON(BASE2,VALUE='save',uvalue='save')
   ; ***** EXIT and HELP *****
   EXITBASE=WIDGET_BASE(mainbase,row=1,frame=1)
   XLOADBUT=WIDGET_BUTTON(exitbase,value='XLOADCT',UVALUE='XLOA')
   HELPBUT=WIDGET_BUTTON(exitbase,value='HELP',UVALUE='HELP')
   EXITBUT=WIDGET_BUTTON(exitbase,value='EXIT',UVALUE='EXIT')
    WIDGET_CONTROL, MAINBASE, /REALIZE,group_leader=mainbase
        ;    0           1             2            3                 4        5        6       7
BUTTONS=[edit_LOG,Yshift,BUT_DEL_FTS,BUT_CREATE_datacube,BUT_BINX,BUT_BINY,BUT_Neta,BUT_VIEW,$
        ;   8       9          10       11     12  13   14    15     16     17     18    19      20	       21		22
  BUT_RUN,BUT_DEFAULT,BUT_LOAD,BUT_SAVE,Xbeg,wx,Ntra,DegTra,lines_atm,reg,Y_pos,Y_width,window_atm];,label_pos,obj_pos,label_w,obj_w,radiusX,radiusY,$
  ;    23      24            25                26        27       28
 ;	shiftPA_W,out_beg_wave_W,out_end_wave_W,width_integr,DEG_DISP];,P_ISM_W,PA_ISM_W]
for k=0,11 do WIDGET_CONTROL,BUTTONS(k),SENSITIVE=0
;for k=15,20  do WIDGET_CONTROL,BUTTONS(k),SENSITIVE=0
WIDGET_CONTROL,BUTTONS(6),SET_DROPLIST_SELECT=ind_Neta
WIDGET_CONTROL,BUTTONS(17),sensitive=0

  XMANAGER, 'WOLL2', MAINBASE
END

COMMON SET_step,step
COMMON obj_pos,objpos,overscan,overscan_val,bias,bias_val,width,yrange_val
 COMMON SET_key,key1,key2,key3,key_1,key_2,key_3
COMMON LOG,dLOG,text_log,log_dir,LOGFILE,idl_dir,gview
COMMON set_COL,background,color
COMMON PAR,par2,par5,par6_1,par6_2,par11_1,par11_2,par11_3

if !version.os_family eq 'unix' then $
idl_dir='/home/elias/SCORPIO/sppol_pipeline_v2023.8/WOLLASTON-2.lib' else $
idl_dir='h:\WOLLASTON-2.lib\'
;if os_family() eq 'unix' then $
;;log_dir='/db1/red_data/scorpio.log/' else $
;log_dir='h:\spectraPOL.log\'
if !version.os_family eq 'unix' then $
gview='gv ' else $
gview='c:\Program Files\Ghostgum\gsview\gsview64.exe '

;!P.background=16777215 &   !P.color=0
step=make_array(11,value=0,/integer)
par11_1='5' & par11_2=' ' & par11_3=' '
key1=[0,1,0,0,0,0,0,0,0,0,0]
key2=[0,0,0,0,0,0,0,0,0,0,0]
key3=[0,0,0,0,0,0,0,1,0,0,0]
objpos=0
WOLL2
print,step
print,key1
print,key2
end
