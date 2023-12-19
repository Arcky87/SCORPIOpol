pro  LOADFILE,REP=rep,DIR=dir,YSZ=ysz,XSZ=xsz
COMMON IMAblock,file,h,frame,ima,WDIR,Nx,Ny,Nz,Nt,Naxis,sc,bin,plane,target,d_lambda,lambda_0,over,x_size,y_size,range,level,avg_ima,rms_ima
COMMON WINblock,win_num,swin,planeW,draw,tra,ind_tra,planeT
if keyword_set(dir) then WDIR=dir
if keyword_set(xsz) then x_size=xsz
if keyword_set(ysz) then y_size=ysz
file=DIALOG_pickfile(PATH=wdir,/read,filter='*.fts')
tmp_frame=readfits(file,h)
Naxis=sxpar(h,'NAXIS')
Nx=sxpar(h,'NAXIS1')
N_y=sxpar(h,'NAXIS2')
Ny=N_y*4
print,'Ny',Ny

Npol=sxpar(h,'NAXIS3')
sc=float(y_size)/float(Ny)
d_lambda=sxpar(h,'CDELT1')
lambda_0=sxpar(h,'CRVAL1')
type=string(Naxis,format='(I1)')
print,x_size,y_size
CASE type OF
'3':	begin
		frame=fltarr(Nx,Ny)
		for k=0,Npol-1 do frame(*,k*N_y:N_y*(k+1)-1)=tmp_frame(*,*,k)
		ima=congrid(frame(*,*),x_size,y_size)
		Nz=0
				end

'5':	begin
		;definition number  exposure of targets
		Nobj=sxpar(h,'NAXIS5')  & Nz=Nobj
		Nt=fltarr(Nobj)
		for k=0,Nobj-1 do Nt(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))

		Nexp=sxpar(h,'NAXIS4')
		frame=fltarr(Nx,Ny,Nexp,Nobj)
		ima=fltarr(x_size,y_size,Nexp,Nz)
		for k=0,Npol-1 do frame(*,k*N_y:N_y*(k+1)-1,*,*)=tmp_frame(*,*,k,*,*)
		for j=0,Nz-1 do ima(*,*,*,j)=congrid(frame(*,*,*,j),x_size,y_size,Nexp)

		end

ENDCASE



;bin=1
;XY=float(Nx)/float(Ny)
;sc=float(y_size)/float(Ny)
;if bin eq 0 then begin
;bin=1
;goto, cont1
;endif
;bin=fix(st;r_sep(bin,'x')) & bin=float(bin(0))/bin(1)
;cont1:
;
;
;x_size=y_size*XY*float(bin)
;if keyword_set(rep) then begin
;WIDGET_CONTROL, draw, GET_VALUE=win_num,DRAW_XSIZE=x_size
;WIDGET_CONTROL, draw, GET_VALUE=win_num,DRAW_YSIZE=y_size
;WIDGET_CONTROL, draw, GET_VALUE=win_num,SET_DRAW_VIEW=[(x_size-y_size)/2,0]
;print,x_size,y_size
;endif

;ima=congrid(frame,x_size,y_size,Nz)

;endif

;if Naxis eq 3 then begin

;tmp=fltarr(Nx,Ny*4)
;for k=0,3 do  tmp(*,k*Ny:Ny*(k+1)-1)=frame(*,*,k)

;ima=congrid(tmp(*,*),x_size,y_size)
;endif


rms_ima=stdev(ima(*,*,0,0),avg_ima)


;чтение файла траектории

;if Nz ne 0 then begin
;file_tra=STRMID(file,0,strlen(file)-strlen(FILE_BASENAME(FILE)))+'tra.fit'
;if FILE_TEST(file_tra) eq 0 then tra=0 ELSE tra=readfits(file_tra)
;endif
;Ny=4*Ny
END

pro REGIONparam
COMMON IMAblock,file,h,frame,ima,WDIR,Nx,Ny,Nz,Nt,Naxis,sc,bin,plane,target,d_lambda,lambda_0,over,x_size,y_size,range,level,avg_ima,rms_ima
COMMON REG,x_min,y_min,x_max,y_max,region,x,y,a,values
COMMON BOXblock,map,x1,x2,y1,y2,x_2,y_2,sw,base
; calculation region coordinates

x_min=min([x1,x2])*float(Nx)/float(x_size) & x_max=max([x1,x2])*float(Nx)/float(x_size)
y_min=min([y1,y2])*float(Ny)/float(y_size) & y_max=max([y1,y2])*float(Ny)/float(y_size)
if x_min lt 0 then x_min=0
if y_min lt 0 then y_min=0
if x_max gt Nx-1 then x_max=Nx-1
if y_max gt Ny-1 then y_max=Ny-1
region=frame(x_min:x_max,y_min:y_max,plane,target)
a=size(region)
x=findgen(a(1))+x_min
y=findgen(a(2))+y_min
END


pro BOX
COMMON IMAblock,file,h,frame,ima,WDIR,Nx,Ny,Nz,Nt,Naxis,sc,bin,plane,target,d_lambda,lambda_0,over,x_size,y_size,range,level,avg_ima,rms_ima
COMMON WINblock,win_num,swin,planeW,draw,tra,ind_tra,planeT
COMMON BOXblock,map,x1,x2,y1,y2,x_2,y_2,sw,base
WSET, win_num

TV,map
plot,[0,x_size],[0,y_size],/noerase,/nodata,position=[0,0,1,1],$
	xtickinterval=x_size,ytickinterval=y_size,xst=1,yst=1,$
	xminor=1,yminor=1,/norm
;if file eq 'map_i' then oplot,[1,1]*slitpos,[0,y_size],color=255
oplot,[x1,x_2,x_2,x1,x1],[y1,y1,y_2,y_2,y1],color=255
for k=0,2 do oplot,[0,x_size-1],[1,1]*y_size/4*(k+1),color=255
end


pro DRAWIMAGE
COMMON IMAblock,file,h,frame,ima,WDIR,Nx,Ny,Nz,Nt,Naxis,sc,bin,plane,target,d_lambda,lambda_0,over,x_size,y_size,range,level,avg_ima,rms_ima
COMMON WINblock,win_num,swin,planeW,draw,tra,ind_tra,planeT
COMMON BOXblock,map,x1,x2,y1,y2,x_2,y_2,sw,base
WSET, win_num

;range=50 & level=50
print,range,level
i=50./range
j=50./level

plane_value=['obj','unpolarized star','polarized star']
WIDGET_CONTROL,base,BASE_SET_TITLE='view frame '+file
if Nz eq 0 then begin
WIDGET_CONTROL,planeW,SENSITIVE=0
plane=0
endif
if Naxis eq 3 then begin
WIDGET_CONTROL,planeT,SENSITIVE=0
target=0
endif
;I_lim=[avg_ima-rms_ima, avg_ima+5*rms_ima]
I_lim=[avg_ima*j-rms_ima*i,avg_ima*j+3*i*rms_ima]
;if FILE_BASENAME(file) eq 'flat.fts' then I_lim=[0.5,1.5]
map=255-bytscl(ima(*,*),I_lim(0),I_lim(1))

if Nz gt 0 then begin
WIDGET_CONTROL,planeW,set_value=plane_value
WIDGET_CONTROL,planeW,SET_DROPLIST_SELECT=plane
WIDGET_CONTROL,planeW,SENSITIVE=1
endif
if Naxis eq 5 then begin
Nexp=Nt(plane)
;target_value=['obj total exposure','unpolarized star','polarized star']
target_value='exposure'+string(findgen(Nt(0))+1,format='(I3)')
;print,target_value
for k=0,Nexp-1 do target_value=[target_value,'obj exposure'+string(k+1,format='(I2)')]

WIDGET_CONTROL,planeT,set_value=target_value(0:Nt(plane)-1)
WIDGET_CONTROL,planeT,SET_DROPLIST_SELECT=target
WIDGET_CONTROL,planeT,SENSITIVE=1
;I_lim=[avg_ima-rms_ima, avg_ima+5*rms_ima]
I_lim=[avg_ima*j-rms_ima*i,avg_ima*j+3*i*rms_ima]
;if FILE_BASENAME(file) eq 'flat.fts' then I_lim=[0.5,1.5]
map=255-bytscl(ima(*,*,target,plane),I_lim(0),I_lim(1))


endif





TV,map
plot,[0,x_size-1],[0,y_size-1],/noerase,/nodata,position=[0,0,1,1],$
	xtickinterval=x_size,ytickinterval=y_size,xst=1,yst=1,$
	xminor=1,yminor=1,/norm
for k=0,2 do oplot,[0,x_size-1],[1,1]*y_size/4*(k+1),color=255


if ind_tra eq 1 AND Nz gt 0 then begin
		plot,[0,Nx-1],[0,Ny-1],/nodata,/noerase,xst=1,yst=1,position=[0,0,1,1],color=0,ticklen=0

		a=size(tra)
		for j=0,a(2)-1 do oplot,tra(*,j,plane),color=0,linestyle=2
		endif

;WSET, swin			; Restore the original window
END


PRO ViewPol_2_event, event

; This is the event handler for a draw widget which tracks
; cursor motion events.

; The COMMON block is used because the event handler needs
; widget ID's of the labels:
COMMON Wbase,base_main,base_ccd,base_platform,base_scorpio,base_view,base_plot
COMMON LABELblock, pixel_label, value_label,mean_label,rms_label,POS_label,FWHM_label,SEE_label,wave_label
COMMON SLIDEdlock,viewrange,viewbackground
COMMON REG,x_min,y_min,x_max,y_max,region,x,y,a,values
COMMON IMAblock,file,h,frame,ima,WDIR,Nx,Ny,Nz,Nt,Naxis,sc,bin,plane,target,d_lambda,lambda_0,over,x_size,y_size,range,level,avg_ima,rms_ima
COMMON WINblock,win_num,swin,planeW,draw,tra,ind_tra,planeT
COMMON BOXblock,map,x1,x2,y1,y2,x_2,y_2,sw,base
; When a widget is selected, put its User Value into 'eventval':

WIDGET_CONTROL, event.id, GET_UVALUE = eventval

; Perform actions based on the User Value of the event:

CASE eventval OF
'plane'      :  BEGIN

				plane=event.index
				rms_ima=stdev(ima(*,*,0,plane),avg_ima)
				DRAWIMAGE
				END
'target'     : BEGIN
				target=event.index
				;print,'target=',target
				;ima=congrid(frame(*,*,*,target),x_size,y_size,Nz)
				DRAWIMAGE
				END

'tra':	BEGIN
		WIDGET_CONTROL,event.id,GET_VALUE=ind_tra
		DRAWIMAGE
		END
'DRAW_WIN_EVENT': BEGIN
			;HELP, /STRUCT, event

            if event.press EQ 1  then begin
            sw=1
            x1=event.x  & y1=event.y

			endif

            if event.release EQ 1 and event.clicks EQ 1 then begin
            sw=0
            x2=event.x  & y2=event.y

			REGIONparam
			DRAWIMAGE
			BOX
            endif

            IF (event.press EQ 0) AND (event.release EQ 0) THEN BEGIN
            x_2=event.x  & y_2=event.y
			if sw eq 1 then  BOX
				WIDGET_CONTROL,pixel_label, $
				   SET_VALUE='Pixel' + STRCOMPRESS('('+STRING(FIX(event.x*float(Nx)/float(x_size)))+','+$
				   STRING(FIX(event.y*float(Ny)/float(y_size)))+')')
				   if d_lambda gt 0 then $
				WIDGET_CONTROL,wave_label, $
				   SET_VALUE='wavelength' + STRCOMPRESS(STRING(FIX((event.x+1)*float(Nx)*d_lambda/float(x_size)+lambda_0))+'A') $
				   else WIDGET_CONTROL,wave_label,SET_VALUE=' '

				if Naxis  eq 3 and event.x ge 0 and event.x le x_size-1 and event.y ge 0 and event.y le y_size-1 then  $
				WIDGET_CONTROL,value_label,$
				SET_VALUE='Value=' + STRING(ima(event.x,event.y)); ELSE $
				;WIDGET_CONTROL,value_label,$
			;	SET_VALUE='Value=' + STRING(ima(event.x,event.y,plane))
			ENDIF
			if event.press eq 4 then $
			 if not (x1 eq x2 and y1 eq y2) then begin
			 rms_ima=stdev(median(region,3),avg_ima)
			  DRAWIMAGE
			endif

		    END
'range': BEGIN
			 WIDGET_CONTROL, viewrange, GET_VALUE=range
			 DRAWIMAGE
			 END
'background': BEGIN
			 WIDGET_CONTROL, viewbackground, GET_VALUE=level

			 DRAWIMAGE
			 END
'plot_region': BEGIN
			if XREGISTERED('Wplot' ) gt 0 then  WIDGET_CONTROL,base_plot,/destroy
			Wplot
			end
'region_statistic':BEGIN
		robomean,region,3,0.5,avg_region,rms_region
		;rms_region=stdev(region,avg_region)
		WIDGET_CONTROL,mean_label,SET_VALUE='mean ' + STRING(avg_region)
		WIDGET_CONTROL,rms_label,SET_VALUE='rms ' + STRING(rms_region)
		end
'region_fwhm': BEGIN
			res=GAUSS2DFIT(region,G)
			FWHM=(sqrt(G(2)^2+G(3)^2))*2.345
			WIDGET_CONTROL,POS_label,SET_VALUE='x=' + STRING(G(4)+x_min,format='(F6.1)')+' y='+ STRING(G(5)+y_min,format='(F6.1)')
			WIDGET_CONTROL,FWHM_label,SET_VALUE='FWHM='+STRING(FWHM,format='(F6.1)')+' px'
			if file eq 'map_i' then begin
			imascale=sxpar(h,'IMSCALE')
			imascale=str_sep(imascale,'x') & imascale=float(imascale(1))
			WIDGET_CONTROL,SEE_label,SET_VALUE='seeing='+STRING(FWHM*imascale,format='(F4.1)')+'"'
			endif
				END
'FILE':     BEGIN
			LOADFILE,/rep
			DRAWIMAGE
			END
   'DONE': BEGIN
   			if XREGISTERED('Wplot' ) gt 0 then  WIDGET_CONTROL,base_plot,/destroy
   			WIDGET_CONTROL, event.top, /DESTROY
   			end
ENDCASE

END



PRO ViewPol_2, GROUP=GROUP

; This is the procedure that creates a draw widget that returns
; motion events.

; The COMMON block is used because the event handler needs
; widget ID's of the labels:
COMMON Wbase,base_main,base_ccd,base_platform,base_scorpio,base_view,base_plot
COMMON LABELblock, pixel_label, value_label,mean_label,rms_label,POS_label,FWHM_label,SEE_label,wave_label
COMMON SLIDEdlock,viewrange,viewbackground
COMMON IMAblock,file,h,frame,ima,WDIR,Nx,Ny,Nz,Nt,Naxis,sc,bin,plane,target,d_lambda,lambda_0,over,x_size,y_size,range,level,avg_ima,rms_ima
COMMON WINblock,win_num,swin,planeW,draw,tra,ind_tra,planeT
COMMON BOXblock,map,x1,x2,y1,y2,x_2,y_2,sw,base
ind_tra=0
Y_size=500  & X_size=1130
swin = !D.WINDOW	; Remember the current window so it can be restored
overscan=20
; This example uses keywords to define the size of the draw area:
range=50 & level=50
plane=0 ;set plane 1
target=0; set target obj total
sw=0 & x1=0& y1=0 & x2=0 & y2=0 & x_2=0 &y_2=0
; A top-level base widget with the title "Motion Event Widget Example"
; will be created.  The size is left unspecified until the draw widget
; is created:

base = WIDGET_BASE(TITLE = 'View frame',/ROW); $
	;/COLUMN)
BASE1 = WIDGET_BASE(base,col=1,MAP=1)
; Setting the managed attribute indicates our intention to put this application
; under the control of XMANAGER, and prevents our draw widgets from
; becoming candidates for becoming the default window on WSET, -1. XMANAGER
; sets this, but doing it here prevents our own WSETs at startup from
; having that problem.
WIDGET_CONTROL, /MANAGED, base


; A widget called 'draw' is created:

;draw = WIDGET_DRAW(base1, X_SCROLL_SIZE=y_size, Y_SCROLL_SIZE=y_size,$
 ;   XSIZE=x_size, YSIZE=y_size, /SCROLL,/FRAME, $

draw = WIDGET_DRAW(base1, XSIZE=x_size, YSIZE=y_size,/FRAME, $
	/MOTION_EVENTS, /BUTTON_EVENTS ,$		;generate LOTS of events
	UVALUE = 'DRAW_WIN_EVENT', 	RETAIN = 1)

BASE2 = WIDGET_BASE(base1,row=1,MAP=1)

base3=WIDGET_BASE(base,col=1,MAP=1)
button1 = WIDGET_BUTTON(base3, $
		UVALUE = 'FILE', $
		VALUE = 'LOAD FILE')
plane_label = WIDGET_LABEL(base3,  $
	VALUE='set target')
planeW=WIDGET_DROPLIST(base3,value=['object','unpolarized star','polarized star'],uvalue='plane',xsize=120)
target_label = WIDGET_LABEL(base3,  $
	VALUE='set target')
planeT=WIDGET_DROPLIST(base3,value=['obj total','unpolarized star','polarized star'],uvalue='target',xsize=120)
traW=CW_BGROUP(base3,['show traectory'],/NONEXCLUSIVE,uvalue='tra',set_value=0,/row,/return_index)
viewrange = WIDGET_SLIDER(base3,$ TITLE = 'contrast', $
                        MINIMUM = 1, MAXIMUM = 100, VALUE = 50, $
			UVALUE = 'range')
contrast_label = WIDGET_LABEL(base3,  $
	VALUE='contrast')
viewbackground = WIDGET_SLIDER(base3,$ TITLE = 'brigthness', $
                        MINIMUM = 1, MAXIMUM = 100, VALUE = 50, $
			UVALUE = 'background')
background_label = WIDGET_LABEL(base3,  $
	VALUE='background')
label1=WIDGET_LABEL(base3,  $
	VALUE='  ')
pixel_label = WIDGET_LABEL(base3, /DYNAMIC_RESIZE, $
		VALUE='')

value_label = WIDGET_LABEL(base3, /DYNAMIC_RESIZE, $
		VALUE='')
wave_label=WIDGET_LABEL(base3, /DYNAMIC_RESIZE, $
		VALUE=' ')
label2=WIDGET_LABEL(base3,  $
	VALUE='  ')
button3 = WIDGET_BUTTON(base3, $
		UVALUE = 'plot_region', $
		VALUE = 'plot region')
label3=WIDGET_LABEL(base3,  $
	VALUE='  ')
button4 = WIDGET_BUTTON(base3, $
		UVALUE = 'region_statistic', $
		VALUE = 'region statistic')
mean_label = WIDGET_LABEL(base3, /DYNAMIC_RESIZE, $
		VALUE='mean value')
rms_label = WIDGET_LABEL(base3, /DYNAMIC_RESIZE, $
		VALUE='rms value')

button5 = WIDGET_BUTTON(base3, $
		UVALUE = 'region_fwhm', $
		VALUE = 'region centroid')
POS_label = WIDGET_LABEL(base3, /DYNAMIC_RESIZE, $
		VALUE='x    y')

FWHM_label = WIDGET_LABEL(base3, /DYNAMIC_RESIZE, $
		VALUE='FWHM')
SEE_label=WIDGET_LABEL(base3, /DYNAMIC_RESIZE, $
		VALUE=' ')

button2 = WIDGET_BUTTON(base3, $
		UVALUE = 'DONE', $
		VALUE = 'DONE')
; Realize the widgets:


WIDGET_CONTROL, base, /REALIZE

; Get the window number from the draw widget.  This can only be done
; after the widget has been realized:

WIDGET_CONTROL, draw, GET_VALUE=win_num,SET_DRAW_VIEW=[(x_size-y_size)/2,0]
DRAWIMAGE
; Use TVSCL to display an image in the draw widget.  Set the window for
; the TVSCL command since there may be other draw windows.

;WIDGET_CONTROL,planeW,SENSITIVE=0
; Hand off control of the widget to the XMANAGER:
XMANAGER, "ViewPol_2", base, GROUP_LEADER=GROUP;, /NO_BLOCK
END

;FILELOG=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path='h:\red_data.pol\AGN\LOGS')
;W_DIR=sxpar(read_table(FILELOG),'w_dir')
log_dir='h:\red_data.pol\AGN\'
LOGFILE=DIALOG_PICKFILE(/read,path=log_dir+'LOGS\',FILTER='*.txt')
;LOGFILE=DIALOG_PICKFILE(/read,path='h:\red_data.pol\LOGS\',FILTER='*.txt')
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
				wdir=log_dir+wdir(N_elements(wdir)-2)+'\'
				print,wdir
LOADFILE,dir=wdir,ysz=500,xsz=1130

ViewPol_2
END