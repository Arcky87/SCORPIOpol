pro DRAWPLOT,OPLOT=oplot
COMMON IMAblock,file,h,frame,ima,LOGFILE,Nx,Ny,Nz,Nt,Naxis,sc,bin,plane,target,d_lambda,lambda_0,over,x_size,y_size,range,level
COMMON WIN,num_win,s_win,xs,ys,xgau,sw_gau,x_pre,y_pre
COMMON WINblock,win_num,swin,planeW,draw,tra,ind_tra,planeT
COMMON BOXblock,map,x1,x2,y1,y2,x_2,y_2,sw
COMMON REG,x_min,y_min,x_max,y_max,region,x,y,a,values
WSET, num_win

imascale=sxpar(h,'IMSCALE')
imascale=str_sep(imascale,'x')
imascale=float(imascale(1))
;imascale=0.178
limit=[x_min,y_min,x_max,y_max]
;!P.color=0 & !P.background=2^24-1
R=FLOAT(x_max-x_min)/float(y_max-y_min)
if R ge 1 then begin

S=x & index=1 & k=2
endif
if R lt 1 then begin

S=y & index=2 & k=1
endif
x_tmp=FIX((limit(index-1)+(limit(index+1)-limit(index-1))/530.*xs))
x_tmp_pre=FIX((limit(index-1)+(limit(index+1)-limit(index-1))/530.*x_pre))
if a(0) eq 2 then begin
vector=total(region,k)/a(k)
;if not(keyword_set(oplot)) then $

plot,S,vector,XST=1,yst=1,position=[0,0,1,1],xmin=1,ymin=1,yrange=[min(vector)-(max(vector)-min(vector))/20.,max(vector)],$
		xtickinterval=1E5,ytickinterval=1E8,psym=10,color=0
;индикация положения  траектории
if index eq 2 and k eq 1 and ind_tra eq 1 then begin
for j=0,10 do begin
oplot,[1,1]*tra((x_max+x_min)/2,j,plane),[min(vector),max(vector)],linestyle=1
xyouts,tra((x_max+x_min)/2,j,plane),min(vector)-(max(vector)-min(vector))/25.,string(tra((x_max+x_min)/2,j,plane),format='(F5.1)'),/data,align=0.5,color=0
endfor
endif
;print,tra((x_max-x_min)/2,*,plane),format='(11F7.1)'
if x_tmp eq limit(index-1) THEN col=2^24-1 ELSE col=0
oplot,x_tmp*[1,1],color=col,$
[min(vector),max(vector)],linestyle=2
if x_tmp_pre ne x_tmp then $
oplot,x_tmp_pre*[1,1],color=2^24-1,$
[min(vector),max(vector)],linestyle=2
oplot,S,vector,color=0,psym=10
if sw_gau ne 1 then begin
str_out=strcompress('pixel '+string(x_tmp))
if d_lambda ne 01 then str_out=str_out+'   '+strcompress('wavelength '+string(FIX(d_lambda*x_tmp+lambda_0))+' A')

WIDGET_CONTROL,values,SET_VALUE=str_out+ $
'    value '+string(vector(x_tmp-limit(index-1)))
endif
if sw_gau eq 1 then begin
;определение уровня фона
fon=(max(vector)-min(vector))*ys/300.+min(vector)
tmp_min=min(vector)-(max(vector)-min(vector))/20.
fon=(max(vector)-tmp_min)*ys/300.+tmp_min
res=MULTIGAUS ( S, vector-fon, x_tmp,FWHM=10,YFIT=FIT)
MM=max(abs(fit),Nmax)

ww=res.FWHM*3
wmin=Nmax-ww & if wmin lt 0 then wmin=0
wmax=Nmax+ww & if wmax gt N_elements(s)-1 then wmax=N_elements(s)-1
FIT=gaussian(S(wmin:wmax),[res.max,res.center,res.FWHM/2.345])
oplot,S(wmin:wmax),fit+fon,color=1e5
PRINT,SIZE(fit)
lamb_gau=''
FWHM_gau=''
flux_gau=''
if index eq 1 and d_lambda ne 0 then begin
lamb_gau=' ('+STRCOMPRESS(string(res.center*d_lambda+lambda_0,format='(F7.1)'),/remove_all)+' A)'
FWHM_gau=' ('+STRCOMPRESS(string(res.FWHM*d_lambda,format='(F7.1)'),/remove_all)+' A)'
flux_gau='  flux '+STRCOMPRESS(string(res.flux,format='(I9)'),/remove_all)
endif
if index eq 2 then see='  seeing '+$
	string(res.FWHM*imascale,format='(F7.1)')+'"' ELSE see=''
WIDGET_CONTROL,values,SET_VALUE='pixel'+string(x_tmp,format='(I5)')+$
'   value '+STRCOMPRESS(string(vector(x_tmp-limit(index-1))))+$
'   center'+string(res.center,format='(F7.1)')+'px '+lamb_gau+$
'   FWHM '+STRCOMPRESS(string(res.FWHM,format='(F7.1)'))+' px '+FWHM_gau+flux_gau+see

endif
	endif
	!P.color=0 & !P.background=2^24-1
;!P.color=255 & !P.background=0
END




Pro Wplot_event, event
COMMON Wbase,base_main,base_ccd,base_platform,base_scorpio,base_view,base_plot
COMMON IMAblock,file,h,frame,ima,LOGFILE,Nx,Ny,Nz,Nt,sc,bin,plane,target,d_lambda,lambda_0,over,x_size,y_size,range,level
COMMON WIN,num_win,s_win,xs,ys,xgau,sw_gau,x_pre,y_pre
COMMON WINblock,win_num,swin,planeW,draw,tra,ind_tra,planeT
COMMON BOXblock,map,x1,x2,y1,y2,x_2,y_2,sw
COMMON REG,x_min,y_min,x_max,y_max,region,x,y,a,values
; When a widget is selected, put its User Value into 'eventval':

WIDGET_CONTROL, event.id, GET_UVALUE = eventval

; Perform actions based on the User Value of the event:

CASE eventval OF
   'DRAW_EVENT': BEGIN

			IF (event.press EQ 0) AND (event.release EQ 0)   THEN BEGIN

			xs=event.x & ys=event.y

			if xs lt 520 then DRAWPLOT,/oplot
			x_pre=xs & y_pre=xs
			ENDIF
			sw_gau=0
			if event.press EQ 1 then begin
			sw_gau=1 & xgau=event.x
			DRAWPLOT,/oplot
			endif


		     END
   'DONE': WIDGET_CONTROL, event.top, /DESTROY
ENDCASE

END



PRO Wplot, GROUP=GROUP
COMMON Wbase,base_main,base_ccd,base_platform,base_scorpio,base_view,base_plot
COMMON IMAblock,file,h,frame,ima,LOGFILE,Nx,Ny,Nz,Nt,sc,bin,plane,target,d_lambda,lambda_0,over,x_size,y_size,range,level
COMMON WIN,num_win,s_win,xs,ys,xgau,sw_gau,x_pre,y_pre
COMMON WINblock,win_num,swin,planeW,draw,tra,ind_tra,planeT
COMMON BOXblock,map,x1,x2,y1,y2,x_2,y_2,sw
COMMON REG,x_min,y_min,x_max,y_max,region,x,y,a,values
s_win = !D.WINDOW	; Remember the current window so it can be restored
sw_gau=0
xs=0 & ys=0
x_pre=1 & y_pre=1
; This example uses keywords to define the size of the draw area:
base_plot = WIDGET_BASE(TITLE = 'Plot region',/COLUMN,xoffset=720,yoffset=490)

WIDGET_CONTROL, /MANAGED, base_plot


;
draw_plot = WIDGET_DRAW(base_plot , XSIZE=530, YSIZE=300,/FRAME, $
	/MOTION_EVENTS,  /BUTTON_EVENTS ,$		;generate LOTS of events
	UVALUE = 'DRAW_EVENT', 	RETAIN = 2)


base_inform=WIDGET_BASE(base_plot,/row,xsize=530,ysize=30)
values=WIDGET_LABEL(base_inform,value=' ',xsize=465)
button = WIDGET_BUTTON(base_inform , $
		UVALUE = 'DONE', $
		VALUE = 'DONE')

WIDGET_CONTROL, base_plot , /REALIZE



WIDGET_CONTROL, draw_plot, GET_VALUE=num_win


DRAWPLOT
; Hand off control of the widget to the XMANAGER:
XMANAGER, "Wplot", base_plot , GROUP_LEADER=GROUP, /NO_BLOCK

END
