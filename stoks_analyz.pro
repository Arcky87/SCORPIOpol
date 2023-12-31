function MOD_value,value
value=value+1
value=value-2*(value/2)
return,value
end


pro DRAWPLOT,OPLOT=oplot
COMMON Wbase,base_plot
COMMON GRAPH,spectra,stoks,wave,dir,z
COMMON SET_widget,BUTTONS
COMMON SPLOT,stoks_in,stoks_out,xrng,plot_ind,plot_items,target,sent_ind,atm_ind,ps
COMMON OPTION,wx,thresh,PA,PA_null,PA_rot,Q_ISM,U_ISM,err,corr

!P.multi=[0,1]
plot,[0,1],[0,1],position=[0,0,1,1],/norm,/nodata

;multiplot;,NUM=num
END

pro multiplot
COMMON GRAPH,spectra,stoks,wave,dir,z
COMMON FILES,LOG_dir,LOGFILE,pro_dir,wdir,woll,head,output_file
COMMON SPLOT,stoks_in,stoks_out,xrng,plot_ind,plot_items,target,sent_ind,range_Stoks,range_Angle
COMMON OPTION,wx,thresh,PA,PA_null,PA_rot,Q_ISM,U_ISM,err,corr

num=plot_ind
print,'plot_ind',plot_ind
RN=where(num eq 1,indN)
xrange=xrng


if not(keyword_set(err)) then err=0
if not(keyword_set(target))then T=0 ELSE T=TARGET
S=N_elements(wave)
R=where(spectra.data(*,0) gt 0, F)
if keyword_set(xrange) then $
R= where(spectra.data(*,0) gt xrange(0) AND spectra.data(*,0) lt xrange(1))
;
!P.charsize=0.8
if not(keyword_set(Nplot)) then Nplot=5
if not(keyword_set(xtitle)) then xtitle=''
if not(keyword_set(ytitle)) then ytitle=strarr(Nplot)

X_margin=[0.1,0.995]
Y_margin=[0.05,0.97]
dY=(Y_margin(1)-Y_margin(0))/Nplot

for k=0,Nplot-1 do begin
if k eq Nplot-1 then xchar=1 ELSE xchar=1e-5
if k eq 0 then noerase=0  ELSE noerase=1

if k eq 0 then begin
;����������� �������� ������ �� ��������
max_spectra=max(spectra.data(0:F-1,T+1))
order=RFIX(ALOG10(max_spectra))
 wavL=[4861,5200,6565]
plot,spectra.data(R,0),spectra.data(R-1,T+1)/10.^order, $
	yrange=[-0.02,1.2]*max_spectra/10.^order,yst=1,psym=10,$

	title=STRCOMPRESS(spectra.date+', object '+spectra.name(T)+', z='+string(z,format='(F7.4)')+$
	', exposure '+string(spectra.exptime(T)*spectra.Nexp(T),format='(I4)')+$
	' s, PA '+string(spectra.PA(T),format='(F6.1)')+' deg'),$
	noerase=noerase,xcharsize=xchar,xst=1,$
	position=[X_margin(0),Y_margin(1)-dY*(k+1),X_margin(1),Y_margin(1)-dY*k],/norm
xyouts,0.05,Y_margin(1)-dY*(k+1)+dY/2,'Flux, 10!U'+string(order-2,format='(I1)')+'!N, ADU/px',/norm,orientation=90,align=0.5
	if target eq 0 then for L=0,2 do oplot,[1,1]*(1+z)*wavL(L),[-0.02,1.2]*max_spectra/10.^order,linestyle=2
endif

if k gt 0 then begin
if strmid(stoks.ytitle(RN(k-1)),0,1) eq 'F' then begin
order=RFIX(ALOG10(stoks.yrange(1,k-1,T)))
stoks.ytitle(0)='Flux, 10!U'+string(order,format='(I1)')+'!N, ADU/px'
endif  ELSE order=0
plot,wave(0:S-1),stoks.value(0:S-1,RN(k-1),T)/10.^order,position=[X_margin(0),Y_margin(1)-dY*(k+1),X_margin(1),Y_margin(1)-dY*k],/norm,$
		noerase=noerase,xcharsize=xchar,xtitle=xtitle,xrange=xrange,xst=1,$
	yrange=stoks.yrange(*,RN(k-1),T)/10.^order,yst=1,psym=10
if target eq 0 then for L=0,2 do oplot,[1,1]*(1+z)*wavL(L),[-0.02,1.2]*max_spectra/10.^order,linestyle=2
if err gt 0  then OPLOTERR, wave(0:S-1),stoks.value(0:S-1,RN(k-1),T)/10.^order,stoks.error(0:S-1,RN(k-1),T)/10.^order,psym=10
xyouts,total(xrange)/2,(stoks.yrange(1,RN(k-1),T)-stoks.yrange(0,RN(k-1),T))*0.8+stoks.yrange(0,RN(k-1),T),stoks.param(RN(k-1),T),align=0.5
xyouts,0.05,Y_margin(1)-dY*(k+1)+dY/2,stoks.ytitle(RN(k-1)),/norm,orientation=90,align=0.5
;ytitle=stoks.ytitle(k-1)
endif

endfor
xyouts,total(X_margin)/2,0.02,'Wavelength, '+STRING(197B),/norm,align=0.5

label=dir+', '+woll+', '+string(wx,format='(I3)')+' px,'+string(thresh,format='(I2)')+' rms,'
if PA_rot ne 0 then label=label+' rotation'+string(PA_rot,format='(I4)')+' deg,'
if target eq 0 and woll eq 'WOLL2' then begin
if corr eq 1 then label=label+' 2order corrected,'
if Q_ISM+U_ISM ne 0  then label=label+' ISM(%)=['+string(Q_ISM,format='(F5.2)')+$
	','+string(U_ISM,format='(F5.2)')+']'
endif
xyouts,0,0.005,label,/norm
NN=N_elements(wave(0:S-1))
openw,1,dir+'test.txt'
for k=0,NN-1 do printf,1,wave(k),stoks.value(k,4,0)
close,1
;wave(0:S-1),stoks.value(0:S-1,4,0),xrange=[4500,8000],xst=1

fin:
end


;******************************************************
pro robust_estimate_stoks,ORDER2=order2
COMMON GRAPH,spectra,stoks,wave,dir,z
COMMON FILES,LOG_dir,LOGFILE,pro_dir,wdir,woll,head,output_file
COMMON SPLOT,stoks_in,stoks_out,xrng,plot_ind,plot_items,target,sent_ind,range_Stoks,range_Angle
COMMON OPTION,wx,thresh,PA,PA_null,PA_rot,Q_ISM,U_ISM,err,corr
dir=wdir


;******************************************************
set_plot,'WIN'

Nx=sxpar(head,'NAXIS1') & Ntarget=3
lambda=findgen(Nx)*sxpar(head,'CDELT1')+sxpar(head,'CRVAL1')

;��������A ����������  SPECTRA
Nmax=5000  & Ntarget=3
spectra={spectra,data:fltarr(Nmax,Ntarget+1),$
     		title:strarr(3),name:strarr(Ntarget),date:strarr(1),$
     		Nexp:intarr(3),exptime:fltarr(Ntarget),PA:fltarr(Ntarget)}  ; j=2
;��������A ����������  STOKS
Smax=2000  & Ntarget=3
stoks={stoks,value:fltarr(Smax,6,Ntarget),error:fltarr(Smax,6,Ntarget),yrange:fltarr(2,6,Ntarget),$
     		mean:fltarr(6,Ntarget),rms:fltarr(6,Ntarget),param:strarr(6,Ntarget),ytitle:strarr(6)}
stoks.ytitle=['Flux, ADU','Q-Stoks, %','U-stoks, %','V-Stoks','Degree polarization, %','Angle polarization, deg']
Npol=4
;������ ���������� ����������

spectra.data(0:Nx-1,0)=lambda
spectra.date=sxpar(head,'DATE-OBS')

for k=0,Ntarget-1 do begin
spectra.name(k)=sxpar(head,'NAME'+string(k+1,format='(I1)'))
spectra.exptime(k)=sxpar(head,'EXPTIME'+string(k+1,format='(I1)'))
spectra.PA(k)=sxpar(head,'PA'+string(k+1,format='(I1)'))
spectra.Nexp(k)=sxpar(head,'NUMEXP'+string(k+1,format='(I1)'))
endfor
PRINT,'WX=',WX
Npos=float(Nx)/wx-1
xpos=findgen(Npos)*wx+wx/2
wave=lambda(xpos)
V=fltarr(Ntarget)
CASE WOLL OF
'WOLL1':begin
		spectra.data(0:Nx-1,1:3)=stoks_in(*,0,0:2)
		Nexp=sxpar(head,'NAXIS3')-Ntarget
		V(*)=total(stoks_in(*,3,0:2),1)/Nx; & print,V
		for k=0,3 do begin
 		for j=0,Npos-1 do begin
 		for i=0,Ntarget-1 do begin
 		IF K EQ 3 and V(i) EQ 0 THEN GOTO,OUT
		if i eq 0 then robomean,stoks_in(xpos(j)-wx/2:xpos(j)+wx/2,k,3:Nexp-1),thresh,0.5,avg_val,rms_val  ELSE $
		robomean,stoks_in(xpos(j)-wx/2:xpos(j)+wx/2,k,i),thresh,0.5,avg_val,rms_val
		stoks.value(j,k,i)=avg_val
		stoks.error(j,k,i)=rms_val
		OUT:
		endfor & endfor & endfor
		;���������� ������� �����������
		stoks.value(*,4,*)=sqrt(stoks.value(*,1,*)^2+stoks.value(*,2,*)^2)
		stoks.error(*,4,*)=sqrt(stoks.error(*,1,*)^2+stoks.error(*,2,*)^2)/sqrt(2)
		stoks.value(*,0,*)=stoks.value(*,0,*)*stoks.value(*,4,*)
		stoks.error(*,0,*)=stoks.value(*,0,*)*stoks.error(*,4,*)
		;���������� ���� ��������� �����������
		;stoks.value(*,5,*)=calc_atan(stoks.value(*,1,*),stoks.value(*,2,*))/2
		stoks.value(*,5,*)=atan(stoks.value(*,2,*)/stoks.value(*,1,*))/2*180/!PI
		for k=0,2 do begin
		;RF=WHERE(stoks.value(*,5,k) gt 180, ind) & if ind gt 0 then stoks.value(RF,5,k)=stoks.value(RF,5,k)-180
		endfor
		stoks.error(*,5,*)=28.4*stoks.error(*,4,*)


		END
'WOLL2': begin
		for j=0,Ntarget-1 do begin
		for k=0,Nx-1 do begin
		robomean,stoks_in(k,0,0:spectra.Nexp(j)-1,j),3,0.5,avg_val,rms_val
		spectra.data(k,j+1)=median(stoks_in(k,0,0:spectra.Nexp(j)-1,j))
		;spectra.data(k,j+1)=avg_val
		endfor & endfor
;������ �������� � ����������������������� ���������
	 ; w;indow,2
	  	!P.multi=[0,1,4]
	  	null=TOTAL(stoks_in(*,*,0:spectra.Nexp(1)-1,1),3)/spectra.Nexp(1)
	  	for k=0,Npol-1 do begin
	  	plot,null(*,k),xst=1
	  	null(*,k)=LOWESS(findgen(Nx),null(*,k),Nx/8,2,3)
	  	oplot,null(*,k),thick=1,color=1e5
	  	endfor
		;���� ���������������� �����������
		for i=0,Ntarget-1 do begin
		for j=0,spectra.Nexp(i)-1 do begin
		for k=1,Npol-1 do stoks_in(*,k,j,i)=stoks_in(*,k,j,i)-null(*,k)
	endfor & endfor


;��������� ������ ���������� ������ � ����
		for k=0,Npol-1 do begin
 		for j=0,Npos-1 do begin
 		for i=0,Ntarget-1 do begin
 		IF k EQ 3 and V(i) EQ 0 THEN GOTO,OUT2
		robomean,stoks_in(xpos(j)-wx/2:xpos(j)+wx/2,k,0:spectra.Nexp(i)-1,i),thresh,0.5,avg_val,rms_val
		stoks.value(j,k,i)=avg_val
		stoks.error(j,k,i)=rms_val
		OUT2:
		endfor & endfor & endfor
			if keyword_set(corr) then begin
		;����������� 2-�� ������� � �������
		IF keyword_set(order2) then begin
		for k=1,2 do begin
		f=goodpoly(ALOG10(wave(0:Npos/2)),ALOG10(stoks.value(0:Npos/2,k,0)),1,3,FIT)
        fit=f(0)+ALOG10(wave)*f(1)
        fit=10.^fit
        stoks.value(0:Npos-1,k,0)=stoks.value(0:Npos-1,k,0)-LOWESS(findgen(Npos),stoks.value(0:Npos-1,k,0),Npos/4,2,3)+fit
        endfor
        	endif
        			endif
        ;������� ��������� �����������
		PA_null=317.3

		for j=0,Ntarget-1 do begin
		TETA=2*(spectra.PA(j)-PA_null+PA_rot)*!PI/180
		Q=stoks.value(*,1,j)*cos(TETA)-stoks.value(*,2,j)*sin(TETA)
		U=stoks.value(*,1,j)*sin(TETA)+stoks.value(*,2,j)*cos(TETA)
		stoks.value(*,1,j)=Q & stoks.value(*,2,j)=U
		endfor
	    ;����������� ����������� ����������� � ������� �������

	    stoks.value(*,1,0)=stoks.value(*,1,0)-Q_ISM/100.
	    stoks.value(*,2,0)=stoks.value(*,2,0)-U_ISM/100.
        ;���������� ������� �����������
		stoks.value(*,4,*)=sqrt(stoks.value(*,1,*)^2+stoks.value(*,2,*)^2)
		stoks.error(*,4,*)=sqrt(stoks.error(*,1,*)^2+stoks.error(*,2,*)^2)
		stoks.value(*,0,*)=stoks.value(*,0,*)*stoks.value(*,4,*)
		stoks.error(*,0,*)=stoks.value(*,0,*)*stoks.error(*,4,*)
		;���������� ���� ��������� �����������
		stoks.value(*,5,*)=calc_atan(stoks.value(*,1,*),stoks.value(*,2,*))/2
		;�������� ���������������� ����


 ;	RF=WHERE(stoks.value(*,5,0) gt 90, ind) & if ind gt 0 then stoks.value(RF,5,0)=stoks.value(RF,5,0)-90
;	RF=WHERE(stoks.value(*,5,0) lt 25, ind) & if ind gt 0 then stoks.value(RF,5,0)=stoks.value(RF,5,0)+90
	;RF=WHERE(stoks.value(*,5,0) gt 90, ind) & if ind gt 0 then stoks.value(RF,5,0)=stoks.value(RF,5,0)-0
		;stoks.value(*,5,*)=atan(stoks.value(*,2,*)/stoks.value(*,1,*))/2*180/!PI
		for k=0,2 do begin
	;	RF=WHERE(stoks.value(*,5,k) gt 180, ind) & if ind gt 0 then stoks.value(RF,5,k)=stoks.value(RF,5,k)-90
		endfor
		stoks.error(*,5,*)=28.4*stoks.error(*,4,*)*100
		end
ENDCASE
;���������� ������� ��������  � ������� ������
stoks_str=['','Q','U','V','P','!9P!3']
stoks_amp=[1,100,100,100,100,1]
min_lim=fltarr(6)-range_Stoks & min_lim(0)=0 & min_lim(5)=-range_Angle
max_lim=fltarr(6)+range_Stoks & max_lim(0)=1 & max_lim(5)= range_Angle
form_str=['','(F5.2)','(F5.2)','(F5.2)','(F5.2)','(F7.1)']

		Wc=5500  & dW=500 & Rc=where(wave gt Wc-dW AND wave lt Wc+dW)
		for j=0,Ntarget-1 do begin

		for k=1,5 do begin
		stoks.value(*,k,j)=stoks_amp(k)*stoks.value(*,k,j)
		stoks.error(*,k,j)=stoks_amp(k)*stoks.error(*,k,j)
		if total(stoks.value(Rc,k,j)) ne 0 then begin
	   	robomean, stoks.value(Rc,k,j),3,0.5,avg_val,rms_val
		stoks.mean(k,j)=avg_val
		stoks.rms(k,j)=rms_val
		stoks.param(k,j)=STRCOMPRESS(stoks_str(k)+'('+string(Wc,format='(I4)')+')='$
			+string(stoks.mean(k,j),format=form_str(k))+'!9 +!3'$
			+string(stoks.rms(k,j),format=form_str(k)))

		endif
		stoks.yrange(*,k,j)=stoks.mean(k,j)+[min_lim(k),max_lim(k)]
		if k eq 4 then 		stoks.yrange(0,k,j)=0
		endfor
		;����������� �������
		RW=where(wave gt xrng(0) and wave lt xrng(1))

		order=ALOG(max(stoks.value(RW,0,j)))
	    rms=stdev(stoks.value(RW,0,j),mean)
	    stoks.yrange(*,0,j)=[0,mean+6*rms]

		endfor




cont:
end



Pro Stoks_analyz_event, event
COMMON Wbase,base_plot
COMMON SET_widget,BUTTONS
COMMON GRAPH,spectra,stoks,wave,dir,z
COMMON FILES,LOG_dir,LOGFILE,pro_dir,wdir,woll,head,output_file
COMMON SPLOT,stoks_in,stoks_out,xrng,plot_ind,plot_items,target,sent_ind,range_Stoks,range_Angle
COMMON OPTION,wx,thresh,PA,PA_null,PA_rot,Q_ISM,U_ISM,err,corr
; When a widget is selected, put its User Value into 'eventval':

;
;IF (event.select EQ 1) THEN PRINT, event.value, ' selected.' $
;ELSE PRINT, event.value, ' de-selected.'
; Perform actions based on the User Value of the event:
WIDGET_CONTROL, event.id, GET_UVALUE = eventval
;print,event.value
CASE eventval OF
;print,eventval


'set_PATH':   BEGIN
  			log_dir=DIALOG_PICKFILE(PATH='/data6/SCORPIO/sppol_pipeline_v2023.8/',/DIRECTORY)
  			WIDGET_CONTROL, BUTTONS(1), SET_VALUE=log_dir
  			end
'LOG_input':BEGIN
   		    CD,log_dir+'LOGS/',CURRENT=pro_dir
   			filename=pickfile(/read,filter='*.txt')
   			CD,pro_dir
			LOGFILE=filename
   			print,logfile
  			name_out=FILE_BASENAME(filename)
  			WIDGET_CONTROL, BUTTONS(2), SET_VALUE=name_out
			wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'/')
			wdir=log_dir+wdir(N_elements(wdir)-2)+'/'
			;mode analyzer
			ANALYZER=sxpar(read_table(LOGFILE),'ANALYZER')
			WIDGET_CONTROL, BUTTONS(3), SET_VALUE='ANALYZER '+ANALYZER
			;��������� ��������� �������� ����������
			wx=5
			WIDGET_CONTROL, BUTTONS(5), SET_VALUE=string(wx,format='(I2)')
			PA_rot=0
			WIDGET_CONTROL, BUTTONS(6), SET_VALUE=string(PA_rot,format='(I2)')
			Q_ISM=0 & U_ISM=0
			WIDGET_CONTROL, BUTTONS(7), SET_VALUE=string(Q_ISM,format='(F5.2)')
			WIDGET_CONTROL, BUTTONS(8), SET_VALUE=string(U_ISM,format='(F5.2)')
			corr=0
			WIDGET_CONTROL, BUTTONS(9), SET_VALUE=corr
			thresh=3
			WIDGET_CONTROL, BUTTONS(10), SET_VALUE=string(thresh,format='(I2)')
			range_Stoks=2
			WIDGET_CONTROL, BUTTONS(11), SET_VALUE=string(range_Stoks,format='(I2)')
			range_Angle=90
			WIDGET_CONTROL, BUTTONS(12), SET_VALUE=string(range_Angle,format='(I2)')
			;spectral coverage
			if   FILE_TEST(wdir+ 'stoks.fit') eq 1 then begin
			stoks_in=readfits(wdir+'stoks.fit',head)
			Nx=sxpar(head,'NAXIS1')
			lambda=findgen(Nx)*sxpar(head,'CDELT1')+sxpar(head,'CRVAL1')
			WIDGET_CONTROL, BUTTONS(4), SET_VALUE=string(min(lambda),format='(I5)')+$
				'-'+string(max(lambda),format='(I4)')+' AA   '+$
				string(sxpar(head,'CDELT1'),format='(F4.1)')+'  A/px'

			ENDIF  ELSE Result = DIALOG_MESSAGE('No file STOKS.FIT!')
			WOLL=sxpar(read_table(LOGFILE),'ANALYZER')
			WOLL=str_sep(WOLL,'-')
			WOLL='WOLL'+WOLL(1)
			if woll eq 'WOLL1' then WIDGET_CONTROL, BUTTONS(9), SENSITIVE=0 $
				 ELSE  WIDGET_CONTROL, BUTTONS(9), SENSITIVE=1
			OUT=FILE_BASENAME(LOGFILE) & OUT=str_sep(out,'.') & OUT=out(0)

			if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot

			output_file=out
			WIDGET_CONTROL,BUTTONS(13),SET_VALUE=output_file
            ;

   			end
'beg_wave':	begin
        	WIDGET_CONTROL,event.id,GET_VALUE=T
        	xrng(0)=FLOAT(T(0))
        	if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
        	end
'end_wave':	begin
        	WIDGET_CONTROL,event.id,GET_VALUE=T
        	xrng(1)=FLOAT(T(0))
        	if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
        	end
'target':		begin
		    target=event.value
		    if TARGET NE 0 then WIDGET_CONTROL,BUTTONS(9),SENSITIVE=0 $
		     		       ELSE WIDGET_CONTROL,BUTTONS(9),SENSITIVE=1
		    if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
		     END
'box':		begin
        	WIDGET_CONTROL,event.id,GET_VALUE=T
        	wx =FLOAT(T(0))
        	if wx le 1 then wx=5
        	WIDGET_CONTROL,BUTTONS(5),SET_VALUE=string(wx,format='(I2)')
        	if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
       		end
'thresh':	begin
        	WIDGET_CONTROL,event.id,GET_VALUE=T
        	thresh =FLOAT(T(0))
            	if thresh lt 1 then thresh=1
        	WIDGET_CONTROL,BUTTONS(10),SET_VALUE=string(thresh,format='(I2)')
        	if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
       		end
'PA_rot':	begin
        	WIDGET_CONTROL,event.id,GET_VALUE=T
        	PA_rot =FLOAT(T(0))
        	if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
            END
'Q_ISM':	begin
        	WIDGET_CONTROL,event.id,GET_VALUE=T
        	Q_ISM=FLOAT(T(0))

        	if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
       		end
'U_ISM':	begin
        	WIDGET_CONTROL,event.id,GET_VALUE=T

        	U_ISM=FLOAT(T(0))

            if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
       		end
'z':		begin
        	WIDGET_CONTROL,event.id,GET_VALUE=T

        	Z=FLOAT(T(0))

            ;if target eq 0 and corr eq 1  then $
        	;robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
       		end

'err':		begin
			err=event.select
        	if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
       		end
'corr':		begin
        	corr=event.select
        	if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
       		end
'Stoks_lim':	begin
        	WIDGET_CONTROL,event.id,GET_VALUE=T
        	range_Stoks=FLOAT(T(0))
        	if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
       		end
'Angle_lim':	begin
        	WIDGET_CONTROL,event.id,GET_VALUE=T
        	range_Angle=FLOAT(T(0))
        	if target eq 0 and corr eq 1  then $
        	robust_estimate_stoks,/order2 ELSE  robust_estimate_stoks
            multiplot
       		end
'menu':		begin
			print,target
			for k=0,5 do IF event.value eq $
			plot_items(k) then plot_ind(k)= event.select
			print,plot_ind
			multiplot
		   	end
'SAVE':		BEGIN
			; ������ �������

			if target eq 0 then begin
			R=where(wave gt xrng(0) and wave lt xrng(1), Nout)
			V=total(stoks.value(*,3,0))
			openw,1,wdir+output_file+'.res'
			out_title='wave'+'    Flux*P  '+' avg_Q'+' rms_Q'+' avg_U'+' rms_U'
			if V ne 0 then  out_title=out_title+' avg_V'+' rms_V'
			out_title=out_title+' avg_P'+' rms_P'+' avg_FI '+' rms_FI'
			printf,1,out_title
			for k=0,Nout-1 do begin
			;            wave         Flux*P

			out_string=string([wave(R(k)),stoks.value(R(k),0,0),$
			;      avg_Q              err_Q
			stoks.value(R(k),1,0),stoks.error(R(k),1,0),$
			;      avg_U            err_U
			stoks.value(R(k),2,0),stoks.error(R(k),2,0)],format='(I4,E12.3,4F6.2)')
			if V ne 0 then out_string=out_string+$
			string([stoks.value(R(k),3,0),stoks.error(R(k),3,0)],format='(2F6.2)')
			out_string=out_string+string([stoks.value(R(k),4,0),stoks.error(R(k),4,0)],format='(2F6.2)')
			out_string=out_string+string([stoks.value(R(k),5,0),stoks.error(R(k),5,0)],format='(F7.1,F6.1)')
			printf,1,out_string
			endfor
			close,1
			endif
			suffics=['obj','null','star']
			;����� �������
			set_plot,'PS'
			device,file=wdir+output_file+'_'+suffics(target)+'.ps',xsize=19,ysize=27,xoffset=1,yoffset=0.5,/portrait
			multiplot
			device,/close
			SPAWN,'gsview64.exe '+wdir+output_file+'_'+suffics(target)+'.ps',/noshell
			set_plot,'WIN'
			multiplot


			END
'output':	BEGIN
			END
'DONE': 	WIDGET_CONTROL, event.top, /DESTROY
ENDCASE


END



PRO Stoks_analyz, GROUP=GROUP
COMMON Wbase,base_plot
COMMON GRAPH,spectra,stoks,wave,dir,z
COMMON SET_widget,BUTTONS
COMMON FILES,LOG_dir,LOGFILE,pro_dir,wdir,woll,head,output_file
COMMON SPLOT,stoks_in,stoks_out,xrng,plot_ind,plot_items,target,sent_ind,range_Stoks,range_Angle
COMMON OPTION,wx,thresh,PA,PA_null,PA_rot,Q_ISM,U_ISM,err,corr

;initial parameters
;����������� ������� ������
DEVICE, GET_SCREEN_SIZE = s
ratio=1.5 & ysz=s(1)*0.9 & xsz=ysz/ratio
Q_ISM=0 & U_ISM=0
target=0
Log_dir='d:\SCORPIO\red_data.pol\'
;Log_dir='h:\red_data.pol\TINATIn\'
pro_dir='h:\WOLLASTON-2.lib\'
xrng=[4500,8000] & wx=5
!p.multi=[0,1,1]
s_win = !D.WINDOW	; Remember the current window so it can be restored
;sw_gau=0
;xs=0 & ys=0
;x_pre=1 & y_pre=1

;x_min=0 & y_min=0 & x_max=1 & y_max=1
; This example uses keywords to define the size of the draw area:
base_plot = WIDGET_BASE(TITLE = 'PLOT POLARIZATION PROPERTIES',/ROW,xoffset=0,yoffset=0)

WIDGET_CONTROL, /MANAGED, base_plot


;
draw_plot = WIDGET_DRAW(base_plot , XSIZE=700, YSIZE=1100,/FRAME);, $
	;/MOTION_EVENTS,  /BUTTON_EVENTS ,$		;generate LOTS of events
	;UVALUE = 'DRAW_EVENT', 	RETAIN = 2)

base_par=WIDGET_BASE(base_plot,/col,xsize=170,ysize=30)
;set PATH
BUT_PATH=WIDGET_BUTTON(base_par,value='SET PATH',uvalue='set_PATH')
text_PATH=WIDGET_TEXT(BASE_par,value=log_dir,XSIZE=25,/align_center,/frame)
;OPEN LOGFILE
BUT_LOG=WIDGET_BUTTON(BASE_par,VALUE='OPEN LOGFILE',uvalue='LOG_input')
text_LOG=WIDGET_TEXT(base_par, XSIZE=25,/align_center)
;INFO ANALYZER
info_POL=WIDGET_LABEL(base_par,value='   ',xsize=160,/align_center,/frame)
INFO=WIDGET_LABEL(base_par,value='INIT SPECTRAL RANGE:',xsize=170,/align_center)
info_WAVE=WIDGET_LABEL(base_par,value='   ',xsize=160,/align_center,/frame)
INFO=WIDGET_LABEL(base_par,value='OUTPUT SPECTRAL RANGE:',xsize=170,/align_center)
;new spectral coverage
base_wave=WIDGET_BASE(base_par,/row,xsize=170)
INFO=WIDGET_LABEL(base_wave,value='OUTPUT',/align_center)
beg_wave=WIDGET_TEXT(base_wave,VALUE=string(xrng(0),format='(I4)'),uvalue='beg_wave',/edit,/align_center,xsize=4)
INFO=WIDGET_LABEL(base_wave,value='- ',/align_center)
end_wave=WIDGET_TEXT(base_wave,VALUE=string(xrng(1),format='(I4)'),uvalue='end_wave',/edit,/align_center,xsize=4)
INFO=WIDGET_LABEL(base_wave,value='AA',/align_center)

;
INFO=WIDGET_LABEL(base_par,value='OUTPUT TARGET',/align_center,/frame,xsize=160)

target_but = CW_BGROUP(base_par,['object','unpolarized star','polarized star'] , /EXCLUSIVE,uvalue='target', set_value=target)
INFO=WIDGET_LABEL(base_par,value='OPTIONS OF THE OUTPUT PLOT',/align_center,/frame,xsize=160)
;
base_wave=WIDGET_BASE(base_par,/row,xsize=170)
INFO=WIDGET_LABEL(base_wave,value='window of robust estimate, px',/align_center)
box= WIDGET_TEXT(base_wave,VALUE=string(wx,format='(I2)'),uvalue='box',/edit,/align_center,xsize=2)

;
thresh=3
base_wave=WIDGET_BASE(base_par,/row,xsize=170)
INFO=WIDGET_LABEL(base_wave,value='  threshold of robust estimate ',/align_center)
threshold= WIDGET_TEXT(base_wave,VALUE=string(thresh,format='(I2)'),uvalue='thresh',/edit,/align_center,xsize=2)
;
PA_rot=0
base_wave=WIDGET_BASE(base_par,/row,xsize=170)
INFO=WIDGET_LABEL(base_wave,value='rotation plane of polarization',/align_center)
corr_angle= WIDGET_TEXT(base_wave,VALUE=string(PA_rot,format='(I3)'),uvalue='PA_rot',/edit,/align_center,xsize=3)
;

base_wave=WIDGET_BASE(base_par,/row,xsize=170)
INFO=WIDGET_LABEL(base_wave,value='ISM(%):  Q=',/align_center)
QISM=WIDGET_TEXT(base_wave,VALUE=string(Q_ISM,format='(F5.2)'),uvalue='Q_ISM',/edit,/align_center,xsize=5)

INFO=WIDGET_LABEL(base_wave,value='  U=',/align_center)
UISM=WIDGET_TEXT(base_wave,VALUE=string(U_ISM,format='(F5.2)'),uvalue='U_ISM',/edit,/align_center,xsize=5)
base_z=WIDGET_BASE(base_par,/row,xsize=170)
INFO=WIDGET_LABEL(base_z,value='redshift of object',/align_center)
z=0
zlab=WIDGET_TEXT(base_z,VALUE=string(z,format='(F10.7)'),uvalue='z',/edit,/align_center,xsize=10)

err=0
ERR_but = CW_BGROUP(base_par,'plot errors bar', /NONEXCLUSIVE,uvalue='err', set_value=err)
corr=0
CORR_but = CW_BGROUP(base_par,'correction 2nd order', /NONEXCLUSIVE,uvalue='corr', set_value=corr)
base_wave=WIDGET_BASE(base_par,/row,xsize=170)

INFO=WIDGET_LABEL(base_par,value='OUTPUT PLOTS',/align_center,/frame,xsize=160)

range_Stoks=2
base_wave=WIDGET_BASE(base_par,/row,xsize=170)
INFO=WIDGET_LABEL(base_wave,value='Output range of Stokes (%)',/align_center)
Stoks_lim=WIDGET_TEXT(base_wave,VALUE=string(range_Stoks,format='(I2)'),uvalue='Stoks_lim',/edit,/align_center,xsize=5)

range_angle=90
base_wave=WIDGET_BASE(base_par,/row,xsize=170)
INFO=WIDGET_LABEL(base_wave,value='0utput range of angle (deg)',/align_center)
Angle_lim=WIDGET_TEXT(base_wave,VALUE=string(range_angle,format='(I2)'),uvalue='Angle_lim',/edit,/align_center,xsize=5)


;******
plot_ind=[0,1,1,0,1,1]
plot_items = [ 'polarized flux','Q-Stoks parameter','U-Stoks parameter','V-Stoks parameter ','degree of polarization','angle of polarization plane']
for k=0,5 do $
menu_plot = CW_BGROUP(base_par, plot_items(k), /NONEXCLUSIVE,uvalue='menu', IDS = button,set_value=plot_ind(k), /RETURN_NAME)
;values=WIDGET_LABEL(base_inform,value=' ',xsize=465)

button = WIDGET_BUTTON(base_par ,UVALUE = 'SAVE',VALUE = 'SAVE RESULT')
base_out=WIDGET_BASE(base_par,/row,xsize=170)
INFO=WIDGET_LABEL(base_out,value='0utput file',/align_center)
OUT_NAME=WIDGET_TEXT(base_out,VALUE=' ',uvalue='output',/edit,/all,/align_center,xsize=17)

button = WIDGET_BUTTON(base_par , $
		UVALUE = 'DONE', $
		VALUE = 'DONE')

 		;    0          1           2        3        4      5      6       7	  8		  9
BUTTONS=[ draw_plot , text_PATH,text_LOG,info_POL,info_WAVE,box,corr_angle,QISM,UISM,CORR_but,$
        ;    10      11        12       13      14   15     16     17       18     19      20	21		22
		threshold,Stoks_lim,Angle_lim,OUT_NAME];,BUT_DEFAULT,BUT_LOAD,BUT_SAVE,Xbeg,wx,Ntra,DegTra,comp,label_pos,obj_pos,label_w,obj_w,radiusX,radiusY,$
		;    23      24            25                26        27       28
		;shiftPA_W,out_beg_wave_W,out_end_wave_W,width_integr];,P_ISM_W,PA_ISM_W]

WIDGET_CONTROL, base_plot , /REALIZE



WIDGET_CONTROL, draw_plot, GET_VALUE=num_win


DRAWPLOT
; Hand off control of the widget to the XMANAGER:
XMANAGER, "Stoks_analyz", base_plot , GROUP_LEADER=GROUP, /NO_BLOCK

END
