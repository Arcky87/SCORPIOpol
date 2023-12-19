;VPHG tests6216.7

;create etalon neon spectra
;!P.background=2^24-1 & !P.color=0
dir='h:\red_data.pol\Sy1\result\pg1700.tmp\'
file_bmp=DIALOG_PICKFILE(/read,filter='*.bmp',PATH=dir,TITLE='input BMP-image')
map_bmp=READ_BMP(file_bmp)
S=size(map_bmp)
print,S

start:
Window,2,xsize=S(2) ,ysize=S(3)
tv,255-map_bmp(1,*,*)
plot,[0,1],[0,1],xst=1,yst=1,position=[0,0,1,1],/nodata,/noerase,/norm

;manual identefication spectra
sw=1
x_out=0 & y_out=0 &W_out=0
xpre=0 & ypre=0
be:
x=0 & y=0
if sw eq 4 then goto,final
cursor,x,y,1,/data
a=!MOUSE
sw=LONG(a.button)
if xpre ne x and ypre ne y then begin
if sw eq 1 then begin
NN=N_elements(x_out)
if NN lt 3 then oplot,[x],[y],psym=6,color=2.2e8,symsize=2
if NN gt 2 then oplot,[x],[y],psym=7,color=2.2e8,symsize=2
;if NN eq 1 then begin
;input=''
;read,input,prompt='input left down corner coordinates (X,Y)'
;input=str_sep(input,',')
;
;XL_in=FLOAT(input(0)) & XL_out=x
;
;YL_in=FLOAT(input(1)) & YL_out=y
;xyouts,XL_out+0.01,YL_out+0.01,'('+input(0)+','+input(1)+')',$
;	charsize=3,charthick=3,color=1e5,align=0
;endif
;if NN eq 2 then begin
;input=''
;read,input,prompt='input rigth up corner coordinates (X,Y)'
;input=str_sep(input,',')
;XR_in=FLOAT(input(0)) & XR_out=x
;YR_in=FLOAT(input(1)) & YR_out=y
;xyouts,XR_out-0.01,YR_out-0.05,'('+input(0)+','+input(1)+')',$
;	charsize=3,charthick=3,color=1e5,align=1
;
;endif


x_out=[x_out,x]
y_out=[y_out,y]

endif
xpre=x & ypre=y
endif
if sw eq 1 then goto,be

final:
xyouts,0.5,0.5,'GAME OVER!',/norm,align=0.5,charsize=10,charthick=10,color=1e5
;

;create output list lines
N=N_elements(x_out)

x_out=x_out(3:N-1)
y_out=y_out(3:N-1)
N=N-3
for j=0,N-1 do print,x_out(j),y_out(j)
;рисование оцифрованного графика
print,XL_in,XR_in, YL_out,YR_in
x_new=XL_in+(XR_in-XL_in)/(XR_out-XL_out)*(x_out-XL_out)
y_new=YL_in+(YR_in-YL_in)/(YR_out-YL_out)*(y_out-YL_out)
window,3
plot,x_new,y_new,psym=10,xrange=[XL_in,XR_in],xst=1
input=''
read,input,prompt='REPEAT? (y/n)'
if input eq 'y' then goto,start
print,FILE_BASENAME(file_bmp)
file_txt=str_sep(FILE_BASENAME(file_bmp),'.')
file_txt=FILE_DIRNAME(file_bmp)+'\'+file_txt(0)+'.txt'
print,file_txt
close,1
openw,1,file_txt
for j=0,N-1 do printf,1,x_new(j),y_new(j),format='(2F10.4)'
close,1
END
