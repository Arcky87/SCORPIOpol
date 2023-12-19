
function traectory_LS,ima,NP=Np,X_BEG=x_beg,X_END=X_END,WX=wx,NDEG=Ndeg,PLOT=plot
;построение траектории для изображения FLAT с маской SCORPIO-2
;Np		- число траекторий
;X_beg	- координата первого строба вдоль щели
;WX		- ширина строба
;Ndeg	- степень полинома аппроксимации траектории
!P.multi=[0,1,1]
a=size(ima) & Nx=a(1) & Ny=a(2)
if not(keyword_set(Np)) then Np=10
if not(keyword_set(X_beg)) then X_beg=wx
if not(keyword_set(X_end)) then X_end=Nx-1
if not(keyword_set(WX)) then wx=80
if not(keyword_set(Ndeg)) then Ndeg=3
Npos=fix((x_end-x_beg )/wx/2)
xpos=findgen(Npos)*2*wx+x_beg
;print,xpos
;print,Npos,Ny
;print,xpos
y=findgen(Ny)
vector_y=fltarr(Npos,Ny)
ypk=fltarr(Np,Npos)
;формирование разрезов (стробов)

for k=0,Npos-1 do vector_y(k,*)=total(ima(xpos(k)-wx:xpos(k)+wx,*),1)
;pk=find_peaks(vector_y(0,*),W=10,/plot)
;return,0
;определение координат пиков среднего разреза
mpk=find_peaks(vector_y(Npos/2,*),W=5,tresh=50);,/plot)
if N_elements(mpk) ne Np then begin
print,'number peaks not equal '+string(Np,format='(I2)')+'!!!'
RETURN,0
endif
;ypk(*,Npos/2)=pk
print,' mean',mpk
 ;определение позиций всех пиков
for k=0,Npos-1 do begin
pk=find_peaks(vector_y(k,*),W=5,tresh=60,/plot)

print,k,N_elements(pk)
;print,pk
index=compare_list(mpk,pk,50)
ypk(index(*,1),k)=pk(index(*,1))
endfor
;построение траекторий  flat
x=findgen(Nx)
tra=fltarr(Nx,Np)
wdelete,2
if keyword_set(plot) then begin
window,1,xsize=1000,ysize=500,xpos=0,ypos=0,title='SPECTRA TRAECTORY'
;set_plot,'PS'

;device,file='g:\traectory.eps',/portrait,/encapsulated,xsize=16,ysize=10
 plot,[0,Nx],[0,Ny],/nodata,xst=1,yst=1,$
;plot,[0,Nx],[50,450],/nodata,xst=1,yst=1,yminor=5,title='Image binning 2x2 px',$
ytitle='Distance along slit, px',xtitle='Distance along dispersion, px',charsize=1
endif
for j=0,Npos-1 do begin
R=where(ypk(*,j) ne 0, ind)
if ind ne 0 and ind ne Np then ypk(*,j)=0
endfor
for k=0,Np-1 do begin
;remove zero values
R=where(ypk(k,*) ne 0,ind)
;print,R
f=goodpoly(xpos(R),ypk(k,R),Ndeg,3,Yfit)

robomean,ypk(k,R)-Yfit,3,0.5,avg_fit,rms_fit

for j=0,Ndeg do tra(*,k)=tra(*,k)+f(j)*x^j
if keyword_set(plot) then begin
;oplot,xpos(R),ypk(k,R),psym=6,symsize=0.6
oplot,xpos(*),ypk(k,*),psym=6,symsize=0.6
oplot,x,tra(*,k)
xyouts,Nx/2,tra(Nx/2,k)+10,'rms ='+string(rms_fit,format='(F5.2)')+' px',/data,align=0.5
endif

endfor
;device,/close
;set_plot,'WIN'
return,tra
end

;LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path='h:\spectraPOL.log')
;wdir=sxpar(read_table(LOGFILE),'w_dir')


;path='h:\red_data.pol\standard.new\'
;LOGFILE=DIALOG_PICKFILE(/READ, FILTER = '*.txt',path=path+'LOGS')

;wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
;wdir=path+wdir(N_elements(wdir)-2)+'\'

;flat=readfits(wdir+'eta_i.fts')

;R=where(flat lt 0) & flat(R)=0
;tra=create_traectory(flat,NP=10,WX=20,NDEG=2,X_beg=200,/plot)
;writefits,wdir+'traectory.fit',tra

;end