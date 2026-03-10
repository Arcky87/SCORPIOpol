pro geometry_WOLL2,neon,tra,DY=dy,X0,Y0,X1,Y1,SCALE=scale,TRESH=tresh,EPS=eps,PLOT=plot,TITLE=title
;neon - изображение спектра сравнения в режиме длинной щели
; tra - екстраполированные трактории спектра плоского поля в режиме точечной маски
; scale - увеличение числа траекторий по высоте щели
;  XX - координаты узлов модели вдоль дисперсии
;  YY - координаты узлов модели вдоль щели
;************************************
;eps  	окрестность для точного определения вершины линии
;tresh	уровень в rms для выделение линии
;************************************

a=size(neon) & Nx=a(1) & Ny=a(2)
if not(keyword_set(scale)) then scale=2
if not(keyword_set(tresh)) then tresh=2 ;10 - mine IY
if not(keyword_set(eps)) then eps=20 ; 40 - mine IY
if not(keyword_set(dy)) then dy=10
if keyword_set(plot) then P=1 else P=0
;Увеличиваем число траекторий по высоте щели
tra=expand_traectory(tra,SCALE=scale)
if keyword_set(TITLE) then titl=title ELSE titl=''
b=size(tra) & Ntra=b(2)

create_repers_WOLL2,neon,tra,xrep,yrep,TRESH=tresh,EPS=eps,plot=P,TITL=titl

a=size(xrep) & Nrep=a(2)
SY=140 & SX=972

if P eq 1 then window,2,xsize=SX,ysize=SY,ypos=620+(SY+40),xpos=0,title='approximation curvature of lines' +titl,retain=2

!P.multi=[0,1,1]
map=neon
;robomean,congrid(map,Nx/10,Ny/10),3,0.5,mean,rms
; 255-bytscl(ALOG10(congrid(neon,SX,SY)),1,3)
if P eq 1 then  begin
tv,255-bytscl(congrid(map,SX,SY),-10,1000);mean-5*rms,mean+rms*50)
plot,[0,Nx-1],[0,Ny-1],xs=1,yst=1,/nodata,position=[0,0,1,1],/noerase
for k=0,Ntra-1 do oplot,tra(*,k)
oplot,xrep,yrep,psym=6,symsize=0.5,thick=2,color=2e6
endif
;построение модели геометрии
Ndeg=2 ; степень полинома аппроксимаци линий вдоль щели
; убираем реперы крайних линий (самую синюю и самую красную) Убрал 
;xrep=xrep(*,1:Nrep-1) & yrep=yrep(*,1:Nrep-1) & Nrep=Nrep-2 ; строка не нужна для woll2 ???
print,'Nrep=',Nrep
ff=fltarr(Ndeg+1,Nrep)
;аппроксимируем линии полиномом Ndeg
for k=0,Nrep-1 do begin
ff(*,k)=goodpoly(yrep(*,k),xrep(*,k),Ndeg,3,Xfit) ; 3 --- default IY
;print,'line=',k,stdev(xrep(*,k)-Xfit)
endfor
;аппроксимация коеффициентов наклона и кривизны линий
ap=fltarr(3,Ndeg+1)
for k=1,Ndeg do begin
ap(*,k)=goodpoly(ff(0,*),ff(k,*),2,3,fit) ; default 2,3 IY
ff(k,*)=fit
endfor
line=fltarr(Ny,Nrep)
y=findgen(Ny)
;формирование линий
for k=0,Nrep-1 do begin
fit=0 & for i=0,Ndeg do fit=fit+ff(i,k)*y^i
line(*,k)=fit
endfor
;добавление линий на краях
;d_tmp=45;
;tmp=fltarr(Ny,Nrep+2)
;;tmp(*,0)=ap(0,1)*y+ap(0,2)*y^2+d_tmp  ; original, works well
;tmp(*,0) = line(*,0) - (min(line(*,0)) - 5)
;tmp(*,Nrep+1)=line(*,Nrep-1)+(Nx-5-max(line(*,Nrep-1)))
;tmp(*,1:Nrep)=line(*,*)
;line=tmp & Nlin=Nrep+2
;добавление линий на краях путем экстраполяции
left_coeffs = fltarr(Ndeg+1)
right_coeffs = fltarr(Ndeg+1)

for i=0,Ndeg do begin
  coeff_trend = goodpoly(ff(0,*), ff(i,*), 1, 2)
  left_coeffs(i) = coeff_trend[0] + coeff_trend[1] * 5.0
  right_coeffs(i) = coeff_trend[0] + coeff_trend[1] * (Nx-5.0)
endfor
  line_left = 0.0
  line_right = 0.0
  for i=0,Ndeg do begin
    line_left = line_left + left_coeffs(i) * y^i
    line_right = line_right + right_coeffs(i) * y^i
endfor
tmp=fltarr(Ny,Nrep+2)
tmp(*,0) = line_left
tmp(*,Nrep+1) = line_right
tmp(*,1:Nrep) = line(*,*)
line=tmp & Nlin=Nrep+2

;for k=0,Nlin-1 do oplot,line(*,k),y,color=3e6
; Исходные линии 
;for k=1,Nlin-1 do oplot,line(*,k),y,color=3e6 ; 
; Левая граница - красный
oplot,line(*,0),y,linestyle=2,color=3e6
; Правая граница - синий  
oplot,line(*,Nlin-1),y,linestyle=2,color=3e6
;двумерная аппроксимация траекторий
coeff_tra=fltarr(3,Ntra) & x=findgen(Nx)
for i=0,Ntra-1 do begin
coeff_tra(*,i)=goodpoly(x,tra(*,i),2,3,Yfit)
endfor
low=fltarr(3)  & high=fltarr(3)
for i=0,2 do begin
ff=goodpoly(tra(Nx/2,*),coeff_tra(i,*),1,2,Yfit) ; default 1,2
low(i)=ff(0)+ff(1)*(tra(Nx/2,0)-dy)
high(i)=ff(0)+ff(1)*(tra(Nx/2,Ntra-1)+dy)
endfor
tra_low=0 & for i=0,2 do tra_low=tra_low+low(i)*x^i
tra_high=0 & for i=0,2 do tra_high=tra_high+high(i)*x^i
oplot,x,tra_low,linestyle=2,color=3e6
oplot,x,tra_high,linestyle=2,color=3e6
tmp=fltarr(Nx,Ntra+2)
tmp(*,0)=tra_low
tmp(*,1:Ntra)=tra
tmp(*,Ntra+1)=tra_high
tra=tmp & Ntra=Ntra+2
;экстраполяция траекторий и линий за края формата
ext=50  ;50  - 2x1
xx=findgen(Nx+2*ext)-ext  ; увеличиваем область вдоль дисперсии
yy=findgen(Ny+2*ext)-ext ; увеличиваем область вдоль щели
line_ext=fltarr(Ny+2*ext,Nlin) ; подготавливаем увеличенные массивы для линий
tra_ext=fltarr(Nx+2*ext,Ntra) ; то же самое для траекторий
for k=0,Nlin-1 do begin
line_ext(*,k)=INTERPOL( line(*,k),y,yy) ; интерполяция на иррегулярной сетке
;ОТЛАДКА Проверяем наличие (0,0) в интерполированных данных
zero_points = where((line_ext(*,k) eq 0) and (yy eq 0), count)
if count gt 0 then begin
    print, 'Line ', k, ': Found ', count, ' points at (0,0) after interpolation'
    print, '  Source line range: ', min(line(*,k)), 'to', max(line(*,k))
endif
;КОНЕЦ ОТЛАДКИ
endfor
for k=0,Ntra-1 do begin
tra_ext(*,k)=INTERPOL( tra(*,k),x,xx)
;ОТЛАДКА Проверяем наличие (0,0) в интерполированных данных
zero_points = where((tra_ext(*,k) eq 0) and (xx eq 0), count)
if count gt 0 then begin
    print, 'Trajectory ', k, ': Found ', count, ' points at (0,0) after interpolation'
      print, '  Source trajectory range: ', min(tra(*,k)), 'to', max(tra(*,k))
endif
;КОНЕЦ ОТЛАДКИ
endfor
;определение координат точек пересечения линий и траекторий
x1=fltarr(Ntra,Nlin) & y1=x1
zero_count = 0


for i=0,Nlin-1 do begin
for j=0,Ntra-1 do begin
 pos=intersection_WOLL2(line_ext(*,i),tra_ext(*,j),3) ; Ищем пересечение для каждой линии с каждой траекторией
;ОТЛАДКА
 if (pos[0] eq 0 and pos[1] eq 0) then begin
    zero_count = zero_count + 1
    print, 'Zero intersection at line=', i, 'traj=', j
    print, '  Line ', i, ' range: ', min(line_ext(*,i)), 'to', max(line_ext(*,i))
    print, '  Traj ', j, ' range: ', min(tra_ext(*,j)), 'to', max(tra_ext(*,j))
endif
;КОНЕЦ ОТЛАДКИ
 R=where(pos lt 0,ind) & if ind ne 0 then print,'negative value',i,j, pos(0),pos(1) ; обработка ошибки отрицательного значения
 x1(j,i)=pos(0) & y1(j,i)=pos(1)
endfor
endfor
print, 'Total zero intersections: ', zero_count
;образование исходной сетки
X0=X1 & for k=0,Ntra-1 do  X0(k,*)=X1(Ntra/2,*)
Y0=Y1 & for k=0,Nlin-1 do Y0(*,k)=Y1(*,Nlin/2)
X0=reform(X0,Ntra*Nlin) & X1=reform(X1,Ntra*Nlin)
Y0=reform(Y0,Ntra*Nlin) & Y1=reform(Y1,Ntra*Nlin)

; ФИНАЛЬНАЯ ПРОВЕРКА НА (0,0)
zero_points = where((X0 eq 0) and (Y0 eq 0) and (X1 eq 0) and (Y1 eq 0), final_zero_count)
if final_zero_count gt 0 then begin
  print, '=== FINAL RESULT: ', final_zero_count, ' POINTS AT (0,0) ==='
  print, 'These points will be saved in geometry.fit'
endif

if P eq 1 then begin
window,4,xsize=SX+800,ysize=SY+500,xpos=0,title='input (cross) and output (square) grid'+titl,retain=2  ;ypos=620+(SY+40)*2,
plot,[min(xx),max(xx)],[min(yy),max(yy)],xst=1,yst=1,/nodata,position=[0,0,1,1],/norm
for k=0,Ntra-1 do oplot,xx,tra_ext(*,k),linestyle=1
for k=0,Nlin-1 do oplot,line_ext(*,k),yy,linestyle=1
oplot,x0,y0,psym=1,symsize=2,color=1e5
oplot,x1,y1,psym=6,symsize=0.5
endif
Nc=N_elements(X0)
end
