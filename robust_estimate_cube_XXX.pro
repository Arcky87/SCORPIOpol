function ROBUST_ESTIMATE_cube,cube,W=w,TRESH=tresh
a=size(cube)
Nx=a(1)
if not(keyword_set(tresh)) then tresh=2
if not(keyword_set(w)) then w=10
if w eq 1 then begin
res=vector
goto,fin
endif
;построение узлов
Npos=Nx/w
xpos=findgen(Npos)*w+w/2
res=fltarr(Npos,3)
;робастная оценка среднего и rms в окне W рх

for k=0,Npos-1 do begin
;print,xpos(k)+w/2,Nx
robomean,cube(xpos(k)-w/2:xpos(k)+w/2-1,*),tresh,0.5,mean,rms
res(k,0)=mean
res(k,1)=rms/5
res(k,2)=xpos(k)
endfor
;подавление выбросов в векторе

;fit=LOWESS(xpos,res(*,0),5,2)
;robomean,res(*,0)-fit,3,0.5,mean,rms
;R=where(abs(res(*,0)-fit) gt rms*3,ind) & if ind gt 0 then res(R,0)=fit(R)
fin:
return,res
end
