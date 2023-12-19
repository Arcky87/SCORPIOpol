function ROBUST_ESTIMATE,vector,W=w,TRESH=tresh
Nx=N_elements(vector)
if not(keyword_set(tresh)) then tresh=2
if not(keyword_set(w)) then w=10
if w eq 1 then begin
res=vector
goto,fin
endif
;���������� �����
Npos=Nx/w
xpos=findgen(Npos)*w+w/2
res=fltarr(Npos,3)
;��������� ������ �������� � rms � ���� W ��

for k=0,Npos-1 do begin
;print,xpos(k)+w/2,Nx
robomean,vector(xpos(k)-w/2:xpos(k)+w/2-1),tresh,0.5,mean,rms
res(k,0)=mean
res(k,1)=rms
res(k,2)=xpos(k)
endfor
;���������� �������� � �������

;fit=LOWESS(xpos,res(*,0),5,2)
;robomean,res(*,0)-fit,3,0.5,mean,rms
;R=where(abs(res(*,0)-fit) gt rms*3,ind) & if ind gt 0 then res(R,0)=fit(R)
fin:
return,res
end
