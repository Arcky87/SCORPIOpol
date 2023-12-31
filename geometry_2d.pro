pro geometry_2D,neon,tra,DY=dy,X0,Y0,X1,Y1,SCALE=scale,TRESH=tresh,EPS=eps,PLOT=plot,TITLE=title
;neon - ����������� ������� ��������� � ������ ������� ����
; tra - ������������������ ��������� ������� �������� ���� � ������ �������� �����
; scale - ���������� ����� ���������� �� ������ ����
;  XX - ���������� ����� ������ ����� ���������
;  YY - ���������� ����� ������ ����� ����
;************************************
;eps  	����������� ��� ������� ����������� ������� �����
;tresh	������� � rms ��� ��������� �����
;************************************

a=size(neon) & Nx=a(1) & Ny=a(2)
if not(keyword_set(scale)) then scale=2
if not(keyword_set(tresh)) then tresh=2
if not(keyword_set(eps)) then eps=12
if not(keyword_set(dy)) then dy=10
if keyword_set(plot) then P=1 else P=0
tra=expand_traectory(tra,SCALE=scale)
if keyword_set(TITLE) then titl=title ELSE titl=''
b=size(tra) & Ntra=b(2)
create_repers,neon,tra,xrep,yrep,TRESH=tresh,EPS=eps,plot=P,TITL=titl

a=size(xrep) & Nrep=a(2)
N_y=140
if P eq 1 then window,2,xsize=Nx/2,ysize=N_y,ypos=N_y+40,xpos=0,title='approximation curvature of lines' +titl
!P.multi=[0,1,1]
map=ALOG10(neon+1000);
map=neon
;robomean,congrid(map,Nx/10,Ny/10),3,0.5,mean,rms
if P eq 1 then  begin
tv,255-bytscl(congrid(map,Nx/2,N_y),-10,1000);mean-5*rms,mean+rms*50)
plot,[0,Nx-1],[0,Ny-1],xs=1,yst=1,/nodata,position=[0,0,1,1],/noerase
for k=0,Ntra-1 do oplot,tra(*,k)
oplot,xrep,yrep,psym=6,symsize=0.5
endif
;���������� ������ ���������
Ndeg=2
ff=fltarr(Ndeg+1,Nrep)
print,'Nrep=',Nrep
for k=0,Nrep-1 do begin
;print,k,yrep(*,k)
;print,k,xrep(*,k)
ff(*,k)=goodpoly(yrep(*,k),xrep(*,k),Ndeg,3,Xfit)
;oplot,Xfit,yrep(*,k)
print,'line=',k,stdev(xrep(*,k)-Xfit)
endfor
;������������� ������������� ������� � �������� �����
ap=fltarr(3,Ndeg+1)
for k=1,Ndeg do begin
ap(*,k)=goodpoly(ff(0,*),ff(k,*),2,3,fit)
ff(k,*)=fit
endfor
line=fltarr(Ny,Nrep)
y=findgen(Ny)
;������������ �����
for k=0,Nrep-1 do begin
fit=0 & for i=0,Ndeg do fit=fit+ff(i,k)*y^i
line(*,k)=fit
endfor
;���������� ����� �� �����
tmp=fltarr(Ny,Nrep+2)
tmp(*,0)=ap(0,1)*y+ap(0,2)*y^2+20
tmp(*,Nrep+1)=line(*,Nrep-1)+(Nx-3-max(line(*,Nrep-1)))
tmp(*,1:Nrep)=line(*,*)
line=tmp & Nlin=Nrep+2
for k=0,Nlin-1 do oplot,line(*,k),y,color=3e6
;��������� ������������� ����������

coeff_tra=fltarr(3,Ntra)  & x=findgen(Nx)
for i=0,Ntra-1 do begin
coeff_tra(*,i)=goodpoly(x,tra(*,i),2,3,Yfit)
endfor
low=fltarr(3)  & high=fltarr(3)
for i=0,2 do begin
ff=goodpoly(tra(Nx/2,*),coeff_tra(i,*),1,2,Yfit)
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
;������������� ���������� � ����� �� ���� �������
ext=80
xx=findgen(Nx+2*ext)-ext
yy=findgen(Ny+2*ext)-ext
line_ext=fltarr(Ny+2*ext,Nlin)
tra_ext=fltarr(Nx+2*ext,Ntra)
for k=0,Nlin-1 do line_ext(*,k)=INTERPOL( line(*,k),y,yy)
for k=0,Ntra-1 do tra_ext(*,k)=INTERPOL( tra(*,k),x,xx)

;����������� ��������� ����� ����������� ����� � ����������
x1=fltarr(Ntra,Nlin) & y1=x1
for i=0,Nlin-1 do begin
for j=0,Ntra-1 do begin
pos=intersection(line_ext(*,i),tra_ext(*,j),5)
R=where(pos lt 0,ind) & if ind gt 1 then print,'negative value',i,j
x1(j,i)=pos(0) & y1(j,i)=pos(1)
endfor
endfor
;����������� �������� �����
X0=X1 & for k=0,Ntra-1 do  X0(k,*)=X1(Ntra/2,*)
Y0=Y1 & for k=0,Nlin-1 do Y0(*,k)=Y1(*,Nlin/2)
X0=reform(X0,Ntra*Nlin) & X1=reform(X1,Ntra*Nlin)
Y0=reform(Y0,Ntra*Nlin) & Y1=reform(Y1,Ntra*Nlin)
if P eq 1 then begin
window,4,xsize=Nx/2,ysize=N_y,ypos=(N_y+40)*2,xpos=0,title='input (cross) and output (square) grid'+titl
plot,[min(xx),max(xx)],[min(yy),max(yy)],xst=1,yst=1,/nodata,position=[0,0,1,1],/norm
for k=0,Ntra-1 do oplot,xx,tra_ext(*,k),linestyle=1
for k=0,Nlin-1 do oplot,line_ext(*,k),yy,linestyle=1
oplot,x0,y0,psym=1,symsize=2,color=1e5
oplot,x1,y1,psym=6,symsize=0.5
endif
Nc=N_elements(X0)
end
