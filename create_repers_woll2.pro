
pro create_repers_WOLL2,neon,tra,XPOS=pos,YPOS=ypos,xrep,yrep,TRESH=tresh,EPS=eps,PLOT=plot,WIN=win,TITL=titl
;формирование реперов по спектру сравнения
R=where(neon lt 1, ind) & if ind gt 1 then neon(R)=1
if not(keyword_set(win)) then win=0
if not(keyword_set(tresh)) then tresh=2
if not(keyword_set(eps)) then eps=12
a=size(neon) & Nx=a(1) & Ny=a(2)
b=size(tra) & Ntra=b(2)
vector=fltarr(Nx)
x=findgen(Nx)
Npos=200
xpos=fltarr(Npos,Ntra) & ypos=fltarr(Npos,Ntra) & Npk=fltarr(Ntra)
;*********************
wy=3  & wx=2
;*********************
SX=1972 & SY=640
if keyword_set(plot) then begin
	window,win,xsize=SX,ysize=SY,xpos=0,ypos=620+(SY+40)*win,title='create array of repers'+TITL,retain=2
	tv,255-bytscl(ALOG10(congrid(neon,SX,SY)),1,3)
	plot,[0,Nx-1],[0,Ny-1],xst=1,yst=1,/noerase,/nodata,$
		position=[0,0,1,1],/norm
		for k=0,Ntra-1 do oplot,tra(*,k),color=1e7
endif

;поиск и выделение линий
for k=0,Ntra-1 do begin
	for kx=0,Nx-1 do vector(kx)=total(neon(kx,tra(kx,k)-wy:tra(kx,k)+wy),2)
	;vector=ALOG10(vector)
	xpk=find_peaks(vector,tresh=1)
	RR=where(xpk gt wx*2 and xpk lt Nx-1 -2*wx) & xpk=xpk(RR)
	Npk(k)=N_elements(xpk)
		;точное определение положения пика
		for j=0,Npk(k)-1 do begin
			P=goodpoly(x(xpk(j)-wx:xpk(j)+wx),vector(xpk(j)-wx:xpk(j)+wx),2,2,fit) ; 2,2 ---default IY
			xpk(j)=-P(1)/P(2)/2
		endfor
	if keyword_set(plot) then oplot,xpk,tra(xpk,k),psym=6,color=2e6,symsize=0.75
	xpos(0:Npk(k)-1,k)=xpk
	ypos(0:Npk(k)-1,k)=tra(xpk,k)
endfor
;отождествление линий в eps-окрестности
xrep=0 & yrep=0
count=intarr(Ntra)
for j=0,Npk(0)-1 do begin
	count(0)=j
		for k=1,Ntra-1 do begin
		R=WHERE(ABS(xpos(*,k)-xpos(j,0)) LT eps,ind)
			if ind gt 0 then count(k)=R(ind-1) else count(k)=R(0)
		endfor
	R=where(count gt -1,ind)
	if ind eq Ntra then begin
		for k=0,Ntra-1 do xrep=[xrep,xpos(count(k) ,k)]
		for k=0,Ntra-1 do yrep=[yrep,ypos(count(k) ,k)]
	endif
endfor
xrep=xrep(1:N_elements(xrep)-1)
yrep=yrep(1:N_elements(yrep)-1)
Nline=N_elements(xrep)/Ntra
xrep=reform(xrep,Ntra,Nline)
yrep=reform(yrep,Ntra,Nline)
index=intarr(Nline)
for k=0,Nline-1 do begin
	f=goodpoly(yrep(*,k),xrep(*,k),1,3,Xfit) ; default 1,3 IY
	err=stdev(xrep(*,k)-Xfit)
	if err lt 0.5 then begin
		if keyword_set(plot) then oplot,xrep(*,k),yrep(*,k),color=3,psym=6
		index(k)=1
		endif
endfor

R=where(index eq 1)
xrep=xrep(*,R)
yrep=yrep(*,R)
end
