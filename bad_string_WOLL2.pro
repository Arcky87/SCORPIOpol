;correction _bad string WOLL-2
function bad_string_WOLL2,ima,bad,PLOT=plot
a=size(ima) & Nx=a(1) & Ny=a(2)
x=findgen(Nx) & y=findgen(Ny) & vector=y
print,total(ima)
for kx=0,Nx-1 do begin
vector=ima(kx,*)
if keyword_set(plot) and kx eq Nx/2 then begin
Window,0
plot,vector,xst=1,xrange=[-1,1]*bad(1)*5+bad(0),psym=10
endif
vector(bad(0)-bad(1):bad(0)+bad(1))=(vector(bad(0)-3*bad(1):bad(0)-bad(1))+vector(bad(0)+bad(1):bad(0)+3*bad(1)))/2
if keyword_set(plot) and kx eq Nx/2 then oplot,vector,psym=10,thick=2
ima(kx,*)=vector
endfor
print,total(ima)
return,ima
end


dir='e:\sbs1419+538_190216\'
suff='obj'
;suff='obj'
ima=readfits(dir+suff+'_i.fts',h)
Nx=sxpar(h,'NAXIS1') & Ny=sxpar(h,'Naxis2')
Npol=sxpar(h,'NAXIS3')
Nexp=sxpar(h,'NAXIS3')
Ncub=sxpar(h,'NAXIS5')
	;for k=0,Ncub-1 do begin
;for j=0,Nexp-1 do begin
;map(*,*)=ima(*,*,0);,j,k)
;ima(*,*,0,j,k)=bad_string_WOLL2(map,[165,5],/plot)
;ima(*,*,0)=bad_string_WOLL2(map,[165,5],/plot)
;endfor
	;ENDFOR

ima_new=bad_string_WOLL2(ima,[45,5],/plot)

writefits,dir+suff+'_avg.fts',ima_new,h

;;ima=readfits(dir+'neon.fts',h)
;a=size(ima) & Ny=a(2)
;dy=3
;window,2,xsize=800,ysize=400
;j=3
;tv,255-bytscl(congrid(ima(*,*,j),800,400),0,100)
;	for j=0,a(3)-1 do begin
;for k=0,Nx-1 do begin
;ima(k,*,j)=LOWESS(findgen(a(2)),ima(k,*,j),a(2)/8,2,2)
;endfor
;	endfor
;window,3,xsize=800,ysize=400
;j=3
;tv,255-bytscl(congrid(ima(*,*,j),800,400),0,100)
;
;;writefits,dir+suff+'_avg.fts',ima,h

end

