function frames_WOLL2,ima,YC=yc,DY=dy,H=h,BIAS=bias

a=size(ima)
if not(keyword_set(dY)) then dY=0
if not(keyword_set(Yc)) then Yc=[82,204,330,453]
if not(keyword_set(H)) then H=140
if keyword_set(bias) then ima=ima-def_bias(ima,/plot)

if dy ne 0 then begin
yc=yc+dy
ima=shift(ima,0,dy)
endif


Nf=4
out=fltarr(a(1),h,Nf)

for j=0,Nf-1 do out(*,*,j)=ima(*,yc(j)-h/2:yc(j)+h/2-1)

return,out
end

dir='h:\red_data.pol\Arp102b_140325\'
;dir='h:\red_data.pol\3C390_140224\'
dir='h:\red_data.pol\3C390_131103\'
eta=readfits(dir+'eta_i.fts')
H_frame=140
y_c=center_frames_WOLL2(eta(*,*,0))
if y_c(0) lt H_frame/2 then d_y=H_frame/2-Y_c(0)+5 ELSE d_y=0
print,y_c
print,d_y

tmp=readfits(dir+'flat_i.fts')
out=frames_WOLL2(tmp(*,*,0),YC=Y_c,DY=d_Y,H=H_frame,/bias)
a=size(out)
Window,2,xsize=a(1)/2,ysize=a(2)*4
wait, 5

for k=0,3 do begin
map=fltarr(a(1),a(2))
map(*,*)=out(*,*,k)
tv,255-bytscl(congrid(map,a(1)/2,a(2)),0,1000),0,a(2)*k
endfor
end
