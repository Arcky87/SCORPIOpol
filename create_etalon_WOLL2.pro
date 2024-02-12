
;VPHG tests6216.7

pro create_etalon_WOLL2,dir
;create etalon neon spectra
;!P.background=2^24-1 & !P.color=0

;dir=def_wdir(LOGFILE)
;goto,cont
;формирование еталонного спектра неона, гамма-коррекция и вычитание фона
cube=readfits(dir+'neon.fts',h)
a=size(cube) & Nx=a(1)  & Ns=a(3)
spectra=fltarr(Nx,Ns)
x=findgen(Nx)
G=0.25
for j=0,Ns-1 do begin
spectra(*,j)=total(cube(*,*,j),2)
R=where(spectra(*,j) lt 1,ind1) & if ind1 gt 1 then spectra(R,j)=1
spectra(*,j)=spectra(*,j)^G
fi_peak,x ,spectra(*,j),0,ipix,xpk,ypk,bkpk,ipk
fon=INTERPOL(bkpk,xpk,x)
spectra(*,j)=spectra(*,j)-fon & R=where(spectra(*,j) lt 0) & spectra(R,j)=0
endfor



window,2,xsize=1600,ysize=500,title=dir
!P.multi=[0,1,1]
plot,spectra(*,0),xst=1,yrange=[-max(spectra(*,0))/20,max(spectra(*,0))*1.2],yst=1;,position=[0,0,1,1];,/norm
oplot,[0,Nx-1],[0,0],linestyle=2

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
if sw eq 1 then begin oplot,[x],[y],psym=1
input=''
read,input,prompt='wavelength?'
w=float(input)
xyouts,x+5,y,string(w,format='(F9.2)'),charsize=0.75,orientation=90,/data
w_out=[w_out,w]
x_out=[x_out,x]
y_out=[y_out,y]

endif
xpre=x & ypre=y
endif
if sw eq 1 then goto,be

final:
;
;create output list lines
N=N_elements(x_out)

x_out=x_out(1:N-1)
w_out=w_out(1:N-1)
N=N-1
x=findgen(Nx)
Ndeg=3
f=goodpoly(x_out,w_out,Ndeg,3)
wav=0 & for j=0,Ndeg do wav=wav+f(j)*x^j
wav_min=10*fix(wav(0)/10) & wav_max=10*fix(wav(Nx-1)/10)   ;comment - rfix?!
disp=fix((wav(Nx/2+1)-wav(Nx/2))/0.1)*0.1 & order=1
print,disp,wav_min,wav_max

xpos=fltarr(5,N) & xpos(0,*)=w_out  & xpos(1,*)=x_out


print,xpos
eps=10
for k=0,N-1 do begin
for j=0,Ns-1 do begin
max_value=max(spectra(xpos(1,k)-eps:xpos(1,k)+eps,j),Nmax)
xpos(1+j,k)=xpos(1,k)-eps+Nmax
endfor
print,xpos(*,k)
endfor
window,3,xsize=1000,ysize=600,title=dir
!P.multi=[0,1,4]
for j=0,Ns-1 do begin
plot,spectra(*,j),xst=1
for k=0,N-1 do begin
oplot,[1,1]*xpos(1+j,k),[0,1e6],linestyle=1
endfor
endfor
;запись эталонных спектров
writefits,dir+'etalon.fit',spectra,h
;запись списка отождествленных линий
openw,1,dir+'etalon.txt'
printf,1,disp,wav_min,wav_max,0,0,format='(F7.2,4I6)'
for j=0,N-1 do printf,1,xpos(*,j),format='(F7.2,4I6)'
close,1
end
log_dir='h:\red_data.pol\'
LOGFILE=DIALOG_PICKFILE(/read,path=log_dir+'LOGS\',FILTER='*.txt')
;LOGFILE=DIALOG_PICKFILE(/read,path='h:\red_data.pol\LOGS\',FILTER='*.txt')
wdir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'\')
    wdir=log_dir+wdir(N_elements(wdir)-2)+'\'
    print,wdir
create_etalon_WOLL2,wdir
end