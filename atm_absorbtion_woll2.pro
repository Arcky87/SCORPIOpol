pro atm_absorbtion_WOLL2,LOGFILE,LINE=line,REGION=region,STAR=star,PLOT=plot
;line   	длина волны абсорбции
;region		ширина интервала проведения сонтинуума
tmp=str_sep(LOGFILE,'/')
path=tmp(0)
for j=1,N_elements(tmp)-3 do begin
path=path+'/'+tmp(j)
endfor
path=path+'/'
dir=str_sep(sxpar(read_table(LOGFILE),'w_dir'),'/')
dir=path+dir(N_elements(dir)-2)+'/'
;
j=star
;

Nline=N_elements(line)
cube=readfits(dir+'obj-sky.fts',h)

Nx=sxpar(h,'NAXIS1')  & Ny=sxpar(h,'NAXIS2') & Npol=sxpar(h,'NAXIS3')  & Ncube= sxpar(h,'NAXIS5')
Nexp=fltarr(Ncube)
for k=0,Ncube-1 do Nexp(k)=sxpar(h,'NUMEXP'+string(k+1,format='(I1)'))
lambda_0= sxpar(h,'CRVAL1')  & d_lambda= sxpar(h,'CDELT1')
lambda=findgen(Nx)*d_lambda+lambda_0
pos_line=(line-lambda_0)/d_lambda
pos_reg=region/d_lambda/2
;формирование спектров
p=1
if keyword_set(plot) then p=1
wy=10
spectra=fltarr(Nx,Npol,Nexp(0),Ncube)
spectra=total(cube(*,Ny/2-wy:Ny/2+wy,*,*,*),2)
if p eq 1 then begin
window,0,title='absorbtion of atmosphere  in spectra star #'+string(star,format='(I1)')
!P.multi=[0,1,N_elements(line)]
endif ELSE !P.multi=[0,1,1]

atm_abs=fltarr(Nx,Npol)+1
spectra_atm=fltarr(2*pos_reg+1,Nline)
spectra_reg=fltarr(2*pos_reg+1,Nexp(j))

for p=0,Npol-1 do begin
for k=0,Nline-1 do begin

if p eq 1 then plot,[0,2*pos_reg],[0,1.2],xst=1,yst=1,/nodata,charsize=2,title= FIX(line(k))
for i=0,Nexp(j)-1 do begin
spectra_reg(*,i)=spectra(pos_line(k)-pos_reg:pos_line(k)+pos_reg,p,i,j)
xpos=findgen(2*pos_reg+1)
cont=LOWESS(xpos,spectra_reg(*,i),2*pos_reg,2,1)
;определение области проведения континуума
robomean,cont-spectra_reg(*,i),1,0.5,avg_S,rms_S
R=where(cont-spectra_reg(*,i) lt rms_S)
f=goodpoly(xpos(R),spectra_reg(R,i),1,2)
cont=f(0)+f(1)*xpos
spectra_reg(*,i)=spectra_reg(*,i)/cont
if p eq 1 then oplot,xpos,spectra_reg(*,i)
endfor
;spectra_atm(*,k)=total(spectra_reg,2)/Nexp(j)
for i=0,2*pos_reg do spectra_atm(i,k)=median(spectra_reg(i,*))
R=where(spectra_atm(*,k) gt 1.01,ind) & if ind gt 1 then spectra_atm(R,k)=1;.02
;spectra_atm(*,k)=spectra_atm(*,k)-0.02
if p eq 1 then oplot,xpos,spectra_atm(*,k)-0.2
if p eq 1 then oplot,[0,2*pos_reg],[1,1]*0.8,linestyle=2
endfor
;формирование выходного спектра
for k=0,Nline-1 do atm_abs(pos_line(k)-pos_reg:pos_line(k)+pos_reg,p)= spectra_atm(*,k)
endfor
writefits,dir+'atm_abs.fit',atm_abs,h
end

;LOGFILE=DIALOG_PICKFILE(/read,path='h:\red_data.pol\AGN\LOGS\',FILTER='*.txt')
LOGFILE=DIALOG_PICKFILE(/read,path='/home/elias/SCORPIO/sppol_pipeline_v2023.8/LOGS/',FILTER='*.txt')
;LOGFILE='h:\red_data.pol\LOGS\3C390_140324.txt\
atm_absorbtion_WOLL2,LOGFILE,LINE=[6890,7200,7620],region=300,STAR=1
;atm_absorbtion_WOLL2,LOGFILE,LINE=[6260,6260,6890],region=200,STAR=1
END

