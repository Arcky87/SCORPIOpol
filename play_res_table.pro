;play res table
dir='h:\red_data.pol\Sy1\Mkn335_131109'
h=headfits(dir+'\avg_spectra.fit')
z=sxpar(h,'z')  & Ha=6563
name_tab=FILE_BASENAME(dir)
tab=read_table(dir+'\'+name_tab+'.res')
wave=FLOAT(strmid(tab,0,5))
F=FLOAT(strmid(tab,6,10))
P=FLOAT(strmid(tab,41,5))
PA=FLOAT(strmid(tab,53,7))
;
RA=where(ABS(wave-Ha*(1+z)) lt 500, ind)
wave=wave(RA)  & Na=N_elements(wave)
flux=F(RA)/P(RA)/100
P=P(RA)
PA=PA(PA)
Q=P*cos(PA*!PI/90)
U=P*sin(PA*!PI/90)
window,1,xsize=400 ,ysize=1100
!P.multi=[0,1,5]
plot,wave,flux,xst=1,charsize=2
plot,wave,Q,xst=1,charsize=2
plot,wave,U,xst=1,charsize=2
plot,wave,P,xst=1,charsize=2
plot,wave,PA,xst=1,charsize=2
end
