function create_continuum,spectra,REP=rep
if not(keyword_set(rep)) then rep=4
N=N_elements(spectra)
for k=0,rep-1 do begin
cont=LOWESS(findgen(N),spectra,N/6,3,3)
robomean,spectra-cont,3,0.5,avg_cont,rms_cont
R=where(spectra-cont lt -rms_cont, ind)
if ind gt 0 then spectra(R)=cont(R)
endfor
return,spectra
end