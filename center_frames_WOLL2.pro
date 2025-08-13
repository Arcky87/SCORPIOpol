function center_frames_WOLL2,eta
a=size(eta)
eta=eta-def_bias(eta)
w=20 ; default = 20
vector=total(eta(a(1)/2-w:a(1)/2+w,*),1)/(2*w+1)
R=where(vector lt 0,ind) & if ind gt 0 then vector(R)=0
xpk=find_peaks(vector,W=20,TRESH=5,/plot)  ;20 100
tmp=reform(xpk,3,4) & yc=fltarr(4)  & yc(*)=tmp(1,*)
return,yc
end
dir='/data6/SCORPIO/sspol_pipeline_v2023.8/Mrk1018/'
;dir='h:\red_data.pol\Arp102b_131103\'
eta=readfits(dir+'eta_i.fts')
print,center_frames_WOLL2(eta(*,*,0))
END
