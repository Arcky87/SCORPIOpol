pro fi_peak,xpl,ypl,pkcut,ipix,xpk,ypk,bkpk,ipk  
;      
;+      
; NAME:      
;        FI_PEAK     
; PURPOSE:      
;        Find peaks in a lineout      
; CATEGORY:      
;        reduction MPFS-data      
; CALLING SEQUENCE:      
;        peak_find      
; INPUTS:      
;        xpl:    x array      
;        ypl:    y array      
;        pkcut:  cutoff value for determining peak      
; KEYWORD PARAMETERS:      
;        no 
; OUTPUTS:      
;        ipix:   pixel location of peaks, array      
;        xpk:    x location of peaks, array      
;        ypk:    y location of peaks, array      
;        bkpk:   background at peak locations, array      
; COMMON BLOCKS:      
;        None      
; SIDE EFFECTS:      
;        None.      
; RESTRICTIONS:      
;        None.      
; PROCEDURE:      
;        no      
; MODIFICATION HISTORY: 
;	Written by Victor Afanasiev, Special Astrophysical Observatory RAS, Jul 1999     
;-      
xpk=fltarr(500)      
ypk=fltarr(500)      
bkpk=fltarr(500)      
ipix=intarr(500)      
asiz=size(xpl)      
ma=asiz(1)      
ipk=0      
isc=7      
nscan=10  
;(isc)>1.          
   mm=1      
   is=isc      
while (is lt ma-isc) and (is ge isc) do begin
  ib1=(is-nscan)>0      
  ib2=(is+nscan)<(ma-1)      
  ych=ypl(ib1:ib2)      
     back2=min(ych)*1.01      
     ymx=max(ypl(is-isc:is+isc))      
     ytmp=ypl(is)-back2      
  if (ypl(is) eq ymx) and (ytmp gt pkcut) then begin
;      
; ....Have found line peak      
;      
     ipix(ipk)=is      
     xpk(ipk)=xpl(is)      
     ypk(ipk)=ypl(is)      
     bkpk(ipk)=back2      
     ipk=ipk+1      
     is=is+mm*(isc+1)      
  endif else is=is+mm*1      
endwhile
if ipk eq 0 then begin
   ipk=1      
   xpk=[0.]      
   ypk=[0.]      
   bpk=[0.]      
   ipix=[0]      
endif else begin
   xpk=xpk(0:ipk-1)      
   ypk=ypk(0:ipk-1)      
   bkpk=bkpk(0:ipk-1)      
   ipix=ipix(0:ipk-1)      
endelse      
return      
end
