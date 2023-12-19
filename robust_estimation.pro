function robust_estimation,vector,WX=wx,TRESH=tresh
if not(keyword_set(tresh)) then tresh=3
if not(keyword_set(wx)) then wx=20
N=N_elements(vector)
Npos=N/wx

xpos=findgen(Npos)*Wx+wx/2
avg_vector=fltarr(Npos,3)
avg_vector(*,0)=xpos
for k=0,Npos-1 do begin
robomean,vector(xpos(k)-wx/2:xpos(k)+wx/2),tresh,0.5,avg_value,rms_value
avg_vector(k,1)=avg_value
avg_vector(k,2)=rms_value
endfor
return,avg_vector
end
