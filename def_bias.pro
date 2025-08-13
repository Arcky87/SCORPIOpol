function def_bias,ima,PLOT=plot
; определение уровня фона (bias)
Nhist=100
bin_hist=3
x_hist=findgen(Nhist+1)*bin_hist+1050-bin_hist*Nhist/2
y_hist=histogram(ima,min=x_hist(0),max=x_hist(Nhist),bin=bin_hist)
gau=gaussfit(x_hist,y_hist,G)
bias=g(1)
print, 'calculated bias level: ', bias
if keyword_set(plot) then begin
window,2
plot,x_hist,y_hist,psym=10
oplot,x_hist,gau
   wait, 1
endif
return,bias
end
