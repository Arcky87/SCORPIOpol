function expand_traectory,tra,SCALE=scale
;Добавление 
a=size(tra)
print,'size tra', size(tra)
tra_new=fltarr(a(1),(a(2)-1)*scale+1)
y_in=findgen(a(2))
y_out=findgen((a(2)-1)*scale+1)/scale
for x=0,a(1)-1 do tra_new(x,*)=INTERPOL(tra(x,*),y_in,y_out)
return,tra_new
end
