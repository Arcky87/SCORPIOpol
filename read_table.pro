;reading table
function READ_TABLE,file
a=' '
openr,UNIT,file,/GET_LUN
readf,UNIT,a
table=a
while NOT EOF(UNIT) do begin
readf,UNIT,a
table=[table,a]
ENDWHILE
close,UNIT
FREE_LUN,UNIT
return,table
end



