function readtable, file, Nc

N=numlines(file)
tab=dblarr(Nc,N)
openr, 1, file
readf, 1, tab
close, 1

return, tab

end