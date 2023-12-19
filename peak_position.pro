function peak_position,vector,xpk,eps
a=size(vector)

x=findgen(a(2))

f=goodpoly(x(xpk-eps:xpk+eps),vector(xpk-eps:xpk+eps),2,3,FIT)

xpk=-f(1)/f(2)/2
return,xpk
end