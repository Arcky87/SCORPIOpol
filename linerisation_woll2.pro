function linerisation_WOLL2,ima,disp,PARAM=param
a=size(ima)
a=size(ima)  & Nx=a(1) & Ny=a(2)
;определение пределов линеаризации
D=total(disp,1)/a(2) & Ndeg=N_elements(D)-1
lambda=0 & for j=0,Ndeg do lambda=lambda+D(j)*findgen(Nx)^j
;d_lambda=fix(D(1))
;N_lin=FIX((lambda(Nx-1)-lambda(0))/d_lambda)
;lambda_0=FIX(lambda(0))
;if not(keyword_set(param)) then param=[lambda_0,d_lambda,N_lin]
lambda_lin=findgen(param(2))*param(1)+param(0)

ima_lin=fltarr(param(2),Ny)
;линеаризация по высоте щели
for y=0,Ny-1 do begin
lambda=0 & for j=0,Ndeg do lambda=lambda+disp(y,j)*findgen(Nx)^j
ima_lin(*,y)=INTERPOL(ima(*,y),lambda,lambda_lin)

endfor
return,ima_lin
end