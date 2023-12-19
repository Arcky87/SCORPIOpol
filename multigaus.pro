;+
; NAME:
;       MULTIGAUS
;
; PURPOSE:
;       Computes a non-linear least-squares fit to a function Y(x)
;       with numerous of gaussians (single and double-components) without continuum.
;        Number of gaussian is equialent   to the number of elements in POS vector.
;
; CALLING SEQUENCE:
;       Result=MULTIGAUS ( X, Y, POS, FWHM=FWHM [,/FIXpos][,/FIXfwhm], [,/PLOT=PLOT])
;                          [,/SILENT] [ABSORP=ABSORP] [,/LIM_FWHM=], [,/DOUBL_LEN]
;                          [,/DOUBL_RATIO=] [,/FIXRATIO=],YFIT=YFIT,sigma=sigma)
;
; INPUTS:
;       X:  A row vector of independent variables.
;       Y:  A row vector containing the dependent variable.
;     POS:  A row vector containing the preliminary values of centers of the gaussian
;
; KEYWORDS:
;    FWHM:  a row vector containing the  preliminary values of  FWHM's value of an each gaissian
;  FIXpos:  if this keyword set then gaissians center will be fixed as POS-values
; FIXfwhm:  if this keyword set then gaissians FWHM's will be fixed as FWHM-keyword value
;    PLOT:  if it set then fitting results will be plotting to a window,10
;  SILENT:  if set then Error messages will be skipping
;  ABSORP:  a row vector, if absorp(i)>0 then   amplitude(i)<0 will be using (for absorption profiles)
;  LIM_FWHM: a two-component vector containing the maximal and minimal limits for FWHM
;
; KEYWORDS for doublets
; DOUBL_LEN:  a row vector with distances of the second gaussian components. The size of this vector
;           must be equalent to the number of elements in POS vector. Set 0 for single gaussian.
; DOUBL_RATIO: a row vector with preliminary amplitude ratio of the secondary components in doublets
;              The size of this vector is  equalent of the POS size.
; FIXRATIO  :if this keyword set then doublet amplitudes will be fixed as DOUBLE_RATIO
; SIGMA : A structure variable that will contain a vector of standard deviations for the elements of the output parameters
;
; OUTPUTS:
;       Returns a structure-vector of fitting values:
;       {GAUS,MAX:float,CENTER:float.,FWHM:float,FLUX:float,RATIO:float}
;       Therefore the parameters of the i-gaussian are:
;       RESULT(i).MAX -- amplitude
;       RESULT(i).CENTER -- center position
;       RESULT(i).FWHM -- FWHM
;       RESULT(i).FLUX -- ingegral value
;       RESULT(i).RATIO-- relative amplitude of the second component of doublets
;
;
; COMMON BLOCKS:
;       CURPRF, PARAMETERS, SIL, FITCOR, DOUBLETS
;
; MODIFICATION HISTORY:
;        A.V. Moiseev SAO RAS, 24 Oct 2002


; function for multi-gaussian curvefit
PRO gaumul, X, A, F, pder
COMMON CURPRF,Ycur
COMMON PARAMETERS,Aold,N,recfit,nfit,fix1,fix2
common SIL,SILENT
common FITCOR, dfw,da,dpos
common DOUBLETS,d_len,dbl_rec,dbl_rat,fix3,N_dbl

 if (nfit-N_dbl)/3. ne N then begin
           An=Aold
           An(recfit)=A
 endif else An=A


;print,an
 ; correction for amplitudes

aamp=an(0:N-1)
bad=where( (aamp-da(*,0)) lt 0,counts)
     if counts gt 0 then aamp(bad)=da(bad,0)
bad=where( (aamp-da(*,1)) gt 0,counts)
if counts gt 0 then aamp(bad)=da(bad,1)
an(0:N-1)=aamp


 ; correction on FWHM limits
 IF fix2 eq 0 then begin
   afw=abs(an(2*N:3*N-1))
   if dfw(0) ne -1 then begin
     bad=where( afw lt dfw(0),counts)
     if counts gt 0 then afw(bad)=dfw(0)
     bad=where( afw gt dfw(1),counts)
     if counts gt 0 then afw(bad)=dfw(1)
   endif
   an(2*N:3*N-1)=afw
ENDIF

; return new values of the fitting parameters
a=an(recfit)


 F=0.
 for i=0,N-1 do F=F+an(i)*exp(-(x-an(N+i))^2/(2*an(2*N+i)^2))

; for doublets
 if N_dbl gt 0 then begin
  for i=0,N_dbl-1 do begin
   j=dbl_rec(i)
   F=F+an(j)*an(3*N+i)*exp(-(x-an(N+j)-d_len(j))^2/(2*an(2*N+i)^2))
  endfor
 endif

    rec=where( finite(f) eq 0,counts)

    ; removing infinite values

    if counts GT 0 then begin
       f(rec)=Ycur(rec)
       if SILENT eq 0 then print, 'Fitting error!'
       a(*)=0
     endif

  IF N_PARAMS() GE 4 THEN BEGIN ;       If the procedure is called with four parameters,
                          ;calculate the partial derivatives.
  pder=make_array(n_elements(x),Nfit)

  for i=0,N-1 do begin
   j=N+i & k=2*N+i & e=exp(-0.5*(x-an(j))^2/an(k)^2)
   pder(*,i)=e
   if fix1 eq 0 then pder(*,j)=an(i)*e*(x-an(j))/an(k)^2
   if (fix2 eq 0)and (fix1 eq 0) then pder(*,k)=an(i)*e*(x-an(j))^2/an(k)^3
   if (fix2 eq 0)and (fix1 eq 1) then pder(*,j)=an(i)*e*(x-an(j))^2/an(k)^3
  endfor

; for doublets
 if N_dbl gt 0 then begin
  for i=0,N_dbl-1 do begin
   i0=dbl_rec(i)
   j=N+i0 & k=2*N+i0 & l=3*N+i
   e2=exp(-0.5*(x-an(j)-d_len(i0))^2/an(k)^2)
   pder(*,i0)=pder(*,i0)+e2*an(l)
   if fix1 eq 0 then pder(*,j)=pder(*,j)+an(i0)*an(l)*e2*(x-an(j)-d_len(i0))/an(k)^2
   if (fix2 eq 0)and (fix1 eq 0) then pder(*,k)=pder(*,k)+an(i0)*e2*an(l)*(x-an(j)-d_len(i0))^2/an(k)^3
   if (fix2 eq 0)and (fix1 eq 1) then pder(*,j)=pder(*,j)+an(i0)*e2*an(l)*(x-an(j)-d_len(i0))^2/an(k)^3

    if fix3 eq 0 then begin
     j=N+i & k=2*N+i
      if (fix2 eq 0)and (fix1 eq 0) then pder(*,l)=an(i0)*e2
      if (fix2 eq 0)and (fix1 eq 1) then pder(*,k)=an(i0)*e2
      if (fix2 eq 1)and (fix1 eq 0) then pder(*,k)=an(i0)*e2
     if (fix2 eq 1)and (fix1 eq 1) then pder(*,j)=an(i0)*e2
    endif
  endfor
endif

  ; removing infinite values
     rec=where( finite(pder) eq 0,counts)
    if counts GT 0 then begin
      pder(rec)=0.
      if SILENT eq 0 then print, 'Fitting error!'
      f=Ycur
      a(*)=0.
    endif
  ENDIF

END


;1D fitting by multi gaussian
; A.V. Moiseev (SAO RAN) May 2001

Function multigaus,x,y,pos,fwhm=fwhm,plot=plot,fixpos=fixpos,fixfwhm=fixfwhm,$
         SILENT=SILENT,absorp=absorp,lim_fwhm=lim_fwhm,SIGMA=SIGMA,$
         doubl_len=doubl_len,doubl_ratio=doubl_ratio,fixratio=fixratio,YFIT=YFIT
COMMON CURPRF,Ycur
COMMON PARAMETERS,Aold,N,recfit,nfit,fix1,fix2
common sil,SILEN
common FITCOR, dfw,da,dpos
common DOUBLETS,d_len,dbl_rec,dbl_rat,fix3,N_dbl
if not(keyword_set(silent)) then silen=0 else silen=1

out=0
mul=2.35482 ; convert from sigma to FWHM

if n_params() ne 3 then begin
                          message,' Syntax - Result=MULTIGAUS( X, Y, POS [, FWHM=FWHM, /FIXpos, /FIXfwhm, /FIXRATIO, /PLOT,',/contin
                          print, '                ,/SILENT, /ABSORP ,  LIM_FWHM=LIM_FWHM, DOUBL_LEN=DOUBL_LEN, DOUBL_RATIO=DOUBL_RATIO, YFIT=YFIT, SIGMA=SIGMA])
                          print
                          print,'RESULT and SIGMA are structures with fileds:.MAX; .CENTER; .FWHM and  .FLUX'
                          return,0.
                        endif

N=n_elements(pos)
nx=n_elements(x)
IF not(keyword_set(fwhm)) then fwhm=make_array(N,value=5.)
if keyword_set(fixpos) then fix1=1 else fix1=0
if keyword_set(fixfwhm) then fix2=1 else fix2=0
if keyword_set(lim_fwhm) then dfw=lim_fwhm/mul else dfw=[-1.,-1]
if n_elements(y)   ne nx then begin message,'X and Y must be equivalent!' & return,0. & end
if N ne n_elements(fwhm) then begin message,'FWHM and POS must be equivalent!' & return,0. & end

; doublet setting
fix3=-1
N_dbl=0
if keyword_set(doubl_len) then begin
  d_len=doubl_len
  dbl_rec=where(doubl_len ne 0,N_dbl)
  if keyword_set(doubl_ratio) then dbl_rat=doubl_ratio else dbl_rat=make_array(N,val=1.)
  if keyword_set(fixratio) then fix3=1 else if N_dbl gt 0 then fix3=0
endif ; doublets

; normalization:
norm=max(abs(y))
if norm eq 0 then norm=1.
y=y/norm

 sig=fwhm/mul
 xi=pos
 Ycur=Y

 ; set preliminary values of parameters
 A=fltarr(3*N+N_dbl)
 ; amplitudes
 if not(keyword_set(absorp)) then absorp=bytarr(N)
 for i=0,N-1 do begin
           rec=where( (x ge pos(i)-fwhm(i)) and  (x le pos(i)+fwhm(i)),count)
           if count gt 0 then $
              if absorp(i) gt 0  then a(i)=min(y(rec)) else a(i)=max(y(rec))

 endfor


 A(N:2*N-1)=float(pos) ; positions
 A(2*N:3*N-1)=FWHM/mul ; sigma
 if N_dbl gt 0 then a(3*N:3*N+N_dbl-1)=dbl_rat(dbl_rec) ; doublets ratio

 a1=0.3*a(0:N-1); set boundaries for amplitude
 a2=3*a(0:N-1)  ;
 da=fltarr(N,2)
 for i=0,N-1 do begin
    if absorp(i) gt 0 then da(i,*)=[a2(i),a1(i)] else da(i,*)=[a1(i),a2(i)]
 endfor
 ;if keyword_set(absorp) then da=reform([a2,a1],N,2) else da=reform([a1,a2],N,2)

 Aold=A
 recfit=indgen(N)
 if not(keyword_set(fixpos)) then recfit=[recfit,N+indgen(N)]
  if not(keyword_set(fixfwhm)) then recfit=[recfit,2*N+indgen(N)]
  if (not(keyword_set(fixratio)) )and (N_dbl gt 0) then recfit=[recfit,3*N+dbl_rec]

 nfit=n*(1+(1 xor fix1)+(1 xor fix2))+N_dbl*(fix3 xor 1)

 W=make_array(nx, value=1.);./nx)
 A=Aold(recfit)
 itmax=50

 Yfit =CURVEFIT(x,y,W,A,Sig1A,function_name='gaumul',itmax=itmax,tol=1e-4,iter=iter,/double)
 Aold(recfit)=A
 A=Aold

 sig1=fltarr(3*N+N_dbl)
 sig1(recfit)=sig1A

y=y*norm
A(0:N-1)=A(0:N-1)*norm
sig1(0:N-1)=sig1(0:N-1);*norm

  if keyword_set(plot) then begin
  window,10,xs=500
  plot,x,y,psym=4,xst=1,yst=1,symsize=0.3
  sum=0.
  for i=0,N-1 do sum=sum+gaussian(x,[a(i),a(N+i),a(2*N+i),0])
 ; for doublets
  if N_dbl gt 0 then begin
   for i=0,N_dbl-1 do begin
    j=dbl_rec(i)

    sum=sum+a(j)*a(3*N+i)*exp(-(x-a(N+j)-d_len(j))^2/(2*a(2*N+j)^2))
    endfor
   endif

   oplot,x,sum
                       endif
  out=replicate({gaus,max:-1.,center:-1.,fwhm:-1.,flux:-1.,ratio:-1.},N)
  sigma=replicate({sigma,max:-1.,center:-1.,fwhm:-1.,flux:-1.,ratio:-1.},N)
  for i=0,n-1 do begin
     out(i).max=a(i)
     out(i).center=a(N+i)
     out(i).fwhm=(a(2*N+i)*mul)
     out(i).flux=a(i)*sqrt(2*!pi)*(a(2*N+i))

     sigma(i).max=sig1(i)
     sigma(i).center=sig1(N+i)
     sigma(i).fwhm=(sig1(2*N+i)*mul)
     s1=sqrt(2*!pi)*(a(2*N+i))*sig1(i)
     s2=sqrt(2*!pi)*a(i)*sig1(2*N+i)
     sigma(i).flux=sqrt(s1^2+s2^2)

   endfor
  out(*).ratio=0
  if (N_dbl gt 0) then out(dbl_rec).ratio=a(3*N:3*N+N_dbl-1)

  return,out
END
