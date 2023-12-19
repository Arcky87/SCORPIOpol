; $Id: poly_fit.pro,v 1.13.4.1 2002/09/23 21:07:44 chris Exp $
;+
; NAME:
;	POLY_FIT
;
; PURPOSE:
;	Perform a least-square polynomial fit with optional error estimates.
;
;	This routine uses matrix inversion.  A newer version of this routine,
;	SVDFIT, uses Singular Value Decomposition.  The SVD technique is more
;	flexible, but slower.
;
; CATEGORY:
;	Curve fitting.
;
; CALLING SEQUENCE:
;   Result = POLY_FIT(X, Y, Degree)
;
; INPUTS:
;   X:  The independent variable vector.
;
;   Y:  The dependent variable vector, should be same length as x.
;
;   Degree: The degree of the polynomial to fit.
;
; OUTPUTS:
;	POLY_FIT returns a vector of coefficients with a length of NDegree+1.
;
; KEYWORDS:
;   CHISQ:   Sum of squared errors divided by MEASURE_ERRORS if specified.
;
;   COVAR:   Covariance matrix of the coefficients.
;
;	DOUBLE:  if set, force computations to be in double precision.
;
;   MEASURE_ERRORS: Set this keyword to a vector containing standard
;       measurement errors for each point Y[i].  This vector must be the same
;       length as X and Y.
;
;     Note - For Gaussian errors (e.g. instrumental uncertainties),
;        MEASURE_ERRORS should be set to the standard
; 	     deviations of each point in Y. For Poisson or statistical weighting
; 	     MEASURE_ERRORS should be set to sqrt(Y).
;
;   SIGMA:   The 1-sigma error estimates of the returned parameters.
;
;     Note: if MEASURE_ERRORS is omitted, then you are assuming that
;           your model is correct. In this case,
;           SIGMA is multiplied by SQRT(CHISQ/(N-M)), where N is the
;           number of points in X. See section 15.2 of Numerical Recipes
;           in C (Second Edition) for details.
;
;   STATUS = Set this keyword to a named variable to receive the status
;          of the operation. Possible status values are:
;          0 for successful completion, 1 for a singular array (which
;          indicates that the inversion is invalid), and 2 which is a
;          warning that a small pivot element was used and that significant
;          accuracy was probably lost.
;
;    Note: if STATUS is not specified then any error messages will be output
;          to the screen.
;
;   YBAND:	1 standard deviation error estimate for each point.
;
;   YERROR: The standard error between YFIT and Y.
;
;   YFIT:   Vector of calculated Y's. These values have an error
;           of + or - YBAND.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	None.
;
; MODIFICATION HISTORY:
;	Written by: George Lawrence, LASP, University of Colorado,
;		December, 1981.
;
;	Adapted to VAX IDL by: David Stern, Jan, 1982.
;       Modified:    GGS, RSI, March 1996
;                    Corrected a condition which explicitly converted all
;                    internal variables to single-precision float.
;                    Added support for double-precision inputs.
;                    Added a check for singular array inversion.
;		     SVP, RSI, June 1996
;                     Changed A to Corrm to match IDL5.0 docs.
;                    S. Lett, RSI, December 1997
;                     Changed inversion status check to check only for
;                     numerically singular matrix.
;                    S. Lett, RSI, March 1998
;                     Initialize local copy of the independent variable
;                     to be of type DOUBLE when working in double precision.
;       CT, RSI, March 2000: Changed to call POLYFITW.
;       CT, RSI, July-Aug 2000: Removed call to POLYFITW,
;                   added MEASURE_ERRORS keyword,
;                   added all other keywords (except DOUBLE),
;                   made output arguments obsolete.
;-

FUNCTION POLY_FIT, x, y, ndegree, $
	yfit_old, yband_old, yerror_old, corrm_old, $     ; obsolete arguments
	CHISQ=chisq, $
	COVAR=covar, $
	DOUBLE=double, $
	MEASURE_ERRORS=measure_errors, $
	SIGMA=sigma, $
	STATUS=status, $
	YBAND=yband, $
	YERROR=yerror, $
	YFIT=yfit

    COMPILE_OPT idl2

	ON_ERROR,2		;RETURN TO CALLER IF ERROR

	n = N_ELEMENTS(x)
	IF (n NE N_ELEMENTS(y)) THEN MESSAGE, $
		'X and Y must have same number of elements.'
	m = ndegree + 1	; # of elements in coeff vec

	double = (N_ELEMENTS(double) GT 0) ? KEYWORD_SET(double) : $
		(SIZE(x,/TNAME) EQ 'DOUBLE') OR (SIZE(y,/TNAME) EQ 'DOUBLE')

	no_weight = (N_ELEMENTS(measure_errors) EQ 0)
	sdev = 1d
	IF NOT no_weight THEN sdev = sdev*measure_errors
	sdev2 = sdev^2

	; construct work arrays
	covar = DBLARR(m,m) ; least square matrix, weighted matrix
	b = DBLARR(m)	; will contain sum weights*y*x^j
	z = DBLARR(n) + 1	; basis vector for constant term
	wy = DOUBLE(y)/sdev2
	yfit = DBLARR(n)
	yband = DBLARR(n)
	yerror = !VALUES.D_NAN


	covar[0,0] = no_weight ? n : TOTAL(1/sdev2)
	b[0] = TOTAL(wy)


	FOR p = 1L,2*ndegree DO BEGIN	; power loop
		z = TEMPORARY(z)*x	; z is now x^p
		IF p LT m THEN b[p] = TOTAL(wy*z)	; b is sum weights*y*x^j
		sum = TOTAL(z/sdev2)
		FOR j = 0 > (p-ndegree), ndegree < p DO covar[j,p-j] = sum
	ENDFOR ; end of p loop, construction of covar and b


	covar = INVERT(TEMPORARY(covar), status)


	IF NOT ARG_PRESENT(status) THEN BEGIN
		CASE status OF
		1: MESSAGE, "Singular matrix detected."
		2: MESSAGE,/INFO, "Warning: Invert detected a small pivot element."
		ELSE:
		ENDCASE
	ENDIF

	if (status eq 1) then begin
	    result = !VALUES.D_NAN
	    goto, done
	endif


	result = (TEMPORARY(b) # covar)  ; construct coefficients


	; compute optional output parameters.

	; one-standard deviation error estimates, init
	yfit = TEMPORARY(yfit) + result[ndegree]
	FOR k = ndegree-1L, 0, -1 DO yfit = result[k] + TEMPORARY(yfit)*x  ; sum basis vectors

	chisq = TOTAL((yfit-y)^2/sdev2)
	variance = covar[lindgen(M)*(M+1)]
	IF no_weight THEN variance = variance*chisq/(n-m)
	sigma = SQRT(ABS(variance))


	; Experimental variance estimate, unbiased
	var = (n GT m) ? TOTAL((yfit-y)^2 )/(n-m) : 0d

	yerror = SQRT(var)
	z = DBLARR(n) + 1d
	yband = TEMPORARY(yband) + covar[0,0]
	FOR p=1L,2*ndegree DO BEGIN	; compute correlated error estimates on y
		z = TEMPORARY(z)*x		; z is now x^p
		sum = 0
		FOR j=0 > (p - ndegree), ndegree  < p DO sum = sum + covar[j,p-j]
		yband = TEMPORARY(yband) + sum * z  ; add in all the error sources
	ENDFOR	; end of p loop

	yband = TEMPORARY(yband)*var
	IF (MIN(yband) LT 0) OR (MIN(FINITE(yband)) EQ 0) THEN BEGIN
		status = 3
		IF NOT ARG_PRESENT(status) THEN MESSAGE, $
			'Undefined (NaN) error estimate encountered.'
	ENDIF ELSE yband = SQRT( TEMPORARY(yband) )

done:

    ; If necessary, convert all results to single precision.
    if (not double) then begin
        chisq = FLOAT(chisq)
        covar = FLOAT(covar)
        result = FLOAT(result)
        sigma = FLOAT(sigma)
        var = FLOAT(var)   ; needed for corrm_old below
        yfit = FLOAT(yfit)
        yband = FLOAT(yband)
        yerror = FLOAT(yerror)
    endif

; fill in obsolete arguments, if necessary
	IF (N_PARAMS() GT 3) THEN BEGIN
		corrm_old = covar*var   ; convert to correlation matrix
		yerror_old = yerror
		yband_old = yband
		yfit_old = yfit
	ENDIF

	RETURN, result
END
