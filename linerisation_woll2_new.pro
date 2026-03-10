;+
; NAME:
;       MAKE_BINS_IDL
;
; PURPOSE:
;       Given a series of wavelength points, find the edges and widths
;       of corresponding wavelength bins. Ported from the Python version
;       of spectres.
;
; CALLING SEQUENCE:
;       result = make_bins_idl(wavs)
;
; INPUTS:
;       wavs: A 1D array of wavelength points.
;
; OUTPUTS:
;       A structure containing the bin edges and widths.
;-
FUNCTION make_bins_idl, wavs
  n = N_ELEMENTS(wavs)
  edges = DBLARR(n + 1)

  edges[1:n-1] = (wavs[0:n-2] + wavs[1:n-1]) / 2.0
  edges[0] = wavs[0] - (wavs[1] - wavs[0]) / 2.0
  edges[n] = wavs[n-1] + (wavs[n-1] - wavs[n-2]) / 2.0

  widths = edges[1:*] - edges[0:*-1]

  RETURN, {edges:edges, widths:widths}
END

;+
; NAME:
;       SPECTRES_IDL
;
; PURPOSE:
;       Resamples a spectrum onto a new wavelength basis, conserving flux.
;       Ported from the Python version of spectres.
;
; CALLING SEQUENCE:
;       new_fluxes = spectres_idl(new_wavs, spec_wavs, spec_fluxes, FILL=fill)
;
; INPUTS:
;       new_wavs:    Array containing the new wavelength sampling desired.
;       spec_wavs:   1D array containing the current wavelength sampling.
;       spec_fluxes: 1D array containing spectral fluxes at spec_wavs.
;
; KEYWORD PARAMETERS:
;       FILL:   Where new_wavs extends outside the wavelength range in
;               spec_wavs, this value will be used as a filler.
;
; OUTPUTS:
;       An array of resampled flux values.
;-
FUNCTION spectres_idl, new_wavs, spec_wavs, spec_fluxes, FILL=fill

  ; Set default fill value if not provided
  IF NOT KEYWORD_SET(fill) THEN fill = 0.0

  ; Make arrays of edge positions and widths for the old and new bins
  old_bins = make_bins_idl(spec_wavs)
  old_edges = old_bins.edges
  old_widths = old_bins.widths

  new_bins = make_bins_idl(new_wavs)
  new_edges = new_bins.edges

  ; Generate output array
  new_fluxes = DBLARR(N_ELEMENTS(new_wavs))

  start = 0L
  stop = 0L

  ; Calculate new flux values, looping over new bins
  FOR j = 0, N_ELEMENTS(new_wavs) - 1 DO BEGIN

    ; Add filler values if new_wavs extends outside of spec_wavs
    IF (new_edges[j] LT old_edges[0]) OR (new_edges[j+1] GT old_edges[N_ELEMENTS(old_edges)-1]) THEN BEGIN
      new_fluxes[j] = fill
      CONTINUE
    ENDIF

    ; Find first old bin which is partially covered by the new bin
    WHILE old_edges[start+1] LE new_edges[j] DO start = start + 1

    ; Find last old bin which is partially covered by the new bin
    WHILE old_edges[stop+1] LT new_edges[j+1] DO stop = stop + 1

    ; If new bin is fully inside an old bin, start and stop are equal
    IF stop EQ start THEN BEGIN
      new_fluxes[j] = spec_fluxes[start]
    ENDIF ELSE BEGIN
      ; Otherwise, multiply the first and last old bin widths by P_ij

      start_factor = (old_edges[start+1] - new_edges[j]) / (old_edges[start+1] - old_edges[start])
      end_factor = (new_edges[j+1] - old_edges[stop]) / (old_edges[stop+1] - old_edges[stop])

      temp_old_widths = old_widths
      temp_old_widths[start] = temp_old_widths[start] * start_factor
      temp_old_widths[stop] = temp_old_widths[stop] * end_factor

      ; Sum the flux from all overlapping old bins
      f_widths = temp_old_widths[start:stop] * spec_fluxes[start:stop]
      new_fluxes[j] = TOTAL(f_widths) / TOTAL(temp_old_widths[start:stop])

    ENDELSE
  ENDFOR

  RETURN, new_fluxes
END


function linerisation_WOLL2,ima,disp,PARAM=param
a=size(ima)
a=size(ima)  & Nx=a(1) & Ny=a(2)
;определение пределов линеаризации
D=total(disp,1)/a(2) & Ndeg=N_elements(D)-1
lambda=0
for j=0,Ndeg do lambda=lambda+D(j)*findgen(Nx)^j
;d_lambda=fix(D(1))
;N_lin=FIX((lambda(Nx-1)-lambda(0))/d_lambda)
;lambda_0=FIX(lambda(0))
;if not(keyword_set(param)) then param=[lambda_0,d_lambda,N_lin]
lambda_lin=findgen(param(2))*param(1)+param(0)

ima_lin=fltarr(param(2),Ny)
;линеаризация по высоте щели
for y=0,Ny-1 do begin
lambda=0 & for j=0,Ndeg do lambda=lambda+disp(y,j)*findgen(Nx)^j
;ima_lin(*,y)=INTERPOL(ima(*,y),lambda,lambda_lin)
ima_lin(*,y) = spectres_idl(lambda_lin, lambda, ima(*, y))

endfor
return,ima_lin
end
