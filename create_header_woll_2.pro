;crea header
function create_header_WOLL_2,Nx,Ny,Nray,Nobj,Ncube,lambda_0,d_lambda
mkhdr,header,0
sxaddpar,header,'BITPIX',-32,after='SIMPLE'
sxaddpar,header,'NAXIS',4,' NUMBER OF AXES'
sxaddpar,header,'NAXIS1',Nx,' LENGTH OF SPECTRA',after='NAXIS'
sxaddpar,header,'NAXIS2',Ny,' HIGH OF SPECTRA',after='NAXIS1'
sxaddpar,header,'NAXIS3',Nray,' POLARIZED RAY',after='NAXIS12
sxaddpar,header,'NAXIS4',Nobj,' NUMBER EXPOSURES OF OBJECT'
sxaddpar,header,'NAXIS5',Ncube,' NUMBER OF DATA CUBE
sxaddpar,header,'CRVAL1',lambda_0,after='NAXIS1',' INITIAL WAVELENGTH ,A'
sxaddpar,header,'CDELT1',d_lambda,after='CRVAL1',' DISPERSION A/PX'
sxaddpar,header,'DATE-OBS','  ',' DATE OF OBSERVATION'
sxaddpar,header,'PROG-ID','  ',' PROGRAMM OF OBSERVATION'
sxaddpar,header,'AUTHOR','  ',' AUTHOR OF PROGRAMM'
sxaddpar,header,'OBSERVER',' ',' OBSERVERS'
sxaddpar,header,'DIR',' ',' DIRECTORY OF INITIAL DATA'
sxaddpar,header,'OBSERVAT','SAO RAS ',' OBSERVATORY

sxaddpar,header,'TELESCOP','BTA 6-meter ',' TELESCOPE NAME
sxaddpar,header,'FOCUS',0,' focus of telescope (mm)
sxaddpar,header,'DETECTOR','E2V CCD42-90 red',' DETECTOR
sxaddpar,header,'BINNING',' ',' BINNING
sxaddpar,header,'RATE',0,' readout rate (KPix/sec)
sxaddpar,header,'GAIN',0,' gain, electrons per adu
sxaddpar,header,'NODE','A',' output node (A, B, AB)
sxaddpar,header,'INSTRUME','SCORPIO-2',' INSTRUMENT
sxaddpar,header,'IMSCALE',' ',' image scale ("/Pix x "/Pix)
sxaddpar,header,'CAMFOCUS',0,' FOCUS OF CAMERA, mm
sxaddpar,header,'COLFOCUS',0,' FOCUS OF COLLIMATOR, mm
sxaddpar,header,'MODE',' ',' mode of instrument
sxaddpar,header,'DISPERSE',' ',' GRATING
sxaddpar,header,'SLITWID',0,' slit width in arcsec
sxaddpar,header,'SLITMASK',' ',' slit mask
sxaddpar,header,'FILTERS',' ',' filter name, filter wheel number 3
sxaddpar,header,'FILTPOS1',0,' position of wheel number 1
sxaddpar,header,'FILTPOS2',0,' position of wheel number 2

for k=1,3 do begin
sxaddpar,header,'NAME'+string(k,format='(I1)'),'  ',' NAME OF TARGET'
if  k eq 1 then sxaddpar,header,'Z','  ',' REDSHIFT'
sxaddpar,header,'CUBE'+string(k,format='(I1)'),'  ',' NAME DATA CUBE'
sxaddpar,header,'START'+string(k,format='(I1)'),'  ',' START OF EXPOSURE'
sxaddpar,header,'NUMEXP'+string(k,format='(I1)'),'  ',' NUMBER OF EXPOSURE'
sxaddpar,header,'EXPTIME'+string(k,format='(I1)'),'  ',' TIME OF EXPOSURE'
sxaddpar,header,'RA'+string(k,format='(I1)'),'  ',' APPARENT RIGTH ASCESSION'
sxaddpar,header,'DEC'+string(k,format='(I1)'),'  ',' APPARENT DECLINATION'
sxaddpar,header,'A'+string(k,format='(I1)'),'  ',' AZIMUT'
sxaddpar,header,'Z'+string(k,format='(I1)'),'  ',' ZENIT DISTANCE '
sxaddpar,header,'PA'+string(k,format='(I1)'),'  ',' POSITION ANGLE OF SLIT'
endfor
fxaddpar,header,'COMMENT','************ DEVICE ***********',before='TELESCOP'
fxaddpar,header,'COMMENT','************ TARGET: OBJECT ***********',before='name1'
fxaddpar,header,'COMMENT','** TARGET: UNPOLARIZED STANDARD STAR **',before='name2'
fxaddpar,header,'COMMENT','*** TARGET: POLARISED STANDARD STAR ***',before='name3'
return,header
end
