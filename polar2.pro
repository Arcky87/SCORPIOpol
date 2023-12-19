;BATH file SCORPIO-2 data  reduction in spectropolarimetry WOLL2 mode
;CALLING SEQUENCE:
;
;	@polar2
;

;definition PATH
!PATH=!PATH+';'+'/home/elias/SCORPIO/sppol_pipeline_v2023.8/WOLLASTON-2.lib/'
.compile poly_fit.pro
;.compile stdev.pro
.compile WOLL2_v1.pro
.compile ViewPol_2.pro
.compile frames_WOLL2.pro
.compile lowess.pro
.compile def_bias.pro
.compile robomean.pro
.compile fi_peak.pro
.compile badpar.pro ; IY Buie
.compile trimrank.pro ; IY Buie
.compile goodpoly.pro
.compile read_cube.pro
.compile moment4.pro
.compile center_frames_WOLL2.pro
.compile create_header_woll_2.pro
.compile create_initial_data_woll2.pro
.compile peak_position.pro
.compile traectory_woll2.pro
.compile geometry_woll2.pro
.compile expand_traectory.pro
.compile interpol.pro
.compile create_repers_woll2.pro
.compile intersection_woll2.pro
.compile find_peaks.pro
.compile find_peaks_woll2.pro
.compile robust_poly_fit.pro
.compile rob_checkfit.pro
.compile robust_sigma.pro
.compile read_cube_woll2.pro
.compile med.pro
.compile polyfitw.pro
.compile poly.pro
.compile create_etalon_WOLL2.pro
.compile correction_image.pro
.compile dispersion_woll2.pro
.compile linerisation_woll2.pro
.compile create_flat_woll2.pro
.compile subsky_woll2.pro
.compile atm_absorbtion_woll2.pro
.compile   extraction_woll2.pro
.run  WOLL2_v1.pro
