pro visualize_spectra, filename=filename, save_fits=save_fits

  ; Загружаем данные и заголовок
  if ~keyword_set(filename) then filename = 'spectra.fit'
  if ~file_test(filename) then begin
	  print, 'Ошибка: файл ' + filename + ' не найден!'
	  return
  endif

  sp = readfits(filename)
  h = headfits(filename)
  
  ; Определяем ось длин волн
  Nx = sxpar(h, 'NAXIS1')
  lam = indgen(Nx) * sxpar(h, 'CDELT1') + sxpar(h, 'CRVAL1')
  
  ; Размеры данных
  dims = size(sp, /dimensions)
  n_dims = n_elements(dims)

  if n_dims eq 3 then begin
	  n_waves = dims[0]
	  n_pol = dims[1]
	  n_obj = dims[2]
	  n_exp = 1
	  has_exposures = 0
  endif else if n_dims eq 4 then begin
	  n_waves = dims[0]
      n_pol = dims[1]    ; 4 вектора поляризации
      n_exp = dims[2]    ; количество экспозиций
      n_obj = dims[3]    ; 3 объекта
      has_exposures = 1
  endif else begin
	  print, 'Unknown data dimensions'
	  return
  endelse
   
 ; print, 'Spectra length: ', n_waves, '\n Number of exposures: ', n_exp, '\n Number of objects: ', n_obj


  ; Цвета для векторов поляризации
  colors = ['red', 'blue', 'green', 'magenta']
  
  ; Для каждого объекта создаем отдельное окно
  for obj_idx = 0, n_obj-1 do begin
    ; Создаем новое окно для каждого объекта
    window, obj_idx, xsize=1200, ysize=1000, title='Object ' + strtrim(obj_idx+1, 2)
    
    ; Массив для хранения оверплотов
    if has_exposures then begin
		oplots = objarr(n_exp * n_pol)
	endif else begin
		oplots = objarr(n_pol)
	endelse
    plot_counter = 0

	if has_exposures then begin    
    ; Для каждой экспозиции
    for exp_idx = 0, n_exp-1 do begin
      ; Для каждого вектора поляризации
      for pol_idx = 0, n_pol-1 do begin
        ; Создаем оверплот для конкретной поляризации и экспозиции
        sp_plot = Obj_New('cgOverPlot', lam, sp(*, pol_idx, exp_idx, obj_idx), $
                          color=colors[pol_idx], $
                          linestyle=exp_idx mod 4) ; разные стили линий для разных экспозиций
                          
        
        ; Сохраняем указатель на оверплот
        oplots[plot_counter] = sp_plot
        plot_counter = plot_counter + 1
      endfor
    endfor
    endif else begin
		for pol_idx = 0, n_pol-1 do begin
			sp_plot = Obj_New('cgOverPlot', lam, sp(*, pol_idx, obj_idx), $
                          color=colors[pol_idx], $
                          linestyle=0, thick=2)
		    oplots[plot_counter] = sp_plot
			plot_counter = plot_counter + 1
		endfor
	  endelse
    
    ; Создаем первый график для масштабирования осей
    if has_exposures then begin
	; Используем данные первой поляризации первой экспозиции
    cgplot, lam, sp(*, 0, 0, obj_idx), $
            xstyle=1, ystyle=1, $          ; включаем все стили осей
            xtitle='Wavelength', $
            ytitle='Flux', $
            title='Object ' + strtrim(obj_idx+1, 2), $
            charsize=1.2, $
            oplots=oplots
    endif else begin
		cgplot, lam, sp(*, 0, obj_idx), $
              xstyle=1, ystyle=1, $
              xtitle='Wavelength', $
              ytitle='Flux', $
              title='Object ' + strtrim(obj_idx+1, 2), $
              charsize=1.2, $
              oplots=oplots
	endelse

    ; Очищаем оверплоты
    for i = 0, plot_counter-1 do begin
      if (obj_valid(oplots[i])) then obj_destroy, oplots[i]
    endfor
    
    ; Добавляем легенду (опционально)
 ;   cglegend, position=[0.15, 0.95], chars=0.8, /norm, $
 ;            text=['Pol 1', 'Pol 2', 'Pol 3', 'Pol 4'], $
 ;            colors=colors
  endfor
; Сохраняем отдельные FITS файлы для каждой экспозиции каждого объекта
    for obj_idx = 0, n_obj-1 do begin
	  if has_exposures then begin
	  for exp_idx = 0, n_exp-1 do begin
		  filename_out = 'object_' + strtrim(obj_idx+1, 2) + '_exp_' + strtrim(exp_idx+1, 2) + '.fits'
		 ; Извлекаем данные для всех 4 поляризаций этой экспозиции и объекта
		 ; data_to_save = reform(sp(*, *, exp_idx, obj_idx), n_waves, n_pol)
		  data_to_save =sp(*, *, exp_idx, obj_idx)
		 ; header = string('SIMPLE  =                    T')
		 ; header = [header, 'BITPIX  =                  -64']  ; float64
		 ; header = [header, 'NAXIS   =                    2']
		 ; header = [header, 'NAXIS1  =            ' + strtrim(n_waves, 2)]
		 ; header = [header, 'NAXIS2  =             ' + strtrim(n_pol, 2)]
		 ; header = [header, 'EXTEND  =                    T']
		 ; header = [header, 'OBJECT  = ' + strtrim(obj_idx+1, 2)]
		 ; header = [header, 'EXPOSURE= ' + strtrim(exp_idx+1, 2)]
		 ; header = [header, 'CDELT1  = ' + strtrim(sxpar(h, 'CDELT1'))]
		 ; header = [header, 'CRVAL1  = ' + strtrim(sxpar(h, 'CRVAL1'))]
		 ; header = [header, 'CRPIX1  =                    1']
		 ; header = [header, 'CTYPE1  = ''Wavelength''']

		  writefits, filename_out, data_to_save, h

		  print, 'Сохранен файл: ' + filename_out
		endfor
	  endif else begin
		  filename_out = 'object_' + strtrim(obj_idx+1, 2) + '_averaged.fits'
		  data_to_save = sp(*, *, obj_idx)
		  writefits, filename_out, data_to_save,h
		  print, 'Saved: ' + filename_out
	 endelse
	endfor
end
