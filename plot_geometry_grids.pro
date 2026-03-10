; Скрипт для визуализации сеток из файла geometry.fit
pro plot_geometry_grids

  ; Чтение данных из файла
  geom = readfits('geometry.fit')
  
  ; Анализ реальных размерностей
  dims = size(geom)
  print, 'Geometry dimensions: ', dims[1], dims[2], dims[3]
  
  ; 1000 точек × 4 координаты × 4 кадра
  n_points = dims[1]  ; 1000
  n_coords = dims[2]  ; 4 (X0, Y0, X1, Y1)
  n_frames = dims[3]  ; 4 кадра
  
  ; Визуализация для каждого кадра
  for frame = 0, n_frames - 1 do begin
    ; Извлекаем данные для текущего кадра
    frame_data = geom[*, *, frame]
    
    ; Разделяем координаты
    X0 = frame_data[*, 0]  ; исходные X координаты
    Y0 = frame_data[*, 1]  ; исходные Y координаты  
    X1 = frame_data[*, 2]  ; преобразованные X координаты
    Y1 = frame_data[*, 3]  ; преобразованные Y координаты

print,N_elements(X0)
print, 'X0 min: ', min(X0), ' X0 max: ', max(X0)
zero_count = total(X0 eq 0)
non_zero_count = total(X0 ne 0)
print, ' Zero values X0: ', zero_count, ' (', float(zero_count)/n_elements(X0)*100, '%)'
print, ' Non-zero values X0: ', non_zero_count, ' (', float(non_zero_count)/n_elements(X0)*100, '%)'
print,N_elements(X1)
print, 'X1 min: ', min(X1), ' X1 max: ', max(X1)
zero_count = total(X1 eq 0)
non_zero_count = total(X1 ne 0)
print, ' Zero values : ', zero_count, ' (', float(zero_count)/n_elements(X1)*100, '%)'
print, ' Non-zero values : ', non_zero_count, ' (', float(non_zero_count)/n_elements(X1)*100, '%)'

print,N_elements(Y0)
print, 'Y0 min: ', min(Y0), ' Y0 max: ', max(Y0)
zero_count = total(Y0 eq 0)
non_zero_count = total(Y0 ne 0)
print, ' Zero values : ', zero_count, ' (', float(zero_count)/n_elements(Y0)*100, '%)'
print, ' Non-zero values : ', non_zero_count, ' (', float(non_zero_count)/n_elements(Y0)*100, '%)'

print,N_elements(Y1)
print, 'Y1 min: ', min(Y1), ' Y1 max: ', max(Y1)
zero_count = total(Y1 eq 0)
non_zero_count = total(Y1 ne 0)
print, ' Zero values : ', zero_count, ' (', float(zero_count)/n_elements(Y1)*100, '%)'
print, ' Non-zero values : ', non_zero_count, ' (', float(non_zero_count)/n_elements(Y1)*100, '%)'

    ; Создаем окно для текущего кадра
    window, frame, xsize=800, ysize=600, $
           title='Geometry Grid - Frame ' + strtrim(string(frame), 2)
    
    ; Определяем диапазоны для осей
    all_x = [X0, X1]
    all_y = [Y0, Y1]
    x_min = min(all_x) - 20
    x_max = max(all_x) + 20
    y_min = min(all_y) - 20
    y_max = max(all_y) + 20
    
    ; Рисуем фон
    plot, [x_min, x_max], [y_min, y_max], /nodata, $
          xtitle='Dispersion Axis (X)', $
          ytitle='Slit Axis (Y)'
    
    ; Рисуем исходную сетку (X0, Y0) - красные крестики
    plots, X0, Y0, psym=1, symsize=1.2, color=10e6, thick=1.5
    
    ; Рисуем преобразованную сетку (X1, Y1) - зеленые квадраты
    plots, X1, Y1, psym=6, symsize=0.8, color=10e6, thick=1.5
    
    ; Добавляем соединительные линии между соответствующими точками
    for i = 0, n_points - 1 do begin
      if (i mod 1 eq 0) then begin ; рисуем каждую 20-ю точку для наглядности
        plots, [X0[i], X1[i]], [Y0[i], Y1[i]], color=9e6, thick=0.75, linestyle=0
      endif
    endfor
    
  endfor
  
  print, 'Geometry grids plotted for ', n_frames, ' frames'

end

; Диагностический скрипт для точного анализа структуры
pro analyze_geometry

  geom = readfits('geometry.fit')
  real_dims = size(geom)
  print, 'Real dimensions: ', real_dims[1], 'x', real_dims[2], 'x', real_dims[3]
  print, 'Total elements: ', n_elements(geom)
  
end

; Запуск анализа
analyze_geometry
print, '---'
; Запуск визуализации
plot_geometry_grids
END
