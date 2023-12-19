; docformat = 'rst'
;+
; This is an example program to demonstrate how to create a 3D scatter plot on a map
; with Coyote Graphics routines. A copy of cgGrid_Map downloaded after 12 Feb 2013
; is required.
;
; :Categories:
;    Graphics
;
; :Examples:
;    Save the program as "scatter3d_on_map.pro" and run it like this::
;       IDL> .RUN scatter3d_on_map
;
; :Author:
;    FANNING SOFTWARE CONSULTING::
;       David W. Fanning
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: david@idlcoyote.com
;       Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 11 February 2013 by David W. Fanning.
;
; :Copyright:
;     Copyright (c) 2013, Fanning Software Consulting, Inc.
;-
PRO Scatter3D_On_Map

   ; Set up variables for the plot. Normally, these values would be
   ; passed into the program as positional and keyword parameters.
    seed = 1L
    x = RANDOMU(seed, 32)
    y = RANDOMU(seed, 32)
    z = EXP(-3 * ((x - 0.5)^2 + (y - 0.5)^2))

    ; Load a color table and create colors for the scatterplot.
    cgLoadCT, 33
    zcolors = BytScl(z)

   ; Scale the data locations into longitude and latitude coordinate space.
   x = cgScaleVector(x, -180, 180)
   y = cgScaleVector(y, -90, 90)

   ; Set the 3D coordinate space with axes.
   cgSurf, DIST(5), /NODATA, /SAVE, XRANGE=[-180,180], $
      YRANGE=[-90,90], ZRANGE=[0, 1], XSTYLE=1, $
      YSTYLE=5, ZSTYLE=5, CHARSIZE=2.0, $
      POSITION=[0.1, 0.1, 0.95, 0.95, 0.0, 1.0]

   ; Save the current plotting system variables (since MAP_SET will change them).
   bangp = !P
   bangx = !X
   bangy = !Y
   bangz = !Z

   ; Draw a map projection in the 3D coordinate space.
   cgMap_Set, /Cylindrical, /T3D, /NoErase, Position=[0.1, 0.1, 0.95, 0.95]
   cgMap_Continents, /Fill, /T3D, Color='sienna'
   cgMap_Grid, /T3D, Color='charcoal'

   ; Restore the system variables.
   !P = bangp
   !X = bangx
   !Y = bangy
   !Z = bangz

   ; Set the 3D coordinate space with axes. Need to repair axis damange here.
   cgSurf, DIST(5), /NODATA, /SAVE, XRANGE=[-180,180], $
      YRANGE=[-90,90], ZRANGE=[0, 1], XSTYLE=1, $
      YSTYLE=1, ZSTYLE=1, CHARSIZE=2.0, $
      POSITION=[0.1, 0.1, 0.95, 0.95, 0.0, 1.0], /NOERASE
   cgAXIS, XAXIS=1, /T3D, CHARSIZE=2.0
   cgAXIS, YAXIS=1, /T3D, CHARSIZE=2.0

    ; Plot the random points in 3D space with a filled circle shape.
    phi = Findgen(32) * (!PI * 2 / 32.)
    phi = [ phi, phi(0) ]
    cgPlotS, x, y, z, PSYM=16, COLOR=zcolors, SYMSIZE=2.5, /T3D

    ; Connect the data points to the XY plane of the plot.
    FOR j=0,31 DO cgPlotS, [x(j), x(j)], [y(j), y(j)], [0, z(j)], $
       COLOR=zcolors(j), /T3D


END ;*****************************************************************

; This main program shows how to call the program and produce
; various types of output.

  ; Display the plot in a graphics window.
  Scatter3D_On_Map

  ; Display the plot in a resizeable graphics window.
  cgWindow, 'Scatter3D_On_Map', WTitle='3D Scatter Plot on a Map in a Resizeable Graphics Window'

  ; Create a PostScript file.
  cgPS_Open, 'scatter3d_on_map.ps'
  Scatter3D_On_Map
  cgPS_Close

  ; Create a PNG file with a width of 600 pixels.
  cgPS2Raster, 'scatter3d_on_map.ps', /PNG, Width=600

END