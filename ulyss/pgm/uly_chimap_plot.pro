;+
; NAME:
;              ULY_CHIMAP_PLOT
;
; PURPOSE:
;              Plot a chi square map
;
; USAGE:
;              uly_chimap_plot, <chimap>
;                       [, /PS] or [PS='filename.ps']  
;                       [, /NOSCALEBAR]
;                       [, XRANGE=xrange][, YRANGE=yrange]
;                       [, /XLOG][, /YLOG]
;              		[, PLOT_VAR=<plot_var>]
;			[, /QUIET]
;			[, ...extra keywords]
;
; ARGUMENTS:
;    <chimap>: An input map structure, or FITS file containing the map.
;              Chi2 maps are generated with the function CHIMAP.
;              The function CHIMAP_READ may be used to read such a file
;              and return a map structure, or the name of the file can 
;              directly be passed here.
; KEYWORDS:
;    /PS or PS='filename.ps'
;   	       Makes a postcript hardcopy of the plot. When using /PS the
;              filename is called 'chimap.ps'. To use a different name,
;              specify it with PS='other_name.ps'.
;
;    /OVERPLOT 
;              Overplots on an existing graphics.
;    
;    /NOSCALEBAR
;              When a halftone map is drawn, setting this keyword suppress
;              the scale bar which is otherwise plotted on the right side.
;
;    /QUIET  Suppress some message printed on the screen.
;
;    XRANGE YRANGE
;              Two-elements arrays specifying the range for respectively
;              the X and Y axes.
; 
;    CHIRANGE
;              Two-elements arrays specifying the range of chi2 values
;              corresponding to the extreme countours and/or colours.
;
;    /XLOG /YLOG
;              Plot X (resp. Y) axis in logarithmic scale.
;
;    [ ...Extra keywords]
;	       All extra keywords are passed to uly_plot_init, if PLOT_VAR
;	       is not defined.
;
; DESCRIPTION:
;    Draw a chi2 map with optionally halftone background and iso-contours
;    overlay. The node of the local and absolute minima are figured by
;    overplotted symbols (note that these minima are the position of
;    the pixels of the map and do not generally coincide with the
;    result of the minimization since the map has a limited resolution).
;
; HISTORY:
;              Antoine Bouchard, 2008/02/08, created
;-
; CATEGORY:    ULY

pro uly_chimap_plot, chimap, PS=ps, OVERPLOT=overplot,   $
                     NOSCALEBAR=noscalebar,              $
                     XRANGE=xrange, YRANGE=yrange,       $
                     CHIRANGE=chirange,                  $
                     XLOG=xlog, YLOG=ylog,               $
                     QUIET=quiet, PLOT_VAR=plot_var, _EXTRA=extra

  if n_elements(chimap) eq 0 then begin ; Check for valid input

     print, 'ULY_CHIMAP_PLOT: No input data supplied.'
     print, 'Usage:   uly_chimap_plot, <chimap> [,PS]'
     return
  endif

  if not keyword_set(plot_var) or size(plot_var, /TYPE) ne 8 then begin
     if not keyword_set(halftonect ) then $
        plot_var = uly_plot_init(HALFTONECT=3, _EXTRA=extra) else $
           plot_var = uly_plot_init(_EXTRA=extra)
  endif

  if size(chimap, /TYPE) eq 7 then begin
     file = chimap
     chimap = uly_chimap_read(chimap)
     readfile = 1
  endif else readfile = 0

;if tag_names(chimap) ne tag_names(uly_chimap_init(1,1)

; Reads the map, if needed.

  if keyword_set(ps) then begin ; Set default parameters for postscript output
     set_plot, 'ps'
     if n_elements(ps) ne 0 then filename=ps else filename = 'chimap.ps'
     device, FILENAME=filename, /COLOR
     if not keyword_set(quiet) then print, 'Writing file: ', filename
  endif else if !d.name ne 'PS' then begin ; Initialize the plotting window
     device, WINDOW_STATE=wstate, DECOMPOSED=0
     if total(wstate) lt 1. then window, 0, XSIZE=plot_var.xsize, YSIZE=plot_var.ysize
  endif


  if n_elements(xrange) eq 2 then begin
     x1 = max([(where(chimap.xnode ge xrange[0]))[0]-1, 0])
     x2 = (where(chimap.xnode ge xrange[1], cx2))[0]
     if cx2 eq 0 then x2 = n_elements(chimap.xnode)-1
  endif else begin
     x1 = 0
     x2 = n_elements(chimap.xnode)-1
  endelse

  if n_elements(yrange) eq 2 then begin
     y1 = max([(where(chimap.ynode ge yrange[0]))[0]-1, 0])
     y2 = (where(chimap.ynode ge yrange[1], cy2))[0]
     if cy2 eq 0 then y2 = n_elements(chimap.ynode)-1
  endif else begin
     y1 = 0
     y2 = n_elements(chimap.ynode)-1
  endelse

  map = chimap.map[x1:x2,y1:y2]

  posit = where(map gt 0, cnt)
  if cnt eq 0 then return
  mn    = alog(min(map[posit], location, /NAN, MAX=mx))
  mx    = alog(mx)
  location = posit[location]
  if mx eq mn then mx += 0.005  ; map is flat, or only 1 px
  mn_i  = array_indices(map, location) + [x1, y1]

  if n_elements(chirange) ge 1 then mn = alog(chirange[0])
  if n_elements(chirange) ge 2 then mx = alog(chirange[1])

  nlevs = exp(((mx-mn)*findgen(22)/21. + mn)[0:19])
  clevs = exp(((mx-mn)*findgen(62)/61. + mn)[0:60])

  xtitle = chimap.xtitle
  if chimap.xunit ne '' then xtitle += ', ' + strtrim(chimap.xunit, 2)
  if keyword_set(xlog) then begin
     zero = where(chimap.xnode le 0, c)
     if c gt 0 then $
        message, 'Cannot use XLOG because some positions are negative or null'
     xnode = chimap.xnode
     if n_elements(xrange) ne 2 then rg = xnode([x1,x2]) else rg = xrange
     rgo = minmax(rg)
     if floor(alog10(rgo[1]))-ceil(alog10(rgo[0])) lt 1 then begin
        xtickv = tickslog(rg, NMIN=4)
        xticks = n_elements(xtickv)-1
     endif
  endif else xnode = chimap.xnode

  ytitle = chimap.ytitle
  if chimap.yunit ne '' then ytitle += ', ' + strtrim(chimap.yunit, 2)
  if keyword_set(ylog) then begin
     zero = where(chimap.ynode le 0, c)
     if c gt 0 then $
        message, 'Cannot use YLOG because some positions are negative or null'
     ynode = chimap.ynode
     if n_elements(yrange) ne 2 then rg = ynode([y1,y2]) else rg = yrange
     rgo = minmax(rg)
     if floor(alog10(rgo[1]))-ceil(alog10(rgo[0])) lt 1 then begin
        ytickv = tickslog(rg, NMIN=4)
        yticks = n_elements(ytickv)-1
     endif
  endif else ynode = chimap.ynode

  lm1_i = where(map lt shift(map,1,0)  and $
                map lt shift(map,-1,0) and $
                map lt shift(map,0,1)  and $
                map lt shift(map,0,-1) and $
                map lt shift(map,1,1)  and $
                map lt shift(map,1,-1) and $
                map lt shift(map,-1,1) and $
                map lt shift(map,-1,-1), c)

  if c gt 0 then begin
     lm2_i = array_indices(map, lm1_i) + rebin([x1,y1], 2, n_elements(lm1_i))
     if not keyword_set(quiet) then begin
        for k = 0, n_elements(lm1_i)-1 do begin
           print, 'Local minimum: x=', xnode(lm2_i[0,k]), $
                  ' y=', ynode(lm2_i[1,k]), ' chi2=', chimap.map[lm2_i[0,k],lm2_i[1,k]]
        endfor
     endif
  endif

                                ; device settings
  if !d.name eq 'PS' then white = 255 else white = !p.color
  x_win_size = !d.x_size/!d.x_ch_size ; find the end of the plot window    

  if not keyword_set(noscalebar) then begin
     if keyword_set(ps) then xmargin_main = [10,18] $
     else xmargin_main = [10,23]
     if keyword_set(ps) then xmargin_side = [68,8] $
     else xmargin_side = [x_win_size-14,8]
  endif

  if plot_var.showhalftone eq 1 then nodata=0 else nodata=1

  if not keyword_set(overplot) then $
     contour, chimap.map, xnode, ynode, NODATA=nodata, /FILL, LEVELS=clevs,      $
              BACK=plot_var.bgcolor, COLOR=plot_var.axiscolor,                   $
              XSTYLE=1, XTITLE=xtitle, XRANGE=xrange, XLOG=xlog,                 $
              XTICKS=xticks, XTICKV=xtickv, XMARGIN=xmargin_main,                $
              YSTYLE=1, YTITLE=ytitle, YRANGE=yrange, YLOG=ylog,                 $
              YTICKS=yticks, YTICKV=ytickv
  
;if plot_var.showhalftone eq 1 then $
; contour, chimap.map, xnode, ynode, levels=clevs, /fill, color=FSC_Color(plot_var.halftonecolor), /overplot

; Draw the contours
  if plot_var.showcontour eq 1 then begin
     contour, chimap.map, xnode, ynode, LEVELS=nlevs, COLOR=plot_var.linecolor, /OVERPLOT
  endif

  if plot_var.showminima eq 1 then begin
;   plot the local minima
     if n_elements(lm2_i) gt 0 then begin
        oplot, [xnode[lm2_i[0,*]]], [ynode[lm2_i[1,*]]], COLOR=FSC_Color(plot_var.sigmacolor), PSYM=4
     endif

;   plot the absolute minimum 
     oplot, [xnode[mn_i[0]]], [ynode[mn_i[1]]], COLOR=FSC_Color(plot_var.sigmacolor), PSYM=1, THICK=plot_var.linethick, SYMSIZE=2
     
  endif
  pchimap = !p & xchimap = !x & ychimap = !y
  if plot_var.showhalftone eq 1 and not keyword_set(noscalebar) then begin ; Draw the graycsale bar
     grayscale = fltarr(2,1000)
     gray_y    = fltarr(1000)
     for j=0, 999 do begin
        grayscale[0:1,j] = exp((j/999.)*(mx-mn)+mn)
        gray_y[j] = exp((j/999.)*(mx-mn)+mn)
     endfor
     
     contour, grayscale, [0,1], gray_y, /FILL, /CLOSED, /NOERASE, $
              COLOR=plot_var.axiscolor, BACK=plot_var.bgcolor,           $
              LEVELS=clevs, $
              XTICKS=1, XMARGIN=xmargin_side, XMINOR=1, $
              YSTYLE=1, XTICKNAME=replicate(' ', 30), YTICKNAME=replicate(' ', 30)
     
     axis, YAXIS=1, YSTYLE=1, YTITLE='!7v!dm!n!3!u2!n', $
           COLOR=plot_var.axiscolor
  endif

  !p = pchimap & !x = xchimap & !y = ychimap ; Restore the plotting
                                ; parameters so that the
                                ; oplot command can be used
                                ; on the main plot (not the
                                ; colorscale)

  if keyword_set(ps) then begin ; Close postscript file output
     device, /CLOSE
     set_plot, 'x'
  endif 

  if readfile eq 1 then chimap = file
  return

end

