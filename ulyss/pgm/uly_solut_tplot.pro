;+
; NAME:    
;                 ULY_SOLUT_TPLOT
;
; PURPOSE:
;                 Plot a result file of ULYSS
; 
; USAGE:
;                 uly_solut_tplot, <filename>[, PS=ps]
;                                 [, SKIPLINE=skipline][, NUMLINE=numline]
;                                 [, XAXIS=xaxis][, YAXIS=yaxis]
;                                 [, XERROR=xerror][, YERROR=yerror]
;                                 [, XRANGE=xrange][, YRANGE=yrange]
;                                 [, /XLOG][, /YLOG]
;                                 [, /CVG]
;                                 [, HISTOGRAM=nbins]
;                                 [, TITLE=title]
;                                 [, QUIET=quiet]
;                                 [, /OVERPLOT]
;
; ARGUMENT:
;   <filename>    Filename or solut structure.
;
; KEYWORDS:
;   PS            Set this keyword to generate Postscript output, otherwise
;                 the procedure makes a X window plot.
;
;   SKIPLINE      Number of lines to skip at the beginning of the .res file.
;                 By default, read from the first line.
;
;   NUMLINE       Number of lines to read from the .res file. By default
;                 read the entire file.
;
;   XAXIS YAXIS   Description of X & Y-Axes. 
;                 These keywords can be either an integer number or arrays
;                 of two integers describing the field number. 
;                 When a single value is given, 0 is for the row index,
;                 and fields in the file are numbered from 1. Each field usually
;                 consists in three columns (value, error and guess)
;                 or (weight, weight error and light-fraction).
;                 When a 2-integers array is given, they are the "coordinates"
;                 of a free parameter, as returned by ULY_CMP_SELECT.
;                 If any of XAXIS or YAXIS is not given, ULY_CMP_SELECT is 
;                 called to select the corresponding axis or axes.
;
;   /XERROR /YERROR 
;                 Instruct to respectively plot error bars for X-, Y-  or
;                 both axes.
;
;   XRANGE YRANGE
;                 Two elements arrays specifying the range for respectively
;                 the X and Y axes.
;
;   /XLOG /YLOG
;                 Plot X (resp. Y) axis in logarithmic scale.
;
;   /CVG          Connect with vectors the guesses with the solutions
;                 (convergence map).
;
;   HISTOGRAM     Make an histogram of the column X.
;                 If given the keywords XERROR, YERROR and CVG are
;                 ignored. The value of HISTOGRAM is the number of bins
;                 to use to build the histogram.
;
;   TITLE:        Title of the plot.
;
;   /QUIET        Set this keyword to eliminate information printed on
;                 the terminal.
;
;   /OVERPLOT     Overplots the result on an existing plot window.
;
; OUTPUT:
;
;   STD_X, STD_Y  Set these keywords to a variable which will receive
;                 the standard deviation of X,Y arrays
;
; DESCRIPTION:
;     This program reads a .res file produced by ULYSS, and produces a scatter
;     plot of one column versus another one or histogram of one column.
;     Error bars can be optionally drawn on any or both axes.
;
;     Vectors linking the guess to the solution can be shown with the keyword
;     /CVG. This modes allows to study the convergence (these graphics are
;     called "convergence maps"), for example to assess the sensitivity to
;     the guesses and give an idea of the convergence basin.
;
; EXAMPLES:
;   
; AUTHOR:
;                 The author in UNKNOWN (probably Mina Koleva)
;                 Modified by Philippe Prugniel, 2008/08/03
;
;-
; CATEGORY:       ULY
;------------------------------------------------------------------------------

function uly_solut_vector, solution, axis, i, j, MESSAGE=mess
; extract a column vector from the solution array

mess = ''

if n_elements(solution) eq 0 then begin
    mess = 'solution array is empty'
    return, 0
endif
if n_elements(i) eq 0 then i = 0 
if n_elements(j) eq 0 then j = n_elements(solution)-1

if n_elements(axis) eq 1 then begin           ;  specify a column by its number
    if axis eq 0 then begin                   ; row index
        return, {value:findgen(j-i+1) + i, $
                 error:fltarr(j-i+1), $
                 guess:findgen(j-i+1)+i, $
                 label:'index', unit:''}
    endif
    
    nax = n_elements(solution[i].losvd)
    n1 = -1
    if axis le nax then n2 = axis-1 else begin
        for n1=0,n_elements(solution[i].cmp)-1 do begin
            nax ++
            n2 = -1
            if nax eq axis then break     ; found weight
            for n2=0,n_elements(*solution[i].cmp[n1].para)-1 do begin
                nax++
                if nax eq axis then break  ; found a para in a cmp
            endfor
            if nax eq axis then break ;
        endfor
        
        if nax ne axis then begin
            mess = 'Invalid axis number (' +strtrim(string(axis),2)+')' 
            return, 0
        endif
    endelse

endif else begin
    n1 = axis[0]
    n2 = axis[1]
    if n1 ge n_elements(solution[i].cmp) then begin
        mess = 'Invalid axis, cmp number is too large (' + $
          strtrim(string(n1),2)+')' 
        return, 0
    endif
    if n1 ge 0 then begin
        if n2 ge n_elements(*solution[i].cmp[n1].para) then begin
            mess = 'Invalid axis, param number is too large (' + $
              strtrim(string(n2),2)+')' 
            return, 0
        endif
    endif
endelse

if n1 ne -1 then begin          ; one of the cmp
    if n2 ge 0 then begin      ; one of the param of this cmp
        label = (*solution[i].cmp[n1].para)[n2].name
        unit = (*solution[i].cmp[n1].para)[n2].unit
        x = dblarr(j-i+1, /NOZERO)
        e_x = dblarr(j-i+1, /NOZERO)
        g_x = dblarr(j-i+1, /NOZERO)
        for k=i,j do begin
            cmp = solution[k].cmp
            x[k-i]   = (*cmp[n1].para)[n2].value
            e_x[k-i] = (*cmp[n1].para)[n2].error
            g_x[k-i] = (*(*cmp[n1].para)[n2].guess)[0]
        endfor
        if (*solution[i].cmp[n1].para)[n2].dispf eq 'exp' then begin
            x = exp(x)
            e_x *= x
            g_x = exp(g_x)
        endif          
    endif else if n2 eq -1 then begin            ; the weight of this cmp  
        label = 'weight' 
        unit = ''
        x = solution[i:j].cmp[n1].weight
        e_x = solution[i:j].cmp[n1].e_weight
        g_x = 0
    endif else  if n2 eq -2 then begin            ; the light of this cmp  
        label = 'light' 
        unit = '%'
        tw = total (solution[i:j].cmp.l_weight, 1) /100
        if n_elements(solution[0].cmp) ne 1 then begin
           x = solution[i:j].cmp[n1].l_weight / tw
           e_x = solution[i:j].cmp[n1].e_weight /solution[i:j].cmp[n1].weight $
                 * solution[i:j].cmp[n1].l_weight / tw
           g_x = 0
        endif else begin ;if there is only one component
           x = replicate(100.,j-i)
           e_x = replicate(0.0,j-i)
           g_x = 0.
        endelse
    endif
endif else begin                 ; one of the LOSVD parameters
    x = solution[i:j].losvd[n2]
    e_x = solution[i:j].e_losvd[n2]
    g_x = 0
    label = (['cz', 'velocity dispersion', 'h3', 'h4', 'h5', 'h6'])[n2]
    unit = (['km/s', 'km/s', '', '', '', ''])[n2]
endelse

return, {value:x, error:e_x, guess:g_x, label:label, unit:unit}

end

; ----------------------------------------------------------------------------
pro uly_solut_tplot, filename, PS=ps,                    $
                     SKIPLINE=skipline, NUMLINE=numline, $
                     XAXIS=xaxis, YAXIS=yaxis,           $
                     XERROR=xerror, YERROR=yerror,       $
                     XRANGE=xrange, YRANGE=yrange,       $
                     XLOG=xlog, YLOG=ylog,               $
                     CVG=cvg,                            $
                     HISTOGRAM=nbins,                    $
		     OVERPLOT=overplot, 		 $
                     QUIET=quiet,                        $
                     TITLE=title,                        $
                     XTITLE=xtitle, YTITLE=ytitle,       $
		     PLOT_VAR=plot_var,                  $
		     STD_X=std_x,STD_Y=std_y,            $
		     MNT_X=mnt_x,MNT_Y=mnt_y,            $
                     _EXTRA=extra

if n_elements(filename) eq 0 then begin
    print, 'Usage:'
    print, 'uly_solut_tplot, filename, PS=<ps>,          $'
    print, '       SKIPLINE=<skipline>, NUMLINE=<numline>, $'
    print, '       XAXIS=<xaxis>, YAXIS=<yaxis>,           $'
    print, '       CVG=cvg,                                $'
    print, '       HISTOGRAM=nbins,                        $'
    print, '       QUIET=quiet,                            $'
    print, '       PLOT_VAR=plot_var'
endif

if not keyword_set(title) then begin
    if size(filename, /TYPE) eq 7 then title = filename[0] else title = ''
endif
if not keyword_set(plot_var) or size(plot_var, /type) ne 8 then begin
  if not keyword_set(psym) then psym=4
  plot_var = uly_plot_init(PSYM=psym, _EXTRA=extra)
endif

if size(filename, /TYPE) ne 8 then $
  solution = uly_solut_tread(filename) $    ; Read the .res ASCII file
else solution = filename                    ; Solut structure in input

if size(solution, /TYPE) ne 8 then $
  message, 'Could not read the input ASCII file'

; Select the axis to plot (we assume that all the records in the file
; are identical ...)
if n_elements(xaxis) eq 0 then begin
    print, 'Select XAXIS now...:'
    xaxis = uly_cmp_select(solution[0].cmp)
endif

if n_elements(nbins) eq 0 then begin
    if n_elements(yaxis) eq 0 then begin
        print, 'Select YAXIS now...:'
        yaxis = uly_cmp_select(solution[0].cmp)
    endif
endif

; Select the range of lines to extract
i = 0
j = (size(solution))[1] - 1
if (n_elements(skipline) ne 0) then i = skipline
if (n_elements(numline) ne 0) then j = min(i + numline - 1,(size(solution))[1] - 1)


; Extract the proper vectors to plot
x = uly_solut_vector(solution, xaxis, i, j, MESS=mess)
if mess ne '' then begin
    message, 'Could not extract X-axis .. '+mess, /INFO
    return
endif
err_x = mean(x.error,/DOUBLE)
mnt_x = moment(x.value, SDEV=std_x)
if not keyword_set(quiet) then $
  print, x.label,' mean: ',mnt_x[0], ' ['+x.unit+']', ' deviation:', std_x,' error:', err_x

if n_elements(nbins) eq 0 then begin
    y = uly_solut_vector(solution, yaxis, i, j, MESS=mess)
    if mess ne '' then begin
        message, 'Could not extract Y-axis .. '+mess, /INFO
        return
    endif
    if n_elements(y) ne n_elements(x) then begin
        message, 'X and Y axis have different length', /INFO
    endif
    err_y = mean(y.error,/DOUBLE)
    mnt_y = moment(y.value, SDEV=std_y)
    if not keyword_set(quiet) then $
      print, y.label,' mean: ',mnt_y[0], ' ['+y.unit+']', ' deviation:', std_y,' error:', err_y
endif else $
  y = {value:fltarr(j-i+1), error:fltarr(j-i+1), guess:fltarr(j-i+1), label:'', unit:''} 

if size(filename, /TYPE) ne 8 then heap_free, solution.cmp 

if keyword_set(xlog) then begin
    ind = where(x.value gt 0, c)
    if c eq 0 then begin
        message, 'Cannot plot in log, because X has no positive value', /INFO
        return
    endif
    val = x.value[ind]
    err = x.error[ind]
    gue = x.guess[ind]
    x = {value:val, error:err, guess:gue, label:x.label, unit:x.unit} 
endif

if keyword_set(ylog) then begin
    ind = where(y.value gt 0, c)
    if c eq 0 then begin
        message, 'Cannot plot in log, because Y has no positive value', /INFO
        return
    endif
    val = y.value[ind]
    err = y.error[ind]
    gue = y.guess[ind]
    y = {value:val, error:err, guess:gue, label:y.label, unit:y.unit} 
endif

if n_elements(xrange) ne 2 then xrange = [0,0]
if n_elements(yrange) ne 2 then yrange = [0,0]

if not keyword_set(quiet) then begin
    print, 'xrange= [', xrange[0], ',', xrange[1], ']'
    print, 'yrange= [', yrange[0], ',', yrange[1], ']'
endif

if keyword_set(ps) then begin
    set_plot,"ps"
    if n_elements(ps) ne 0 then psfilename=ps else psfilename = filename+'.ps'
    DEVICE, FILENAME=psfilename, /COLOR, /ENCAPSULATED,$
      XS=20, YS=20,/BOLD, FONT_SIZE=12  
endif else begin
;    window, XS=900, YS=450
;    device, DECOMPOSE=0
endelse

; define colors for the background and the forground of the plot
blue = FSC_Color('Blue')

if not keyword_set(xtitle) then begin
   xtitle = x.label
   if x.unit ne '' then xtitle += ', '+x.unit
   if keyword_set(xlog) then xtitle = 'log('+xtitle+')'   
endif

if not keyword_set(ytitle) then begin
   ytitle = y.label
   if y.unit ne '' then ytitle += ', '+y.unit
   if keyword_set(ylog) then ytitle = 'log('+ytitle+')'
endif

if n_elements(nbins) gt 0 then begin     ; Plot an histogram
    histoplot, x.value, NBINS=nbins
;   we could use more options, see: 
;   http://www.dfanning.com/graphics_tips/histoplot.html
endif else begin
;   The log plot mode of the IDL function PLOT makes bad axis when the
;   data range is less than a decade... in this case we make the ticks 
;   ourselves
    if keyword_set(xlog) then begin
        rg = minmax(xrange)
        if floor(alog10(rg[1]))-ceil(alog10(rg[0])) lt 1 then begin
            xtickv = tickslog(xrange, NMIN=4)
            xticks = n_elements(xtickv)-1
        endif
    endif
    if keyword_set(ylog) then begin
        rg = minmax(yrange)
        if floor(alog10(rg[1]))-ceil(alog10(rg[0])) lt 1 then begin
            ytickv = tickslog(yrange, NMIN=4)
            yticks = n_elements(ytickv)-1
        endif
    endif

    if not keyword_set(overplot) then                              $
       plot, x.value, y.value, /NODATA,                            $
             BACK=plot_var.bgcolor,                                $
             COLOR=plot_var.axiscolor,                             $
             TITLE=title, CHARSIZE=plot_var.charsize,              $
             XTITLE=xtitle, XSTYLE=3, XLOG=xlog,                   $
             XTICKS=xticks, XTICKV=xtickv, XRANGE=xrange,          $
             YTITLE=ytitle, YSTYLE=3, YLOG=ylog,                   $
             YTICKS=yticks, YTICKV=ytickv, YRANGE=yrange
    
    if keyword_set(xerror) or keyword_set(yerror) then begin
       if not keyword_set(xerror) then x.error = 0*x.value
       if not keyword_set(yerror) then y.error = 0*y.value
       oploterror, x.value, y.value, x.error, y.error, PSYM=4, SYMS=0.3, $
                   ERRCOLOR=blue, COLOR=plot_var.linecolor
    endif else $
       oplot, x.value, y.value, PSYM=plot_var.psym, SYMS=0.3, $
              COLOR=plot_var.linecolor, $
              LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick
    if keyword_set(cvg) then begin
       for k=0,n_elements(x.value)-1 do begin
          oplot, [x.value[k], x.guess[k]], [y.value[k], y.guess[k]], $
                 COLOR=plot_var.linecolor
       endfor
    endif
endelse

if keyword_set(ps) then begin
    Device,/CLOSE
    set_plot,"x"
endif

end
