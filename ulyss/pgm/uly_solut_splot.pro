;+
; NAME:
;               ULY_SOLUT_SPLOT
;
; PURPOSE:
;               Plot the result of a fit done with ULY_FIT
;
; USAGE:
;	        uly_solut_splot, <solution>
;                       [, WAVERANGE=<waverange>]
;                       [, SMOOTH_RESID=<smooth_resid>]
;                       [, SMOOTH_FLUX=<smooth_flux>]
;                       [, RESI_YR=<resi_yr>]
;                       [, FLUX_YR=<flux_yr>]
;                       [, SCAN=<scan>]
;                       [, PLOT_VAR=<plot_var>]
;                       [, ...extra keywords]
;
; DESCRIPTION:
;               This task plots the result of a fit performed by ULY_FIT. It
;               is also called by ULYSS if the /PLOT keyword is set. The
;               plot is either on an xwindow or in a postscript file. 
;
; ARGUMENTS:
;   <solution>	Input.
;               The <solution> structure is part of the normal output from
;               the ULY_FIT task. Alternatively, if ULY_FIT was used to
;               generate a file containing the output, then the file name
;               can also be given here. 
;               A <solution> structure is also outputed by the
;               ULY_SOLUT_SREAD task. 
;
; KEYWORDS:
;   WAVERANGE:
;               Array of one or two values specifying the wavelength range
;               to display (values in Angstrom). By default the full range is
;               displayed. If a single value is specified, it is taken as
;               the beginning of the range and the interval extends to
;               the end of the available wavelength range.
;
;   SMOOTH_RESID:
;               Velocity dispersion used to smooth the residuals (convolve),
;               in km/s. No smoothing by default.
;
;   SMOOTH_FLUX:
;               Velocity dispersion used to smooth the flux (convolve),
;               in km/s. No smoothing by default.
;
;   RESI_YR:
;               Scalar, or two-elements array. Define the Y range for the
;               residuals panel. If a scalar or one-element array is 
;               specified, it is the semi-amplitude of the range, centered
;               on zero.
;
;   FLUX_YR:
;               Define the flux range for the flux panel.
;
;   SCAN:       Choose which scan to plot in case of 2D file, by
;               default plot the first dim.
;
;   /XLOG /YLOG
;               Plot X (resp. Y) axis in logarithmic scale.
;
;   TITLE:
;               Title of the plot.
;   PLOT_VAR:
;               An input structure containing all plot variables.
;               This structure is returned by ULY_PLOT_INIT.
;
;   ...Extra keywords:
;               If PLOT_VAR is undefined, all extra keywords are
;               passed to ULY_PLOT_INIT.
; 
; EXAMPLE:
;      The following example shows how to prepare a publication quality
;      figure in postscript format (this example produces an acceptable
;      result for a single-column graphic in a journal, eg. MNRAS).
;
;      !P.FONT = 0
;      !P.THICK = 1.5
;      !X.THICK = 1.5
;      !Y.THICK = 1.5
;      set_plot,'ps'
;      device, FILENAME='fig.ps', /COLOR, FONT_SIZE=6, XSIZE=9.5, YSIZE=5, $
;              /TIMES, /ENCAPSUL
;      uly_solut_splot, 'output'
;      device, /CLOSE
;      set_plot, 'x'
;
;      The device command allows to set the size, position, font ...
;      of the plot
;
; AUTHOR:
;               Antoine Bouchard
; HISTORY:
;               2008/02/05
;-
; CATEGORY:     ULY
;------------------------------------------------------------------------------

pro uly_solut_splot, solution, $
                     WAVERANGE=waverange, RESI_YR=resi_yr, FLUX_YR=flux_yr, $
                     SMOOTH_RESID=smooth_resid, $
                     SMOOTH_FLUX=smooth_flux, $
                     SCAN=scan, $
                     XLOG=xlog, YLOG=ylog,         $
                     TITLE=title,                  $
                     PLOT_VAR=plot_var, _EXTRA=extra

; Check to see if a filename is given
if size(solution,/TYPE) eq size('str',/TYPE) then $
 solut = uly_solut_sread(solution) $
else $
 solut = solution

; Sanity Check
if size(solut, /TYPE) ne 8 then begin
   message, 'No input' ,/INFO
   return
endif
if n_elements(solut) ne 1 then begin
   message, /INFO, 'Cannot plot, because <solution> is an array'
   return
endif

; Initialise plotting variables
if not keyword_set(plot_var) or size(plot_var, /TYPE) ne 8 then begin
  if not keyword_set(linecolor) then linecolor='Black'
  plot_var = uly_plot_init(LINECOLOR=linecolor, _EXTRA=extra)
endif
if n_elements(title) eq 0 then title = solut.title

; choose which scan to plot (if 2D file)
if n_elements(scan) ne 0 then nscan = scan  $
else nscan = 0

if not ptr_valid(solut.data) then  begin
   message, /INFO, 'No data to plot (spectrum is empty)'
   return
endif
if (size(solut.bestfit,/N_DIM)) eq 0 then  begin
   message, /INFO, 'No data to plot (spectrum is empty)'
   return
endif
if size(solut.bestfit,/N_DIM) gt 1 then begin
   if nscan ge (size(solut.bestfit,/DIM))[1] then begin
       message, 'SCAN is out of range' ,/info
       return
   endif
endif else if nscan ne 0 then begin
   message, 'SCAN is out of range' ,/info
   return
endif


; Extracting usefull info
crv       = solut.start
cdv       = solut.step
npix      = (size(solut.bestfit,/DIM))[0]
wave      = findgen(npix)*cdv+crv
if solut.sampling eq 1 then wave = exp(wave)  ; input is in log(wave)
range     = [wave[0], wave[n_elements(wave)-1]]

if n_elements(waverange) gt 0 then begin   ; must extract a range
   if n_elements(waverange) eq 1 then $
     waverange = [waverange,wave[n_elements(wave)-1]]
   if waverange[0] lt range[0] then waverange[0] = range[0]
   if waverange[1] gt range[1] then waverange[1] = range[1]
   range = waverange
   p0 = min(where(wave ge range[0]))
   p1 = max(where(wave le range[1]))
   wave = wave[p0:p1]

endif else begin
   p0 = 0l
   p1 = npix-1
endelse

data      = dblarr(5, p1-p0+1)
data[0,*] = solut.bestfit[p0:p1,nscan]
if (size(*solut.data,/DIM))[0] ge p1 then begin
   data[1,*] = (*solut.data)[p0:p1,nscan]
   factor    = abs(1./median(data[1,*]))
endif else begin
   factor = 1
   message, 'No observed spectrum to plot', /INFO
endelse
if (size(solut.mulcont,/DIM))[0] ge p1 then $
 data[2,*] = solut.mulcont[p0:p1,nscan]
data[3,*] = solut.addcont[p0:p1,nscan]
data[4,*] = (*solut.err)[p0:p1,nscan]


; determining values to ignore when setting plotting ranges
goodpix   = where(solut.mask[p0:p1,nscan] eq 1)
badpix    = where(solut.mask[p0:p1,nscan] eq 0)

spec      = reform(factor*data[1,*]/data[2,*])
model     = reform(factor*data[0,*]/data[2,*])
resid     = spec - model

if n_elements(smooth_resid) eq 1 then begin
   chi = resid / data[4,*]
   velscale = cdv * 299792.458d
   uly_slit, 0d, smooth_resid, 0d, 0d, velscale, kernel
   status = check_math(MASK=32)
   resid = convol(chi,kernel,/EDGE_TRUNCATE) * data[4,*]
   data[4,*] = sqrt(convol(reform(data[4,*])^2,kernel^2,/EDGE_TRUNCATE))
   if status eq 0 then status = check_math(MASK=32)
endif
if n_elements(smooth_flux) eq 1 then begin
   uly_slit, 0d, smooth_flux, 0d, 0d, cdv * 299792.458d, kernel
   spec = convol(spec, kernel, /EDGE_TRUNCATE)
   model = convol(model, kernel, /EDGE_TRUNCATE)
   status = check_math(MASK=32)
endif

if array_equal(data[2,*], 0.*data[2,*]+1.) then begin
  if array_equal(data[3,*], 0.*data[3,*]) then $
    polyn = !VALUES.F_NAN*lindgen(p1-p0+1) else polyn = reform(data[3,*])
endif else begin
  polyn = reform(data[2,*])
endelse

badresid         = resid*!VALUES.F_NAN
goodresid        = resid*!VALUES.F_NAN
badspec          = spec*!VALUES.F_NAN
goodspec         = spec*!VALUES.F_NAN

goodspec[goodpix]  = spec[goodpix]

if n_elements(badpix) gt 1 then begin
   badpix = uly_plot_padpixelrange(badpix)
   if n_elements(smooth_flux) eq 1 then $
      badpix = uly_plot_padpixelrange(badpix)
   badspec[badpix]    = spec[badpix]

   residbadpix = badpix
   if n_elements(smooth_resid) eq 1 then $
;; I think this 'for ... kernel/10' is useless MK 2009/07/28
;     for i=0, n_elements(kernel)/10 do $
   ; When using the smooth_resid option, the residbadpix gets padded with a
   ; representative number of pixels (for simplicity this is the size of
   ; the kernel/10)
     residbadpix = uly_plot_padpixelrange(residbadpix)
      
   badresid[residbadpix]   = resid[residbadpix]
   residgoodpix = uly_plot_padpixelrange(where(finite(badresid, /nan)))
   goodresid[residgoodpix] = resid[residgoodpix]
endif else begin
   goodresid = resid
   residgoodpix = goodpix
endelse



; define colors for the background and the forground of the plot
if !d.name ne 'PS' then device, DECOMPOSED=0

; Setting the plotting ranges
if n_elements(flux_yr) ne 0 then begin
    yr_flux = flux_yr 
    ys_flux = 1
endif else begin
;    yr_flux = [min([spec[goodpix], polyn[goodpix], 0.]), max([spec[goodpix], polyn[goodpix], 1.5])]
    sorted = sort(spec[goodpix])
    n = n_elements(goodpix)
    yr_flux = (spec[goodpix])[sorted[[ceil(0.002*n), floor(0.998*n)]]]
endelse

if n_elements(resi_yr) ne 0 then begin
    yr_resi = resi_yr 
    if n_elements(resi_yr) eq 1 then yr_resi=[-abs(resi_yr), abs(resi_yr)]
    ys_resi = 1         ; the axis limit will be exactly as specified by the user
endif else begin
;    yr_resi = [min(goodresid, MAX=mx, /NAN), mx]
    sorted = sort(goodresid[residgoodpix])
    n = n_elements(residgoodpix)
    yr_resi = (goodresid[residgoodpix])[sorted[[ceil(0.003*n), floor(0.997*n)]]]
endelse


!p.multi=[0,1,2]

; This is the spectra and model plot
plot, wave, spec,                                                           $
  XRANGE=range, XLOG=xlog, XSTYLE=1, XTICKNAME=replicate(' ', 30),          $
  YRANGE=yr_flux, YLOG=ylog, YSTYLE=ys_flux, YMARGIN=[0., !Y.MARGIN[1]],    $
  TITLE=title, CHARSIZE=plot_var.charsize,                                  $
  BACK=plot_var.bgcolor, COLOR=plot_var.axiscolor,                          $
  /NODATA

; Compute the position of the Y-title with respect to the edge of the plotting
; region (axis+title), in order to aligne the title of the 2 panels
; Change the definition of !X.MARGIN[0] to modify the distance of the
; title to the axis ...
ytpos = float(!D.y_ch_size)/!D.x_size * 2 ; leave a 1-char margin
; Set ytpos = float(!D.y_ch_size)/!D.x_size to have the top of chars touch the edge
xyouts, ytpos, mean(!y.window), 'Relative flux', CHARSIZE=plot_var.charsize, /NORMAL, ALIGN=0.5, ORIENT=90, COLOR=plot_var.axiscolor

oplot, wave, goodspec, COLOR=plot_var.linecolor, $
  PSYM=plot_var.psym, LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick
oplot, wave, polyn, COLOR=FSC_Color(plot_var.polycolor), $
  PSYM=plot_var.polypsym, LINESTYLE=plot_var.polystyle, THICK=plot_var.polythick
oplot, wave, model, COLOR=FSC_Color(plot_var.modelcolor), $
  PSYM=plot_var.modelpsym, LINESTYLE=plot_var.modelstyle, THICK=plot_var.modelthick
if n_elements(badpix) gt 1 then $
   oplot, wave, badspec, $
          COLOR=FSC_Color(plot_var.badcolor), PSYM=plot_var.psym, $
          LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick

ymarg = [!Y.MARGIN[0],1.]
xtitle = 'Wavelength, '+'!3' + STRING(197B) + '!X' 

; Zoom on residual plot
plot, wave, resid, /NODATA, $
  BACK=plot_var.bgcolor, COLOR=plot_var.axiscolor,                       $
  XRANGE=range, XLOG=xlog, XSTYLE=1, XTITLE=xtitle,                      $
  YRANGE=yr_resi, YLOG=ylog, YMARGIN=ymarg,                              $
  YSTYLE=ys_resi, $
  CHARSIZE=plot_var.charsize
xyouts, ytpos, mean(!y.window), 'Residual', CHARSIZE=plot_var.charsize, /NORMAL, ALIGN=0.5, ORIENT=90, COLOR=plot_var.axiscolor

oplot, wave, goodresid, COLOR=FSC_Color(plot_var.residcolor), $
  PSYM=plot_var.residpsym, LINESTYLE=plot_var.residstyle, THICK=plot_var.residthick
oplot, wave, badresid, COLOR=FSC_Color(plot_var.badcolor), $
  PSYM=plot_var.residpsym, LINESTYLE=plot_var.residstyle, THICK=plot_var.residthick
oplot, wave, reform(factor*data[4,*]/data[2,*]), COLOR=FSC_Color(plot_var.sigmacolor), $
  PSYM=plot_var.sigmapsym, LINESTYLE=plot_var.sigmastyle, THICK=plot_var.sigmathick
oplot, wave,-reform(factor*data[4,*]/data[2,*]), COLOR=FSC_Color(plot_var.sigmacolor), $
  PSYM=plot_var.sigmapsym, LINESTYLE=plot_var.sigmastyle, THICK=plot_var.sigmathick
oplot, wave, wave*0., LINESTYLE=2, COLOR=FSC_Color(plot_var.sigmacolor)
if !d.name ne 'PS' then  if n_elements(smooth_resid) eq 1 then begin
;  do not put this decoration in PS output (should be in caption)
   x = mean(range)
   y = yr_resi[1] - 0.1 * (yr_resi[1]-yr_resi[0])
   xyouts, x, y, 'Smoothed: ' + strtrim(string(smooth_resid),2) + ' km/s', $
           COLOR=FSC_Color(plot_var.residcolor)
endif
if !d.name ne 'PS' then  if n_elements(smooth_flux) eq 1 then begin
   x = mean(range)
   y = yr_resi[1] + 0.3 * (yr_resi[1]-yr_resi[0])
   xyouts, x, y, 'Smoothed: ' + strtrim(string(smooth_flux),2) + ' km/s', COLOR=fsc_color('black')
endif

!p.multi=0

if !d.name ne 'PS' then device, DECOMPOSED=1

; free the heap variables allocated by ULY_SOLUT_SREAD
if size(solution, /TYPE) ne 8 then heap_free, solut, /PTR

end

