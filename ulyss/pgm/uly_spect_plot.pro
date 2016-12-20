;+
; NAME:          ULY_SPECT_PLOT
;
; PURPOSE:
;                Plot 1d or 2d fits files
;
; USAGE:
;		 uly_spect_plot, <spectrum>
;                [, /OVERPLOT]
;                [, WAVERANGE=<wr>][, YRANGE=<yr>]
;                [, NORM=<norm>]
;                [, ONED=<oned>]
;                [, /XLOG][, /YLOG]
;                [, TITLE=title]
;                [, PLOT_VAR=plot_var]
;                [, /QUIET]
;		 
;
; DESCRIPTION:    
;      This program plots 1d or 2d fits files. It can be used to visualize 
;      the data that you may further want to analyse with ULySS package.
;       
;      Note: for making POSTSCRIPT files you need to setup the device 'ps'. 
;      It is recommended that you plot what you would like to appear in
;      your .ps file, then type "set_plot, 'ps'" and "device, fileout, /COLOR"
;      (more options can be added to the command "device"), recall
;      all your commands, type "device,/close" and "set_plot, 'x'".
;
;
; ARGUMENTS:
;      <spectrum>	Input
;                       FITS file name or 'spect' structure containing either a
;                       1d or 2d spectra.
;                       A 'spect' structure is returned by, e.g. ULY_SPECT_READ.
; KEYWORDS:
;      /OVERPLOT
;                       When this keyword is set, plots the spectra over an
;                       existing plot. If color is set to a color name (eg.
;                       'Red', 'Green', 'Blue'). The colors are defined as
;                       in fsc_color.pro. Default is 'Blue'.
;      WAVERANGE        	
;                       Wavelength range to be plotted
;      YRANGE        	
;                       Flux range to be plotted
;      NORM      	
;                       2 element vector indicating the wavelength
;                       range in which the flux will be normalised. Eg.
;                       NORM=[3990.,4010.] will result in having the
;                       average flux within this region set to 1.
;      ONED
;                       Position of the line to plot from a multidimensional 
;                       dataset is given, set this value to the line number 
;                       to extracted with ULY_SPECT_EXTRACT. Default - 1st line
;      /XLOG /YLOG
;                       Plot X (resp. Y) axis in logarithmic scale.
;      TITLE
;                       Title of the plot. Default: Name of the FITS file
;      PLOT_VAR
;                       An input structure containing all plot variables.
;                       This structure is returned by ULY_PLOT_INIT
;      ...Extra keywords
;                       If PLOT_VAR is undefined, all extra keywords are
;                       passed to ULY_PLOT_INIT.
;      QUIET
;                       Verbosity control.
;
; DEPENDENCIES:
;                fsc_color.pro by David Fanning
;                uly_spect_read.pro from ULySS
;
; HISTORY:        
;     Creation   Mina Koleva, 22/07/2008
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
pro uly_spect_plot, filein, OVERPLOT=over,                 $
                    WAVERANGE=waverange, YRANGE=fluxrange, $
                    NORM=norm, ONED=oned,                  $
                    XLOG=xlog, YLOG=ylog,                  $
                    TITLE=title,                           $
                    PLOT_VAR=plot_var, QUIET=quiet, _EXTRA=extra


; checks
if size(filein,/TYPE) eq 8 then begin
    sp = uly_spect_extract(filein, ONED=oned)
    file = '<spectrum>'
endif else begin
    if size(filein, /TYPE) ne 7 then $
      message,'The input must be a file name or a spect structure'
    file = filein
endelse

if file ne '<spectrum>' then begin
    testinp = file_test(file)
    if testinp eq 0 then $
      message,'The file:'+file+' does not exist. Give correct name for the input file.'
    sp = uly_spect_read(filein, QUIET=quiet) ; read the spectrum
    sp = uly_spect_extract(sp, ONED=oned, /OVERWRITE)
endif 
if keyword_set(over) then if size(over, /type) ne 7 then over='Blue'
  ;message, '"Over" should be string with a color name...'
if keyword_set(waverange) then if size(waverange, /n_elem) ne 2 then $
  message, '"WAVERANGE" should be 2 elements array...'
if keyword_set(fluxrange) then if size(fluxrange, /n_elem) ne 2 then $
  message, '"YRANGE" should be 2 elements array...'
if keyword_set(norm) then if size(norm, /n_elem) ne 2 then $
  message, '"NORM" should be 2 elements array...'
;if keyword_set(line) then if size(line, /n_elem) ne 1 then $
;  message, '"LINE" should be an 1 element integer...'
  
if not keyword_set(plot_var) or size(plot_var, /type) ne 8 then begin
  if not keyword_set(linecolor) then linecolor='Black'
  plot_var = uly_plot_init(LINECOLOR=linecolor, _extra=extra)
endif

; compute wl array
if sp.sampling le 1 then begin
   wl = (sp.start+sp.step*findgen(n_elements(*sp.data)))
   if sp.sampling eq 1 then wl = exp(wl)
endif else if sp.sampling eq 2 then wl = *sp.wavelen

; check if a bad pixel mask is defined
if n_elements(*sp.goodpix) gt 0 then goodpix = *sp.goodpix $
else goodpix = where(finite(*sp.data))

; check what keywords are set
if keyword_set(norm) then $
  *sp.data = *sp.data/mean((*sp.data)[where(wl ge norm[0] and wl le norm[1] )])
if not keyword_set(waverange) then waverange = minmax(wl)
if not keyword_set(fluxrange) then fluxrange = minmax((*sp.data)[goodpix])


if not keyword_set(over) then begin                     ;to over plot a spectrum
    if n_elements(title) eq 0 then title=sp.title
    plot, wl, *sp.data, XR=waverange, YR=fluxrange, XSTYLE=3, YSTYLE=3, $
      XTITLE='Wavelength, '+'!3' + STRING(197B) + '!X', YTITLE='Flux', TITLE=title, $
      CHARSIZE=plot_var.charsize, $
      BACK=plot_var.bgcolor, COLOR=plot_var.axiscolor, $
      XLOG=xlog, YLOG=ylog, /NODATA
endif

goodspec = *sp.data*!VALUES.F_NAN
badspec  = *sp.data*!VALUES.F_NAN

goodspec[goodpix] = (*sp.data)[goodpix]
badpix = uly_plot_padpixelrange(where(finite(goodspec, /NAN), badcount) )

if badcount gt 0 then begin
  badpix = uly_plot_padpixelrange(badpix)	
  badspec[badpix] = (*sp.data)[badpix]
endif

oplot, wl, goodspec, COLOR=plot_var.linecolor, $
  PSYM=plot_var.psym, LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick
if badcount gt 1 then $
  oplot, wl, badspec, COLOR=FSC_Color(plot_var.badcolor), $
    PSYM=plot_var.psym, LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick

uly_spect_free, sp

end
