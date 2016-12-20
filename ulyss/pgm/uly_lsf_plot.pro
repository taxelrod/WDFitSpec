;+
; NAME:          
;                 ULY_LSF_PLOT
;
; PURPOSE: 
;                 Plot a LSF file
;
; USAGE:         
;                 uly_lsf_plot, <lsf_file> [, CONNECT=color][, /OVERPLOT]
;                               [, XTITLE=xtitle][, YSTYLE=ystyle][, TITLE=title]
;                               [, YMODE=ymode]
;                               [, WAVERANGE=<wr>] [, VRANGE=vrange][, SRANGE=srange][, PLOT_VAR=<plot_var>]
;
; DESCRIPTION:
;     Plot the relative velocity shift and instrumental velocity dispersion
;     as a function of the wavelength.
;     LSF (line-spread function) files are created by ULY_LSF.
;       
; ARGUMENTS:
;   <lsf_file>: 
;                 Name of the file containing the LSF.
;
; KEYWORDS:
;   CONNECT:
;                 Set this keyword to a color name to connect the data
;                 points. Example: connect='Red'
;                 Default is not to connect the data points.
;                 
;   PLOT_VAR:
;                 An input structure containing all the plot variables.
;                 This structure is returned by ULY_PLOT_INIT.
;
;   /OVERPLOT:
;                 Make plot over an existing one.
;
;   XTITLE:
;                 See the IDL documentation on "Graphic Keywords". 
;                 This keyword  is passed to the PLOT commands, it is defaulted to "Wavelength, A".
;
;   YMODE:
;                 Quantities represented on the Y axes:
;                 0 : cz and sigma in km/s
;                 1 : shift and FWHM in Angstrom 
;
;   YSTYLE:
;                 See the IDL documantation on "Graphic Keywords". 
;                 This keyword  is passed to the PLOT commands, it is defaulted to 16.
;
;   TITLE:
;                 Title of the plot.
;
;   /XLOG /YLOG
;                 Plot X (resp. Y) axis in logarithmic scale.
;
;   VRANGE SRANGE
;                 Y-axis range for respectively the cz and the sigma
;                 panels
;   WAVERANGE        	
;                 Wavelength range to be plotted
;   _EXTRA 
;                 including all the keyword present in ULY_PLOT_INIT
;
; EXAMPLE:
;     Plot a LSF with blue dots, and overplot a connected black line
;                 uly_lsf_plot, 'lsf.txt'
;                 uly_lsf_plot, 'lsfs.txt', /OVER, LINECOL='black', PSYM=0

;
; AUTHOR:         Mina Koleva, Philippe Prugniel, Antoine Bouchard 
;
;-
; CATEGORY:       ULY
;------------------------------------------------------------------------------
function uly_lsf_read, lsf_file

lsf = uly_solut_tread(lsf_file, STATUS=status)
if status ne 0 then begin
   message, '... try to use a simple tabular format', /INFO
   readcol,lsf_file, v, sig, h3, h4, centrWL, F='D,D,D,D,I'
   return, {v:v, sig:sig, h3:h3, h4:h4, centrWL:centrWL}
endif

if (size(lsf.losvd, /DIM))[0] ge 2 then begin
   v = reform(lsf.losvd[0,*])
   sig = reform(lsf.losvd[1,*])
endif
if (size(lsf.losvd, /DIM))[0] ge 3 then h3 = reform(lsf.losvd[2,*]) $
else h3 = 0*v
if (size(lsf.losvd, /DIM))[0] ge 4 then h4 = reform(lsf.losvd[3,*]) $
else h4 = 0*v
centrWL = 0*v
for k=0,n_elements(lsf)-1 do begin
   if n_elements(*(lsf[k]).hdr) gt 0 then begin
      wave_min = sxpar(*(lsf[k]).hdr, 'WAVE_MIN')
      wave_max = sxpar(*(lsf[k]).hdr, 'WAVE_MAX')
      centrWL[k] = (wave_min + wave_max) / 2
   endif else message, 'Invalid LSF file, cannot find wavelengths'
endfor

; Check if the LSF file was read sucessfully
if n_elements(v) eq 0 then begin
    message, 'Invalid LSF file '+lsf_file
endif

return, {v:v, sig:sig, h3:h3, h4:h4, centrWL:centrWL}
end

pro uly_lsf_plot, lsf_file, CONNECT=connect, PLOT_VAR=plot_var, OVERPLOT=overplot, $
                  XTITLE=xtitle, YSTYLE=ystyle, $
                  YMODE=ymode,                  $
                  WAVERANGE=waverange,          $
                  VRANGE=vrange, SRANGE=srange, $
                  TITLE=title,                  $
                  XLOG=xlog, YLOG=ylog,         $
                  _EXTRA=extra


lsf = uly_lsf_read(lsf_file)

if size(plot_var, /TYPE) ne 8 then begin
   if not keyword_set(psym) then psym=8
   plot_var = uly_plot_init(PSYM=psym, _STRICT_EXTRA=extra)
endif
if n_elements(ystyle) eq 0 then ystyle=18
if n_elements(ymode) eq 0 then ymode=0
if keyword_set(waverange) then if size(waverange, /n_elem) ne 2 then $
  message, '"WAVERANGE" should be 2 elements array...'

if n_elements(xtitle) eq 0 then xtitle = 'Wavelength, '+'!3' + STRING(197B) + '!X'
plotsym, 0, 1, /FILL

c = 299792.458d                 ; [km/s] Speed of light

; multiplot is not used there, but it might be because it is a nice program 
;DEPENDENCE:     multiplot.pro (by Fred Knight, astrolib)
;multiplot, [1,2]

plt = [0,0,0,0]
if min(lsf.v) lt max(lsf.v) then plt[0] = 1
if min(lsf.sig) lt max(lsf.sig) then plt[1] = 1
if min(lsf.h3) lt max(lsf.h3) then plt[2] = 1
if min(lsf.h4) lt max(lsf.h4) then plt[3] = 1
n = where(plt eq 1, nplt)

common plot_init, p1, x1, y1, p2, x2, y2
!p.multi=[0,1,nplt]

margin = 6./nplt ; total margin/ nb plot

case ymode of
    1 : begin
        ytitv='Shift, '+'!3' + STRING(197B) + '!X'
        vp = lsf.v/c * lsf.centrwl
        ytits='FWHM, '+'!3' + STRING(197B) + '!X'
        sigp = lsf.sig/c * lsf.centrwl * 2.355
    end
    else : begin  ; 0 or undefined cases
        ytitv='cz, km/s' 
        vp = lsf.v
        ytits='sigma, km/s'
        sigp = lsf.sig
    end
endcase

if not keyword_set(waverange) then waverange = minmax(lsf.centrwl)

ymargin=[-2,2+margin]
xtickn=replicate(' ', 30)
undefine, xtit

if plt[0] eq 1 then begin
   ymargin +=[margin,-margin]
   nplt -= 1
   if nplt eq 0 then begin
      undefine, xtickn
      xtit = xtitle
   endif
   if not keyword_set(overplot) then begin ;erase
      plot, lsf.centrwl, vp, $
            XLOG=xlog, XTICKNAME=xtickn, XTITLE=xtit, XR=waverange, XST=3, $
            YLOG=ylog, YRANGE=vrange, YSTYLE=ystyle, YTITLE=ytitv, YMARG=ymargin, $
            CHARSIZE=plot_var.charsize, $
            TITLE=title, BACK=plot_var.bgcolor, COLOR=plot_var.axiscolor, $
            /NODATA, FONT=0
      p1 = !P & x1 = !X & y1 = !Y ; save the graphic state in common variables
   endif else  begin
      if n_elements(p1) ne 0 then $
         !P = p1 & !X = x1 & !Y = y1 ; restore the graphic state (needed for /OVER)
   endelse
   if keyword_set(connect) then $
      oplot, lsf.centrwl, vp, COLOR=FSC_Color(connect), PSYM=0, LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick
   
   oplot, lsf.centrwl, vp, COLOR=plot_var.linecolor, PSYM=plot_var.psym, LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick
endif

if plt[1] eq 1 then begin
   ymargin += [margin,-margin]
   nplt -= 1
   if nplt eq 0 then begin
      undefine, xtickn
      xtit = xtitle
   endif
   if not keyword_set(overplot) then begin
      plot, lsf.centrwl, sigp, $
            XLOG=xlog, XTICKNAME=xtickn, XTITLE=xtit, XR=waverange, XST=3, $
            YLOG=ylog, YRANGE=srange, YSTYLE=ystyle, YTITLE=ytits, YMARG=ymargin, $
            CHARSIZE=plot_var.charsize, $
            BACK=plot_var.bgcolor, COLOR=plot_var.axiscolor, $
            /NODATA, FONT=0
      p2 = !P & x2 = !X & y2 = !Y ; save the graphic state in common variables
   endif else begin
      if n_elements(p2) ne 0 then $
         !P = p2 & !X = x2 & !Y = y2 ; restore the graphic state (needed for /OVER)
   endelse
   
   if keyword_set(connect) then $
      oplot, lsf.centrwl, sigp, COLOR=FSC_Color(connect), PSYM=0, LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick
   
   oplot, lsf.centrwl, sigp, COLOR=plot_var.linecolor, PSYM=plot_var.psym, LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick
endif   

if plt[2] eq 1 then begin
   ymargin +=[margin,-margin]
   nplt -= 1
   if nplt eq 0 then begin
      undefine, xtickn
      xtit = xtitle
   endif
   if not keyword_set(overplot) then begin
      plot, lsf.centrwl, lsf.h3, $
            XLOG=xlog, XTICKNAME=xtickn, XTITLE=xtit, XR=waverange, XST=3, $
            YLOG=ylog, YRANGE=srange, YSTYLE=ystyle, YTITLE='h3', YMARG=ymargin, $
            CHARSIZE=plot_var.charsize, $
            BACK=plot_var.bgcolor, COLOR=plot_var.axiscolor, $
            /NODATA, FONT=0
      p2 = !P & x2 = !X & y2 = !Y ; save the graphic state in common variables
   endif else begin
      if n_elements(p2) ne 0 then $
         !P = p2 & !X = x2 & !Y = y2 ; restore the graphic state (needed for /OVER)
   endelse
   
   if keyword_set(connect) then $
      oplot, lsf.centrwl, lsf.h3, COLOR=FSC_Color(connect), PSYM=0, LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick
   
   oplot, lsf.centrwl, lsf.h3, COLOR=plot_var.linecolor, PSYM=plot_var.psym, LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick
endif   

if plt[3] eq 1 then begin
   ymargin +=[margin,-margin]
   nplt -= 1
   if nplt eq 0 then begin
      undefine, xtickn
      xtit = xtitle
   endif
   if not keyword_set(overplot) then begin
      plot, lsf.centrwl, lsf.h4, $
            XLOG=xlog, XTICKNAME=xtickn, XTITLE=xtit, XR=waverange, XST=3, $
            YLOG=ylog, YRANGE=srange, YSTYLE=ystyle, YTITLE='h4', YMARG=ymargin, $
            CHARSIZE=plot_var.charsize, $
            BACK=plot_var.bgcolor, COLOR=plot_var.axiscolor, $
            /NODATA, FONT=0
      p2 = !P & x2 = !X & y2 = !Y ; save the graphic state in common variables
   endif else begin
      if n_elements(p2) ne 0 then $
         !P = p2 & !X = x2 & !Y = y2 ; restore the graphic state (needed for /OVER)
   endelse
   
   if keyword_set(connect) then $
      oplot, lsf.centrwl, lsf.h4, COLOR=FSC_Color(connect), PSYM=0, LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick
   
   oplot, lsf.centrwl, lsf.h4, COLOR=plot_var.linecolor, PSYM=plot_var.psym, LINESTYLE=plot_var.linestyle, THICK=plot_var.linethick
endif   

!p.multi=0


end
;--- end ---
