;+
; NAME:          ULY_PLOT_SAVE
;
; PURPOSE:       Make a hardcopy of a plot.
; 
; USAGE:
;                uly_plot_save, <filename>
;                [, /JPG]
;                [, /PNG]
;                [, /TIFF]
;
; DESCRIPTION:
;                After making a plot, use this command to keep a
;                hardcopy of it.
;
; ARGUMENTS:
;    <filename>  Input
;                Name of the file to be writing.
;
; KEYWORDS:
;    /JPG, /PNG, /TIFF
;                Selects the output format. Default is JPG.
;
; HISTORY:        
;    Creation    Antoine Bouchard, 5/09/2008
;
;------------------------------------------------------------
; CATEGORY:     ULY

pro uly_plot_save, filename, JPG=jpg, PNG=png, TIFF=tiff

  device, DECOMPOSED=1
  image = tvrd(TRUE=1)

  if not keyword_set(jpg) $
    and not keyword_set(gif) $
    and not keyword_set(png) $
    and not keyword_set(tiff) then jpg=1

  if keyword_set(jpg) then write_jpeg, filename, image, QUALITY=100, TRUE=1
  if keyword_set(png) then write_png, filename, image
  if keyword_set(tiff) then write_tiff, filename, image

  device, DECOMPOSED=0

end
