function uly_plot_padpixelrange, pixelrange
;   ULY_PLOT_PADPIXELRANGE  (internal routine)
;   Takes a pixel range that was generated using the where() command and
;   adds 1 extra pixel on either side of each pixels found in the array. 
;   There is a check for sequence of pixels to avoid having redundancy.
;   Example:
;     If pixelrange contains [3,4,5,42] then this function computes a new
;     array by appending the edges, ie the values 2 and 6 (the edges of the
;     first serie) and 41 and 43 (the edges of the last value). The new
;     array becomes: [2,3,4,5,6,41,42,43]
;
;   The reason for this is that when the residuals are smoothed in
;   uly_solut_splot, then we want the masked regions (plotted in yellow by
;   default) to also be 'smoothed'.
;
; AUTHOR:
;               Antoine Bouchard
; HISTORY:
;               2009/05/27  created

limits = where(shift(pixelrange,1)-pixelrange ne -1 )
newpixelrange = [pixelrange,pixelrange[limits]-1, pixelrange[(limits-1)]+1]
  ; After detecting any sequences of pixel, the last value+1 and the next
  ; value-1 are appended.

return, newpixelrange[sort(newpixelrange)]
  ; The sort is required to mimick the output of the where() function; the
  ; appended values were done at the end of the array, so it needs to be
  ; reordered.
end

