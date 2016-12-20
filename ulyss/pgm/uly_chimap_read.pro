;+
; NAME:
;              ULY_CHIMAP_READ
; PURPOSE:
;              Read in a chi square map file
;
; USAGE:       chi = uly_chimap_read(<file>)
;
; ARGUMENTS:
;    FILE:     Name of the file saved by ULY_CHIMAP_WRITE (string)
;
; RETURN VALUE:
;              Chi-square map structure. The structure containing:
;	       .map   2d array with the computed chi2 values over xnode/ynode
;	       .xnode nodes on the abscissa where the chi2 is computed 
;	       .ynode nodes on the ordinate where the chi2 is computed 
;	       .xtitle  name of the x variable (for examp. age, tempreture...)
;	       .ytitle  name of the y variable (for examp. velocity, Fe/H...)
;	       .xunit units of the x variable (for examp. Myr, C...)
;	       .yunit units of the y variable (for examp. km/s, dex...)

;
; DEPENDENCE:  ULY_CHIMAP_INIT
;
; HISTORY:
;     creation: Antoine Bouchard 2008/02/08
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
function uly_chimap_read, file

fits_read, file, map, hdr, EXTEN_NO=0, MESSAGE=mess0
fits_read, file, xnode, EXTEN_NO=1, MESSAGE=mess1
fits_read, file, ynode, EXTEN_NO=2, MESSAGE=mess2

if mess0 ne '' or mess1 ne '' or mess2 ne ''  then begin
    message, 'Failed to read FITS file...'+mess0+mess1+mess2, /INFO
    return, 0
endif

chimap     	= uly_chimap_init(n_elements(xnode), n_elements(ynode))
chimap.map 	= map
chimap.xnode 	= xnode
chimap.ynode   	= ynode

chimap.xtitle 	= sxpar(hdr, 'CTYPE1')
chimap.ytitle 	= sxpar(hdr, 'CTYPE2')

chimap.xunit 	= sxpar(hdr, 'CUNIT1')
chimap.yunit 	= sxpar(hdr, 'CUNIT2')


return, chimap

end
;--- end ---
