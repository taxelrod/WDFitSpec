;+
; NAME:
;              ULY_CHIMAP_INIT
; PURPOSE:
;              Initialize the chi square map structure
;
; USAGE:
;              map = uly_chimap_init(<n_xnode>, <n_ynode>)
;
; ARGUMENTS:
;	X_NODE, Y_NODE
;	       The number of x and y nodes of the chi square map structure.
;
; RETURN VALUE:
;              Initialized chi square map structure. Which contains:
;       .map   	2d array with the computed chi2 values over xnode/ynode
;       .xnode 	nodes on the abscissa where the chi2 is computed 
;       .ynode 	nodes on the ordinate where the chi2 is computed 
;       .xtitle	name of the x variable (for example age, tempreture...)
;       .ytitle	name of the y variable (for example velocity, Fe/H...)
;       .xunit 	units of the x variable (for example Myr, C...)
;       .yunit 	units of the y variable (for example km/s, dex...)
;
; AUTHOR:
;              Antoine Bouchard
; HISTORY:
;              2008/02/08
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
function uly_chimap_init, n_xnode, n_ynode
return, {map:fltarr(n_xnode,n_ynode)*!VALUES.F_NAN, 	$
         xnode:fltarr(n_xnode)*!VALUES.F_NAN, 		$
         ynode:fltarr(n_ynode)*!VALUES.F_NAN, 		$
         xtitle:'', 					$
	 ytitle:'',					$
         xunit:'', 					$
	 yunit:''}
end
