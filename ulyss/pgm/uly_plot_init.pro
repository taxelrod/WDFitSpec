;+
; NAME:         ULY_PLOT_INIT
; 
; PURPOSE:      Initialize plotting variables
;
; USAGE:
;		uly_plot = uly_plot_init([  XSIZE=xsize]	   $
;                                        [, YSIZE=ysize]	   $
;                                        [, BGCOLOR=bgcolor]       $
;                                        [, AXISCOLOR=axiscolor]   $
;                                        [, LINECOLOR=linecolor]   $
;                                        [, MODELCOLOR=modelcolor] $
;                                        [, POLYCOLOR=polycolor]   $
;                                        [, RESIDCOLOR=residcolor] $
;                                        [, SIGMACOLOR=sigmacolor] $
;                                        [, PSYM=psym]             $
;                                        [, MODELPSYM=modelpsym]   $
;                                        [, POLYPSYM=polypsym]     $
;                                        [, RESIDPSYM=residpsym]   $
;                                        [, SIGMAPSYM=sigmapsym]   $
;                                        [, LINESTYLE=linestyle]   $
;                                        [, MODELSTYLE=modelstyle] $
;                                        [, POLYSTYLE=polystyle]   $
;                                        [, RESIDSTYLE=residstyle] $
;                                        [, SIGMASTYLE=sigmastyle] $
;                                        [, LINETHICK=linethick]   $
;                                        [, MODELTHICK=modelthick] $
;                                        [, POLYTHICK=polythick]   $
;                                        [, RESIDTHICK=residthick] $
;                                        [, SIGMATHICK=sigmathick] $
;                                        [, CHARSIZE=charsize]     $
;                                        [, CHARTHICK=charthick]   $
;                                        [, HALFTONECT=halftonect] $
;                                        [, /HIDEHALFTONE]         $
;                                        [, /HIDECONTOUR]          $
;                                        [, /SHOWSIGMA]            $
;                                        [, /HIDEMINIMA]           
;
; DESCRIPTION:
;		Creates a structure containing all plotting variables, such
; 		as colors, character sizes, etc.
;		This task is called by the plotting tasks directly so there
;		is no need to explicitly call it unless you want to change
;		the plotting behaviour.
;               The list of color names can be obtained as: 
;                       print, FSC_COLOR(/NAMES)
; 
; ARGUMENTS:
;		XSIZE, YSIZE
;			If no plotting window exists, create a window of this
;			size, in pixels. If a plotting window already
;			exists, these keywords have no effect, the
;			opened window will be used. 
;
; 		BGCOLOR
;                       Color index (integer) or name (character string).
;			Background color, default is 'White'.
;
;		AXISCOLOR
;                       Color index (integer) or name (character string).
;			Color for the axes and decoration, default is 'Black'.
;
;		LINECOLOR
;                       Color index (integer), or name (character string).
;			Main plotting line color, default is 'Blue'.
;
;		MODELCOLOR, POLYCOLOR, RESIDCOLOR, and SIGMACOLOR
;			Plotting colors for the model, polynomial, residual
;			and 1 sigma level respectively. 
;			Defaults are: Blue for model, red for polynomial,
;			black for residual, green for sigma.
;
;               PSYM, MODELPSYM, POLYPSYM, RESIDPSYM, SIGMAPSYM
;			Plotting symbols for the line, model, polynomial, residual
;			and 1 sigma level respectively.
;			Defaults are 0.
;
;		LINESTYLE, MODELSTYLE, POLYSTYLE, RESIDSTYLE, SIGMASTYLE
;			Line style parameter for data, model, polynomial,
;			residual and sigma level respectively, default is 0.
;                       0	Solid
;                       1	Dotted
;                       2	Dashed
;                       3	Dash Dot
;                       4	Dash Dot Dot
;                       5	Long Dashes 
;
;		LINETHICK, MODELTHICK, POLYTHICK, RESIDTHICK, SIGMATHICK
;			Line thickness parameter for data, model, polynomial,
;			residual and sigma level respectively, default is 1.
;
;		CHARSIZE
;			Axis character size, default is 1.4.
;        
;               CHARTHICK
;                       Axis character thick, default is 1.
;
;		HALFTONECT
;			Color table for halftone overlay in Chi2 maps.
;			Default is 3 (Red Temperature). Use XLOADCT to get
;			all available color tables. 
;
;		/HIDEHALFTONE
;			Set this keyword to hide the halftone overlay in
;			Chi2 maps.
;
;		/HIDECONTOUR
;			Set this keyword to hide contour levels in Chi2
;			maps.
;
;		/SHOWSIGMA
;			Set this keyword to show sigma levels.
;
;		/HIDEMINIMA
;			Set this keyword to hide the local and global
;			minima in Chi2 maps.
;
; EXAMPLE:
;	
;		Define where to get the data:
;			galaxy = uly_root+'/data/VazMiles_z-0.40t07.94.fits'
;		
;		Initialize the plotting variables and setting a green on
;		red plot and log y axis:
;			plot_var = uly_plot_init(linecolor='Blue', bgcolor='Red', /ylog)
;
;		Set the linestyle to be hashed:
;			plot_var.linestyle=2
;
;		Set the axis color to be purple.
;			plot_var.axiscolor=FSC_Color('Purple')
;
;		Make the plot.
;			uly_spect_plot, galaxy, plot_var=plot_var
; 		
;		Overplot a convolved version of that plot, in blue.
;			plot_var.linecolor=FSC_Color('Green')
;			plot_var.linestyle=0
;			uly_spect_plot, galaxy, conv=[10000., 300., 0., 0.], /overplot, plot_var=plot_var
;
;		Contemplate a truly artistic choice of color!
;
; DEPENDENCE:   fcs_color.pro by David Fanning
;
; AUTHOR:
; 		Antoine Bouchard
; 
; HISTORY:
;		2008/08/15	Creation
; -
; CATEGORY:     ULY
;----------------------------------------------------------------------------
function uly_plot_init, XSIZE=xsize,		     $
			YSIZE=ysize,		     $
			BGCOLOR=bgcolor,             $
                        AXISCOLOR=axiscolor,         $
                        HALFTONECT=halftonect,       $
                        LINECOLOR=linecolor,         $
                        MODELCOLOR=modelcolor,       $
                        POLYCOLOR=polycolor,         $
                        RESIDCOLOR=residcolor,       $
                        SIGMACOLOR=sigmacolor,       $
			BADCOLOR=badcolor,	     $
			PSYM=psym,                   $
			MODELPSYM=modelpsym,         $
			POLYPSYM=polypsym,           $
			RESIDPSYM=residpsym,         $
			SIGMAPSYM=sigmapsym,         $
                        LINESTYLE=linestyle,         $
                        MODELSTYLE=modelstyle,       $
                        POLYSTYLE=polystyle,         $
                        RESIDSTYLE=residstyle,       $
                        SIGMASTYLE=sigmastyle,       $             
                        LINETHICK=linethick,         $
                        MODELTHICK=modelthick,       $
                        POLYTHICK=polythick,         $
                        RESIDTHICK=residthick,       $
                        SIGMATHICK=sigmathick,       $
                        CHARSIZE=charsize,           $
                        CHARTHICK=charthick,         $
                        HIDETITLE=hidetitle,         $
                        HIDEHALFTONE=hidehalftone,   $
                        HIDECONTOUR=hidecontour,     $
                        SHOWSIGMA=showsigma,         $
                        HIDEMINIMA=hideminima,       $
			_EXTRA=extra

common plot_init, p1, x1, y1, p2, x2, y2

if not keyword_set(xsize) or size(xsize, /TYPE) ne 2 then xsize = 700
if not keyword_set(ysize) or size(ysize, /TYPE) ne 2 then ysize = 500
if not keyword_set(halftonect) or size(halftonect, /TYPE) ne 2 then halftonect=39
loadct, halftonect, /SILENT

if not keyword_set(bgcolor) then bcolor=FSC_Color("White") $
else if size(bgcolor, /TYPE) eq 7 then bcolor = FSC_Color(bgcolor) $
else bcolor = bgcolor

if not keyword_set(axiscolor) then acolor=FSC_Color("Black") $
else if size(axiscolor, /TYPE) eq 7 then acolor = FSC_Color(axiscolor) $
else acolor = axiscolor

if not keyword_set(linecolor) then lcolor = FSC_Color("Blue") $
else if size(linecolor, /TYPE) eq 7 then lcolor = FSC_Color(linecolor) $
else lcolor = linecolor

if not keyword_set(modelcolor) or size(modelcolor, /TYPE) ne 7 then modelcolor="Blue"
if not keyword_set(polycolor)  or size(polycolor,  /TYPE) ne 7 then polycolor="Cyan"
if not keyword_set(residcolor) or size(residcolor, /TYPE) ne 7 then residcolor="Black"
if not keyword_set(sigmacolor) or size(sigmacolor, /TYPE) ne 7 then sigmacolor="Green"
if not keyword_set(badcolor)   or size(badcolor,   /TYPE) ne 7 then badcolor="Red"
if not keyword_set(psym)       or size(psym,       /TYPE) ne 2 then psym=0
if not keyword_set(modelpsym)  or size(modelpsym,  /TYPE) ne 2 then modelpsym=10
if not keyword_set(polypsym)   or size(polypsym,   /TYPE) ne 2 then polypsym=10
if not keyword_set(residpsym)  or size(residpsym,  /TYPE) ne 2 then residpsym=10
if not keyword_set(sigmapsym)  or size(sigmapsym,  /TYPE) ne 2 then sigmapsym=0
if not keyword_set(linestyle)  or size(linestyle,  /TYPE) ne 2 then linestyle=0
if not keyword_set(modelstyle) or size(modelstyle, /TYPE) ne 2 then modelstyle=0
if not keyword_set(polystyle)  or size(polystyle,  /TYPE) ne 2 then polystyle=0
if not keyword_set(residstyle) or size(residstyle, /TYPE) ne 2 then residstyle=0
if not keyword_set(sigmastyle) or size(sigmastyle, /TYPE) ne 2 then sigmastyle=0
if not (keyword_set(linethick)  or $
        size(linethick,  /TYPE) eq 2 or size(linethick,  /TYPE) eq 4) then $
           if !p.charthick eq 0 then linethick=1 else linethick=!p.charthick
if not (keyword_set(modelthick) or $
        size(modelthick, /TYPE) eq 2 or size(modelthick, /TYPE) eq 4) then $
           if !p.charthick eq 0 then modelthick=1 else modelthick=!p.charthick
if not (keyword_set(polythick)  or $
        size(polythick,  /TYPE) eq 2 or size(polythick,  /TYPE) eq 4) then $
           if !p.charthick eq 0 then polythick=1 else polythick=!p.charthick
if not (keyword_set(residthick) or $
        size(residthick, /TYPE) eq 2 or size(residthick, /TYPE) eq 4) then $
           if !p.charthick eq 0 then residthick=1 else residthick=!p.charthick
if not (keyword_set(sigmathick) or $
        size(sigmathick, /TYPE) eq 2 or size(sigmathick, /TYPE) eq 4) then $
           if !p.charthick eq 0 then sigmathick=1 else sigmathick=!p.charthick
if not (keyword_set(charthick)  or $
        size(charthick,  /TYPE) eq 2 or size(charthick,  /TYPE) eq 4) then $
           if !p.charthick eq 0 then charthick=1 else charthick=!p.charsize
if not (keyword_set(charsize)   or $
        size(charsize,   /TYPE) eq 2 or size(charsize,   /TYPE) eq 4) then $
           if !p.charsize eq 0 then charsize=1 else charsize=!p.charsize
if keyword_set(hidehalftone) then showhalftone = 0 else showhalftone = 1
if keyword_set(hidecontour)  then showcontour = 0  else showcontour = 1
if keyword_set(showsigma)    then showsigma = 1    else showsigma = 0
if keyword_set(hideminima)   then showminima = 0   else showminima = 1

return, { $
          xsize:xsize,			$
          ysize:ysize,			$
          bgcolor:bcolor,		$
          axiscolor:acolor,		$
          halftonect:halftonect,	$
          linecolor:lcolor,		$
          modelcolor:modelcolor,	$
          polycolor:polycolor,		$
          residcolor:residcolor,	$
          sigmacolor:sigmacolor,	$
          badcolor:badcolor,		$
          psym:psym, 			$
          modelpsym:modelpsym,		$
          polypsym:polypsym,		$
          residpsym:residpsym,		$
          sigmapsym:sigmapsym,		$
          linestyle:linestyle,		$
          modelstyle:modelstyle,	$
          polystyle:polystyle,		$
          residstyle:residstyle,	$
          sigmastyle:sigmastyle,	$
          linethick:linethick,		$
          modelthick:modelthick,	$
          polythick:polythick,		$
          residthick:residthick,	$
          sigmathick:sigmathick,	$
          charsize:charsize,		$
          charthick:charthick,          $
          showhalftone:showhalftone,	$
          showcontour:showcontour,	$
          showsigma:showsigma,		$
          showminima:showminima	        $
        }
end
