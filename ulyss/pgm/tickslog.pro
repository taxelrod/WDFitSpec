;-------------------------------------------------------------
;+
; NAME:
;        TICKSLOG 
;
; PURPOSE:
;        Compute values for the ticks of a logarithmic axis
;
; USAGE:
;        result = LOGLEVELS([range | MIN=min,MAX=max], [NMIN=nmin])
;
; DESCRIPTION:
;        Compute the position of the ticks for a graphic with a
;        logarithmic scale (using /XLOG in PLOT). Rounded values,
;        approximately logaritmically spaced, are generated.
;
;        The minimum number of ticks positions to generate may be
;        defined (by default it is 3), and the maximum number is
;        about twice the minimum.
;
; INPUTS:
;        RANGE -> A 2-element vector with the minimum and maximum
;            value to be returned. Only levels _within_ this range
;            will be returned. If RANGE contains only one element,
;            this is interpreted as MAX and MIN will be assumed as
;            3 decades smaller. RANGE superseeds the MIN and MAX
;            keywords. Note that RANGE must be positive definite
;            but can be given in descending order in which case
;            the labels will be reversed.
;
; KEYWORD PARAMETERS:
;        MIN, MAX -> alternative way of specifying a RANGE. If only
;            one keyword is given, the other one is computed as
;            3 decades smaller/larger than the given parameter.
;            RANGE superseeds MIN and MAX.
;
;        NMIN -> the minimum number of ticks to gnerate
;
; OUTPUTS:
;        A vector with "round" logarithmic values within the given
;        range. 
;
; EXAMPLE:
;        data = 5 + findgen(101)/100. *995
;        range = [min(data), max(data) ]
;        ticks = tickslog(range)
;        nticks = n_elements(ticks)
;        plot, data, /YLOG, YRANGE=range, YTICKS=nticks-1, YTICKV=ticks
;
; HISTORY:
;        2008/08/29, Philippe Prugniel, created from LOGLEVELS by Martin Schultz
;        mgs, 17 Mar 1999: VERSION 1.00
;
;-
; CATEGORY:  ULY_UTIL
;
; This program s a modification of LOGLEVELS (see copyright notice below)
; 
; Copyright (C) 1999, Martin Schultz, Harvard University
; This software is provided as is without any warranty
; whatsoever. It may be freely used, copied or distributed
; for non-commercial purposes. This copyright notice must be
; kept with any copy of this software. If this software shall
; be used commercially or sold as part of a larger package,
; please contact the author to arrange payment.
; Bugs and comments should be directed to mgs@io.harvard.edu
; with subject "IDL routine loglevels"
;-------------------------------------------------------------
function tickslog, range, MIN=mind, MAX=maxd, NMIN=nmin

if (n_elements(nmin) eq 0) then nmin = 3

; Make sure to have a valid range which is positive definite
; NOTE that range does not need to be sorted!
if (n_elements(mind) gt 0) then begin
    mind = mind[0]
    if (n_elements(maxd) eq 0) then maxd = mind*1000.
endif
if (n_elements(maxd) gt 0) then begin
    maxd = maxd[0]
    if (n_elements(mind) eq 0) then mind = maxd*0.001
endif

; still not defined, i.e. neither mind nor maxd given
if (n_elements(mind) eq 0) then begin
    mind = 0.1
    maxd = 100.
endif
; RANGE superseeds min and max
if (n_elements(range) eq 0) then range = [ mind,maxd ]
; one element for RANGE is interpreted as MAX
if (n_elements(range) eq 1) then range = [ range*0.001, range ]

thisrange = double(range) > 1.D-100
thisrange = thisrange(sort(thisrange))

; set lower range to 3 decades below upper range if it is zero
if (thisrange[0] lt 1.D-36) then thisrange[0] = thisrange[1]/1000.

; get log of ranges and decadal log
lrange = alog10(thisrange)

; set mode according to following rules:
; - range outside limits -> return decades
; - nmin exceeded -> return decades
; - [automatic] -> return decades if more than nmin decades
;                 otherwise 1,2,5,..

mode = 0                        ; return decades
if ((lrange[1]-lrange[0]) lt nmin) then mode = 1
if (thisrange[0] lt 1.D-15 OR thisrange[1] gt 5.D16) then mode = 0

if (mode) then begin
    dl = alog10(1D + 1d-7)
    lrange[0] -= dl
    lrange[1] += dl
; normalization
    deca = floor(min(double(lrange)))
    lrange -= deca 
    
    labels = [ 1.D+00, 2.D+00, 5.D+00] # [1d, 10d, 100d, 1000d, 10000d]
    llabels = alog10(labels)
    ind = where(llabels ge lrange[0] AND llabels le lrange[1], cnt)

    if cnt lt nmin then begin   ; try with 1-2-4
        lr = lrange +alog10(5d)
        labels /= 5d
        ind = where(llabels ge lr[0] AND llabels le lr[1], cnt)
    endif
    
    if cnt lt nmin then begin  ; try a denser more or less log spaced
        labels = [1D, 1.5D, 2D, 3D, 5D, 7D, 10D, 14D, 20D, 30D, 40D]
        llabels = alog10(labels)            
        ind = where(llabels ge lrange[0] AND llabels le lrange[1], cnt)
    endif
    if cnt lt nmin then begin   ; even denser
        labels = [ 1D, 1.2D, 1.5D, 2D, 3D, 4D, 5D, 6D, 7D, 8D, 10D, 12D, 15D, 20D, 25D, 30D, 40D]
        llabels = alog10(labels)            
        ind = where(llabels ge lrange[0] AND llabels le lrange[1], cnt)
    endif
    if cnt lt nmin and nmin gt 5 then begin ; even denser
        labels = [ 1D, 1.1D, 1.2D, 1.3D, 1.6D, 2D, 2.5D, 3D, 3.5D, 4D, 5D, 6D, $
                   7D, 8D, 9D, 10D, 11D, 13D, 15D, 18D, 21D, 25D, 30D, 35D, 40D]
        llabels = alog10(labels)            
        ind = where(llabels ge lrange[0] AND llabels le lrange[1], cnt)
    endif

    factor = 1d                 ; use a linear spacing
    while cnt lt nmin do begin
        if cnt lt nmin then begin
            labels = 1D + dindgen(20*factor) / factor
            llabels = alog10(labels)            
            ind = where(llabels ge lrange[0] AND llabels le lrange[1], cnt)
        endif
        if cnt lt nmin then begin
            labels = 1D + dindgen(20*factor)/2d / factor
            llabels = alog10(labels)            
            ind = where(llabels ge lrange[0] AND llabels le lrange[1], cnt)
        endif
        if cnt lt nmin then begin
            labels = 1D + dindgen(50*factor) / 5d / factor
            llabels = alog10(labels)            
            ind = where(llabels ge lrange[0] AND llabels le lrange[1], cnt)
        endif
        factor *= 10
        if factor gt 1d25 then break
    endwhile

; if range values are too close, return original range
    if (cnt eq 0) then return,range
    
; return reversed labels if range[0] gt range[1]
    if (range[0] gt range[1]) then $
      return,reverse(labels[min(ind):max(ind)]) *10d^deca $
    else $
      return,labels[min(ind):max(ind)] *10d^deca
    
endif else begin

    if (lrange[0] lt 0.) then lrange[0] = lrange[0] - 1.0D-6
    if (lrange[1] lt 0.) then lrange[1] = lrange[1] - 1.0D-6
    
    maxtick = max([2*nmin, 8])

    lfactor = 1d
    while lrange[1]-lrange[0] gt maxtick do begin ; at max maxtick ticks
        lrange *= lfactor/(lfactor+1)
        lfactor ++
    endwhile
    
    drange = fix(lrange)
    if (lrange[1] lt 0.) then drange[1] = drange[1] - 1
    
    exponent = (drange[0]+indgen(drange[1]-drange[0]+1)) * lfactor

    if (range[0] gt range[1]) then $
      return,reverse(10.D0^exponent)  $
    else $
      return,10.D0^exponent
endelse

end
