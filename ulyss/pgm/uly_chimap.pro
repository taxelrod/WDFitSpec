;+
; NAME:
;               ULY_CHIMAP
;
; PURPOSE:
;               Create a chi2 map
;
; USAGE:
;               chimap = uly_chimap(<spectrum>, <cmp>, [proj]
;                                   [, VELSCALE=<velscale>]
;                                   [, LMIN=<lmin>][, LMAX=<lmax>]
;                                   [, SG=<sg>][, DG=<dg>]
;                                   [, KMOMENT=<kmoment>][, KFIX=<kfix>]
;                                   [, ERR_SP=<err_sp>]
;                                   [, SNR=<snr>]
;                                   [, XNODE=<xnode>][, YNODE=<ynode>]
;                                   [, ADEGREE=<adegree>][, MDEGREE=<mdegree>]
;                                   [, PLOT=<plot>][, /PS]
;                                   [, /CLEAN]
;                                   [, QUIET=<quiet>]
;                                   [, FILE_OUT=<file_out>])
;
; DESCRIPTION:
;      Explore the the parameter space in search for the best fitting model. 
;      At each node of the grid the other parameters are optimized.
;
; ARGUMENTS:
;   <spectrum>  Input
;               File name, or 'spectrum' structure containing 1d spectrum
;               A 'spectrum' structure is returned by,
;               e.g. ULY_SPECT_READ.
;
;   <cmp>       Input
;               Structure returned from user function for defining
;               component (see for exampl ULY_SSP/ ULY_TGM)
;               it can be a vector of components (like [cmp1, cmp2, ...])
;
;   [proj]      Optional Input
;		The 2x2 matrix defining the projection of the chimap
;               in the space of n parameters to be fitted, for example
;               to plot the age of the first
;               comp. (so it is component = 0, parameter = 0), and the
;               fraction of the second component (burst) (comp = 1,
;               param=-1), the projection should be [[0,0],[1,-1]]. 
;               In a case that components and parameters are
;               unknown leave this parameter free.
;
;
;
;
; KEYWORDS:
;   [VELSCALE]    [km/s]
;                 Velocity scale for which the chi square is computed.
;                 Default: The model velocity scale.
;
;   [LMIN]        [ANGSTROM]
;                 Standard 'good region' descriptor. No default.
;
;   [LMAX]        [ANGSTROM]
;                 Standard 'good region' descriptor. No default.
;
;   [SG]          [km/s]
;                 Guess for the redshift.
;                 The input spectrum is shifted by this quantity when
;                 it is read (the guess for the minimization is set
;                 to 0, and cannot by changed).
;                 This value must be quite precise: An error by more than
;                 3 or 5 times the velocity dispersion may prevent the
;                 fit to converge.
;
;   [DG]          [km/s] 
;                 Guess or fixed value for the velocity dispersion.
;                 This guess is not very critical, and the default value 
;                 of 1 pixel is generally satisfactory.
;
;   [ERR_SP]      
;                 1d file containing the error of the spectrum to analyse;
;                 relevant only in the case where a file name is given for
;                 <spectrum>.
;
;   [KMOMENT]        
;                 Number of terms of the LOSVD.
;                 The terms are in the order: [cz, sigma, h3, h4, h5,
;                 h6]. By default, KMOMENT=2, i.e. a Gaussian LOSVD
;                 is fitted.
;
;   [KFIX]        Default: The kinematics parameters are free.
;                 Array used to fix some of the parameters of the LOSVD,
;                 0 means that the parameter is free; 1 that it is fixed.
;                 The parameters are specified in the following order:
;                 [cz, sigma, h3, h4, h5, h6]. cz and sigma must be
;                 given in km/s.
;
;   [SNR]         default=100
;                 Mean signal to noise ratio the of analyzed spectrum.
;                 This parameter is used to derive the errors if the
;                 error spectrum is not attached to the input spectrum
;                 and if ERR_SP is not given. 
;                 It will generate a constant error spectrum.
;
;   [XNODE]       
;                 X values at which the chi square is computed.
;                 If an integer scalar is given, take is as the number of nodes
;                 and use the limits of the grid.
;                 Default: a 10 points regular spaced grid between the limits
;
;   [YNODE]       
;                 Y at which the chi square is computed.
;                 If an integer scalar is given, take is as the number of nodes
;                 and use the limits of the grid.
;                 Default: a 10 points regular spaced grid between the limits
;
;   [ADEGREE]     
;                 degree of the *additive* Legendre polynomial
;                 Default: -1
;
;   [MDEGREE]     
;                 degree of the *multiplicative* Legendre polynomial
;                 Default: 10
;
;   [/PLOT]       
;                 Plot the solution of a fit with ULY_SOLUT_SPLOT as well
;                 as the chi square map with ULY_CHIMAP_PLOT at each
;                 grid node. Default is not to plot at each node.
;
;   [/PS] or [PS='filename.ps']
;      	          Makes a postcript hardcopy of the plot. When using /PS the
;                 filename is called 'chimap.ps'. To use a different name,
;                 specify it with PS='other_name.ps'. 
;
;   [/CLEAN]      set this keyword to use the iterative kappa-sigma clipping.
;
;   [QUIET]       
;                 Verbosity control
;  
;   [FILE_OUT]    
;                 Ouput filename in which the Chisq map is written.
;                 Formally this file is in fits format. 
;                 Default is to not write the file.
;
; RETURN VALUE:
;                 Returns the computed chi square map structure.
;
; EXAMPLE:
;    file_test = uly_root+'/data/VazMiles_z-0.40t07.94.fits'
;    cmp = uly_ssp()
;    map = uly_chimap(file_test, cmp, [[0,0],[0,1]], /quiet)
;
; AUTHOR:
;               Mina Koleva, Antoine Bouchard
;
; HISTORY:
;          2008/07/31 creation, using earlier version from AB
;
; KNOWN BUGS:   If you want to fix the fraction better give
;               nodes. Usually the values of the fraction are not known
;               apriori. So it is better to fit once your spectrum to see
;               the order of the magnitude of the fraction and then
;               make chi2 maps.
;-
; CATEGORY:    ULY

;########################################################################
;+
;
; NAME:
;               ULY_CMP_SELECT
;
; PURPOSE:
;               Make a 2d vector describing an axis of the parameters space
;
; USAGE:
;               axis = uly_cmp_select(cmp)
;
; DESCRIPTION:
;      Called from chi2/convergence map program it is used to select a 
;      projection.
;
;      An axis of the parameters space is described by a two-elements
;      vector. The first element indicates the number of the component
;      or is -1 if the parameter is one of those of the LOSVD.
;      The second element is the order number of the parameter in the 
;      component or LOSVD; it is -1 for the weight.
;
;      The user is prompted for information describing the axis.
;
; ARGUMENT:
;          CMP  previously defined component
;
; RETURN VALUE:
;          AXIS a 2 elements vector: code for the projection of the component  
;
; EXAMPLE:
;               cmp = uly_ssp()
;               axis = uly_cmp_select(cmp)
;
; AUTHOR:
;               Mina Koleva
;
; HISTORY:
;          2008/07/31 creation
;-
; CATEGORY:    ULY

function uly_cmp_select, cmp

ncomp = n_elements(cmp)
if ncomp eq 0 then begin
    print, 'Usage: uly_cmp_select, <cmp>'
    message, 'Miss argument <cmp>', /INFO
    return, 0
endif

print, 'For LOSVD type: -1 '
for i = 0, ncomp-1 do $
  print, 'For ', cmp[i].name, ' type: ', i, FORMAT='(3a, i1)'

comp = ' '
read, comp, prompt = 'Chose a number of component: '

if comp ne -1 then begin
    npars = n_elements(*cmp[comp].para)
    print, 'For fraction type: -1 '
    for i = 0, npars-1 do begin
        if (*cmp[comp].para)[i].name ne '' then $
          print, 'For ', (*cmp[comp].para)[i].name, ' type: ', i, FORMAT='(3a, i1)'
    endfor
endif else begin
    losvd_name = ['velocity ', 'sigma ']
    for i = 0, 1 do print, 'For ', losvd_name[i], 'type: ', i, FORMAT='(3a, i1)'
endelse

param = ' '  
read, param, prompt = 'Chose a number of parameter: '

return, fix([comp, param])

end

;########################################################################
;+
; NAME:
;              ULY_CHIMAP_WRITE
; PURPOSE:
;              Write a chi square map in a file
;
; USAGE:
;              uly_chimap_write, <chimap>, <outfile>
;
; ARGUMENTS:
;    <chimap>   Chimap structure created by ULY_CHIMAP
;    <outfile>  FITS file where the result is saved
;
; AUTHOR:
;              Antoine Bouchard
; HISTORY:
;              2008/02/08
;
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
pro uly_chimap_write, chimap, outfile

mkhdr, hdr, chimap.map
fxaddpar, hdr, 'CTYPE1', chimap.xtitle
if chimap.xunit ne '' then fxaddpar, hdr, 'CUNIT1', chimap.xunit
fxaddpar, hdr, 'CTYPE2', chimap.ytitle
if chimap.yunit ne '' then fxaddpar, hdr, 'CUNIT2', chimap.yunit
fxaddpar, hdr, 'EXTEND', 'T'
writefits, outfile, chimap.map, hdr
writefits, outfile, chimap.xnode, /append
writefits, outfile, chimap.ynode,  /append

end
;########################################################################
function range, x1, x2, n  ; return an array of n values sampling the range x1 to x2
  compile_opt idl2, hidden
  return, x1 + (x2-x1)*lindgen(n)/(n-1d)
end

function uly_chimap, spect, cmp, proj,                        $
                     VELSCALE=velscale, LMIN=lmin, LMAX=lmax, $
                     ERR_SP=err_sp,                           $
                     SG=sg, DG=dg, SNR=snr,                   $
                     KMOMENT=kmoment, KFIX=kfix,              $
                     XNODE=xnode, YNODE=ynode,                $
                     ADEGREE=adegree, MDEGREE=mdegree,        $
                     PLOT=plot, PS=ps,                        $
                     CLEAN=clean,                             $
                     QUIET=quiet, FILE_OUT=file_out


if n_elements(spect) eq 0 or n_elements(cmp) eq 0 then begin
    print,'Usage:'
    print,'uly_chimap, <spectrum: file or structure>, <cmp>, ...'
    return, 0
endif else if size(spect,/TYPE) ne size('str',/TYPE) then begin
    spectrum = spect
endif else begin
    testinp = file_test(spect)
    if testinp eq 0 then $
      message,'The file:'+spect+' does not exist. Give correct name for the input file.'
endelse

if n_elements(cmp) eq 0 or size(cmp, /TYPE) ne 8 $
  then message, 'You must give a component array!'

if n_elements(proj) eq 0 then proj = 0
if not array_equal(size(proj), size([[1,2],[3,4]])) then begin
    print, 'YOU DID NOT GIVE VALID PROJECTIONS, SELECT THE FIRST (X):'
    proj1 = uly_cmp_select(cmp)
    print, '... SELECT THE SECOND (Y) :'
    proj2 = uly_cmp_select(cmp)
    proj = [[proj1],[proj2]]
endif


;; -------------------------------------------------------------------------
;; read the observations
if n_elements(spectrum) ne 0 then begin     ; if we pass a spectrum struct
    SignalLog = uly_spect_extract(spectrum, ONED=0, STATUS=status)
    if status ne 0 then message, 'Could not extract 1D spectrum'
    SignalLog = uly_spect_logrebin(SignalLog, velscale, /OVER)

endif else begin                            ; if we pass a filename
    SignalLog = uly_spect_read(spect, lmin, lmax, VELSCALE=velscale, $
                               SG=shift_guess, ERR_SP=err_sp, $
                               QUIET=quiet $ 
                              )
    SignalLog = uly_spect_extract(SignalLog, ONED=0, STATUS=status, /OVER)
    if status ne 0 then message, 'Could not extract 1D spectrum'
    if SignalLog.sampling ne 1 then $
      SignalLog = uly_spect_logrebin(SignalLog, velscale, /OVER)
endelse

if not uly_spect_get(SignalLog, /VALID) then $
  message, 'Input spectrum is invalid'

; make the good pixels list (if it does not exist) ... needed for the
; tests below
if n_elements(*SignalLog.goodpix) gt 0 then gp = *SignalLog.goodpix $
else gp = lindgen((size(*SignalLog.data,/DIM))[0])

; compute the error if not read from the file (use SNR)
if n_elements(*SignalLog.err) eq 0 then begin
    if n_elements(snr) eq 0 then snr = 100
    mean_error = mean((*SignalLog.data)[gp], /NAN) / snr 
    if not finite(mean_error) then $
      message, 'Cannot compute the mean of the signal!'
    *SignalLog.err = (*SignalLog.data) * 0 + mean_error
endif

; Check if the errors are positive 
negerr = where((*SignalLog.err)[gp] le 0, cnt, COMPLEM=poserr)
if cnt eq n_elements((*SignalLog.err)[gp]) then $
  message, 'The noise is negative!!!'
if cnt gt 0 then $
  (*SignalLog.err)[gp[negerr]] = min((*SignalLog.err)[gp[poserr]])

; Check if some pixels have exagerated weight
weight = 1D / ((*SignalLog.err)[gp])^2
large_weight = where(weight gt 100*mean(weight), cnt)
if cnt gt 0 then message, /INFO, $
  'Some pixels have more than 100 times the average weight ... '+ $
  'it may be an error (up to ' + strtrim(max(weight)/mean(weight),2)+')'

;; -------------------------------------------------------------------------
;; read the models

lamrange = uly_spect_get(SignalLog, /WAVERANGE, STATUS=status)
if status ne 0 then message, 'Could not obtain WAVERANGE for input spectrum'

velscale = SignalLog.step * 299792.458d

; uly_fit_init read the models and place them in the cmp.eval_model
status = uly_fit_init(cmp, WAVERANGE=lamrange, VELSCALE=velscale, QUIET=quiet)
if status ne 0 then begin
    message, 'Fit initialization failed, abort', /CONT
    return, 0
endif

; The signal spectrum may be longer than the model, lets cut it
;    (have to cut the obs narrower in case the 2 grids are shifted)
model_range = exp(cmp[0].start + [0.5,cmp[0].npix-1.5] * cmp[0].step)
SignalLog = uly_spect_extract(SignalLog, WAVERANGE=model_range, /OVER)

if min(cmp.npix) gt (size(*SignalLog.data))[1]+1 then begin
    print, 'The number of pix in the model is greater than in Observation'
    print, '* THIS IS PROBABLY A BUG IN THE PROGRAM *'
    print, 'Please report it'
    print, cmp.npix, ' vs. ', (size(*SignalLog.data))[1]
endif
;; -------------------------------------------------------------------------
;; set some defaults
if not keyword_set(mdegree) then mdegree = 10
if not keyword_set(adegree) then adegree = -1
if n_elements(kmoment) eq 0 then kmoment = 2; moments of the gaussian to be fit

if n_elements(kfix) gt 0 then begin
    if n_elements(kfix) gt kmoment then $
      message, 'The number of elements of KFIX should not exceed '+$
      strtrim(string(kmoment),2)
endif; else kfix = intarr(kmoment)

if n_elements(dg) eq 0 then sigma_guess = velscale else sigma_guess = dg

kguess = [0, sigma_guess]

;; -------------------------------------------------------------------------
;; fix the correct parameters and set the guesses

xx = proj[*, 0] & yy = proj[*, 1]
xname = ''      & yname = ''
xunit = ''      & yunit = ''
xdispf = ''     & ydispf = ''

if xx[0] ne -1 then begin    ;if we want to project the cmp parameters
    if xx[1] ne -1 then begin
        fixed_x = (*cmp[xx[0]].para)[xx[1]].fixed   ; save the state
        guess_x = *(*cmp[xx[0]].para)[xx[1]].guess  ; save the state
        (*cmp[xx[0]].para)[xx[1]].fixed = 1
        if n_elements((*cmp[xx[0]].para)[xx[1]].name) gt 0 then begin
            xname  = (*cmp[xx[0]].para)[xx[1]].name
            xunit  = (*cmp[xx[0]].para)[xx[1]].unit
            xdispf = (*cmp[xx[0]].para)[xx[1]].dispf
        endif
; chose the values of the X axis on which you want to compute chi2
        if n_elements(xnode) le 1 then begin
            nnode = 10
            if n_elements(xnode) eq 1 then if xnode ge 2 then nnode = xnode
            if n_elements((*cmp[xx[0]].para)[xx[1]].limits) eq 2 then begin
                xgr0 = ((*cmp[xx[0]].para)[xx[1]].limits)[0]
                xgr1 = ((*cmp[xx[0]].para)[xx[1]].limits)[1]
                xgrid = xgr0 + dindgen(nnode)*(xgr1-xgr0)/(nnode-1)
                xgrid[0] += (xgrid[1]-xgrid[0])*0.1  
                xgrid[nnode-1] -= (xgrid[1]-xgrid[0])*0.1
            endif else message, 'There are no limits in the cmp, give XNODE'
        endif else xgrid = xnode
    endif else begin
        if n_elements(xnode) le 1 then $
          xgrid = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100] $
        else xgrid = xnode
        guess_x = cmp[xx[0]].lim_weig   ; save the state
    endelse
endif else begin           ;if we want to project the LOSVD parameters
    kfix[xx[1]] = 1
    nnode = 10
    if n_elements(xnode) eq 1 then if xnode ge 2 then nnode = xnode
    if n_elements(xnode) le 1 then xgrid = 10.+190.*dindgen(nnode)/(nnode-1) $
    else xgrid = xnode
    case xx[1] of
        0 : xname = 'cz'
        1 : xname = 'Velocity dispersion'
        else : message,'Not valid LOSVD value, we can fix only v and sigma'
    endcase
    xunit = 'km/s'
endelse

if yy[0] ne -1 then begin    ;if we want to project the cmp parameters
    if yy[1] ne -1 then begin
        fixed_y = (*cmp[yy[0]].para)[yy[1]].fixed   ; save the state
        guess_y = *(*cmp[yy[0]].para)[yy[1]].guess  ; save the state
        (*cmp[yy[0]].para)[yy[1]].fixed = 1
        if n_elements((*cmp[yy[0]].para)[yy[1]].name) gt 0 then begin
            yname = (*cmp[yy[0]].para)[yy[1]].name
            yunit = (*cmp[yy[0]].para)[yy[1]].unit
            ydispf = (*cmp[yy[0]].para)[yy[1]].dispf
        endif
; chose the values of the Y axis on which you want to compute chi2
        if n_elements(ynode) le 1 then begin
            nnode = 10
            if n_elements(ynode) eq 1 then if ynode ge 2 then nnode = ynode
            if n_elements((*cmp[yy[0]].para)[yy[1]].limits) eq 2 then begin
                ygr0 = ((*cmp[yy[0]].para)[yy[1]].limits)[0]
                ygr1 = ((*cmp[yy[0]].para)[yy[1]].limits)[1]
                ygrid = ygr0 + dindgen(nnode)*(ygr1-ygr0)/(nnode-1)
                ygrid[0] += (ygrid[1]-ygrid[0])*0.1  
                ygrid[nnode-1] -= (ygrid[1]-ygrid[0])*0.1
            endif else message, 'No limits in CMP found, give YNODES'
        endif else ygrid = ynode
    endif else begin
        if n_elements(ynode) le 1 then $
          ygrid = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100] $
        else ygrid = ynode
        guess_y = cmp[yy[0]].lim_weig   ; save the state
    endelse
endif else begin ;if we want to project the LOSVD parameters
    kfix[yy[1]] = 1
    nnode = 10
    if n_elements(ynode) eq 1 then if ynode ge 2 then nnode = ynode
    if n_elements(ynode) le 1 then ygrid = 10.+190.*dindgen(nnode)/(nnode-1) $
    else ygrid = ynode
    case yy[1] of
        0 : yname = 'cz'
        1 : yname = 'Velocity dispersion'
        else : message,'Not valid LOSVD value, we can fix only v and sigma'
    endcase
    yunit = 'km/s'
endelse

;; -------------------------------------------------------------------------
;; prepare the call of the fit

; The signal spectrum may be longer than the model, lets cut it
model_range = exp(cmp[0].start + [0,cmp[0].npix-1] * cmp[0].step)
SignalLog = uly_spect_extract(SignalLog, WAVERANGE=model_range, /OVER)

i_elem = n_elements(xgrid)
j_elem = n_elements(ygrid)
chimap = {xnode:xgrid, ynode:ygrid, map: dblarr(i_elem, j_elem), $
          xtitle:xname, ytitle:yname, $
          xunit:xunit, yunit:yunit}
undefine, xgrid & undefine, ygrid
if not keyword_set(quiet) then $
  print, 'Computing Chi square map on', n_elements(chimap.xnode), $
  ' X', n_elements(chimap.ynode), ' elements.'
totnode = n_elements(chimap.xnode)*n_elements(chimap.ynode)
node    = 0

;; -------------------------------------------------------------------------
;; call the fitting routine

for i = 0, i_elem - 1 do begin
    for j = 0, j_elem - 1 do begin
        if xx[0] ne -1 then begin ;if we want to project the cmp parameters
            if xx[1] eq -1 then  $
              cmp[xx[0]].lim_weig = [chimap.xnode[i], chimap.xnode[i]] $
            else $
              *(*cmp[xx[0]].para)[xx[1]].guess = chimap.xnode[i]
        endif else begin   ;if we want to project the LOSVD parameters
            kguess[xx[1]] = chimap.xnode[i]
        endelse
        
        if yy[0] ne -1 then begin ;if we want to project the cmp parameters
            if yy[1] eq -1 then  $
              cmp[yy[0]].lim_weig = [chimap.ynode[j], chimap.ynode[j]] $
            else $
              *(*cmp[yy[0]].para)[yy[1]].guess = chimap.ynode[j] 
        endif else begin   ;if we want to project the LOSVD parameters
            kguess[yy[1]] = chimap.ynode[j]
        endelse
        solution = uly_fit (SignalLog, CMP=cmp,                        $
                            KMOMENT=kmoment, KGUESS=kguess, KFIX=kfix, $
                            MDEGREE=mdegree, ADEGREE=adegree, /QUIET, $
                            CLEAN=clean)    
        chimap.map[i,j] = solution.chisq
        node = node+1

        if not keyword_set(quiet) then $
          print, $
          FO='(%"X = %7.1f, Y = %5.2f, chisq = %7.4f   (%5.1f%% done)")', $
          chimap.xnode[i], chimap.ynode[j], chimap.map[i,j], 100.*node/totnode
        
        if keyword_set(plot) then begin
            uly_chimap_plot, chimap, PS=ps, /QUIET
        endif
    endfor
endfor

; restore the initial state in cmp
if xx[0] ne -1 then begin 
    if xx[1] ne -1 then begin
        (*cmp[xx[0]].para)[xx[1]].fixed = fixed_x
        *(*cmp[xx[0]].para)[xx[1]].guess = guess_x 
    endif else cmp[xx[0]].lim_weig = guess_x
endif 
if yy[0] ne -1 then begin 
    if yy[1] ne -1 then begin
        (*cmp[yy[0]].para)[yy[1]].fixed = fixed_y
        *(*cmp[yy[0]].para)[yy[1]].guess = guess_y 
    endif else cmp[yy[0]].lim_weig = guess_y
endif


;; -------------------------------------------------------------------------
;; save, plot, etc...

if xdispf eq 'exp' then chimap.xnode = exp(chimap.xnode)
if ydispf eq 'exp' then chimap.ynode = exp(chimap.ynode)

if keyword_set(plot) then uly_chimap_plot, chimap, PS=ps, QUIET=quiet

if n_elements(file_out) gt 0 then uly_chimap_write, chimap, file_out

return, chimap

end
