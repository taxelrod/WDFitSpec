;+
; NAME:
;		ULY_FIT_LIN
; PURPOSE:
;		Fit the linear coefficients and return the best fit
;
; USAGE:
;    bestfit = uly_fit_lin(ADEGREE=adegree, 
;                      PAR_LOSVD=par_losvd, CMP=cmp,         
;                      SPECTRUM=SignalLog,                   
;                      GOODPIXELS=goodPixels, 
;                      VOFF=voff,                          
;                      MPOLY=mpoly, ADDCONT=addcont,         
;                      POLPEN=polpen [,/LUM_WEIGHT] [,/QUIET])
;
; DESCRIPTION:
;    Evaluate the different components of the model to fit and convolve
;    them with the LOSVD for a given set of the free parameters.
;    Determine the optimal weight of each component and the coefficients
;    of the multiplicative and additive continuum (linear fit).
;    Return the best fitting model.
; 
;    This function is called at each step of the Levenberg-Marquart
;    minimization performed in ULY_FIT to compute the values of the
;    non-linear parameters of the components and LOSVD.
;
;    For performance reasons, the large data arrays are passed through 
;    commons (not as arguments).
;
; KEYWORDS:
;               See also the documentation of uly_fit
;   ADEGREE (input)
;      Degree of additive polynomial, -1 if no additive polynomial
;   POLPEN (input)
;      Bias level of the multiplicative polynomial. This keyword can be 
;      used to reduce the effect of insignificant terms of the multiplicative
;      polynomial.
;      Terms whose coefficients are smaller, in absolute value, than <polpen> 
;      times their statistical error are amorted by the factor 
;      (abs(coef)/(<polpen>*err))^2.
;      This feature is active only if MDEGREE>0. <polpen>=2 seems a 
;      reasonable choice.
;      By default,  no bias is applied.
;   PAR_LOSVD (input)
;      Parameters of the LOSVD, array of 0 to 6 numbers ordered as: 
;      [cz, sigma, h3, h4, h5, h6], cz and sigma are in pixels.
;   SPECTRUM  (input)
;      Structure containing the spectrum to analyse and its associated error.
;      It is used to determine the weight of each pixel (.err) and to
;      fit the multiplicative continuum (.data)
;   GOODPIXELS (input)
;      List of pixels to be used in the fit.
;   VOFF (input)
;      Velocity offset between the spectrum to analyse and the model. km/s.
;   CMP       (Input & output)
;      Array of components to fit (see documentation in ULY_FIT)
;      In output CMP.WEIGHT is updated with the weight of each component 
;      computed by this routine.
;   MPOLY  (Input/output)
;      multiplicative continuum Legendre polynomial
;      structure {lmdegree, mpolcoefs, poly}
;      .LMDEGREE (integer), Input, maximum degree of legendre polynomials.
;         Unchanged in output.
;      .MPOLCOEFS (double 1D array), LMDEGREE+1 terms.
;         Output: updated values of the coefficients
;         The input values of the coefficients are ignored.
;      .POLY (1D array)
;         Input: polynomial determined at a previous iterration
;         Output: updated value of the polynomial (consistent with .MPOLCOEFS)
;      Usually MPOLY is initialized by ULY_FIT, but if it is NULL in input,
;      it is initialized to contain a degree 0 polynomial (scaling factor).
;   ADDCONT (output)
;      Additive polynomial (array).
;   /LUM_WEIGHT
;      Instruct to compute the luminosity weights and the errors on the
;      weights.
;   /QUIET
;      Set this keyword to run silently.
;
; AUTHOR:
;		Mina Koleva, Philippe Prugniel
;               with credit to Michele Cappellari & Igor Chilingarian
;               for initial development
; HISTORY:
;    The principle of solving the linear dependences within the Levenberg-Marquart
;    function evaluation is borrowed from M. Cappellari in his ppxf package.
;    v.1.0      From early development by M. Cappelleri & I. Chilingarian
;    v.2.0      Separated from the main program, Mina Koleva Sep 2005
;    v.2.3      Changing the interpolation (improve performance): mk 2007/03/01
;    v.3.0      Interpolating between different Mg/Fe: mk 2007/04/17
;    v.4.0      Ph. Prugniel 2007/12/14: linear fit of mult pol. Cache Legendre
;    v.5.0      Ph. Prugniel, 2008/05/27: general n-components fit
;    v.6        Ph. Prugniel, 2011/02: Implement Tikhonov regularisation
;    v.7        Ph. Prugniel, 2015/09: Major debugging and improvements
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
;
; uly_fit_fractsolve:
;  This routine determines the weight of the different 'components' of the 
;  model: one or several components and 'npoly' terms of the additive
;  polynomial.
;  'a' can be a vector (if there is a single template and no polynomial)
;  or an array whose first 'npoly' lines are the polynomial basis, and
;  the next line(s) the templates.
;  'b' is the 'observed' spectrum
;  bvls_state is a I/O array used to store the state of BVLS_PS, for
;  the warm start 
;
;  The function returns an array containing the weight of each component.
;  It is called by 'uly_fit_lin_weight'.

function uly_fit_fractsolve, a, b, npoly, BOUNDS=bounds, STATE=bvls_state
compile_opt idl2, hidden
s = size(a)

if s[0] eq 1 then begin                   ; Fit a single cmp, no add polynomial
    soluz = total(a*b,/DOUBLE)/total(a^2,/DOUBLE)
endif else if s[2] eq npoly+1 then begin ; Fitting a single cmp (no bound)
    soluz = la_least_squares(transpose(a), b, /DOUBLE)
endif else begin                         ; Fitting multiple templates (bounds)
    bnd = replicate(1d,2,s[2]) ; faster than dblarr
;    bnd = dblarr(2,s[2], /NOZERO)
    mx = (machar(/DOUBLE)).xmax
    if npoly gt 0 then begin
       bnd[0,0:npoly-1] = -mx ; No bounds on additive polynomials
       bnd[1,0:npoly-1] = mx  ; No bounds on additive polynomials
    endif
    bnd[*,npoly:*] = double(bounds)
    bvls_ps, a, b, bnd, soluz, bvls_state, ITMAX=15*s[2], STATUS=ierr
;    bvls_lh, a, b, bnd, soluz, ITMAX=15*s[2], STATUS=ierr
    if ierr ne 0 then message, 'BVLS returned with error #' + strtrim(ierr,2)
endelse

return, soluz
end


;-----------------------------------------------------------------------------
; Determination (or refining) of the multiplicative Legendre polynomial
function uly_fit_lin_mulpol, bestfit, CMP=cmp, MPOLY=mpoly,           $
                             SPECTRUM=SignalLog,                      $
                             GOODPIXELS=goodPixels,                   $
                             QUIET=quiet,                             $
                             STATUS=status

status = 0
npix = (size(*SignalLog.data, /DIM))[0]        ; number of wavelength bins
galaxy = double((*SignalLog.data)[goodPixels]) ; signal for non-masked points
noise = double((*SignalLog.err)[goodPixels])   ; error for non-masked points

if mpoly.lmdegree eq 0 then return, bestfit

; The polynomial may get crazy in the unconstrained regions, and this
; may affect the edges of these regions...
; It is critical when we want to fit spectral segment separated by a large
; range without any data.
; To prevent this, we may fit the polynomial everywhere, replacing the
; missing data with values in a correct range, and giving them low weight.
; like:
;    galtmp = galaxy*0 + 1e
;    galtmp[goodPixels] = galaxy[goodPixels]
;    meas_err = galaxy*0 + (mean(noise[goodPixels]) * 100)
;    meas_err[goodPixels] = noise[goodPixels]
; and then fit as:
;    coefs_pol = mregress(cmul, galtmp, MEASURE_ERRORS=meas_err, INVERSE=inverse, STATUS=status)
; This should work in most places, and replacing the '1' by an
; interpolation of the data should make this solution quite general.
; Though, for the moment, we prefer to keep the conservative solution.
; So, the problem still exists.

bestfitp = bestfit / mpoly.poly
cmul = mpoly.leg_array * rebin(bestfitp,n_elements(bestfitp),mpoly.lmdegree+1)

coefs_pol = mregress(cmul[goodPixels,*], galaxy, MEASURE_ERRORS=noise, STATUS=status_mregress)
if status_mregress eq 1 then status = 1
; note that la_least_squares or bvls_ps are about 50% slower than mregress
;  coefs_pol = la_least_squares(transpose(cmul[goodPixels,*]), galaxy, /DOUBLE, STATUS=status)
;  coefs_pol = transpose(coefs_pol)

if n_elements(polpen) ne 0 then begin
;  penalization of unsignificant terms in the pol.
   penal_pol = abs(coefs_pol)/(polpen*inverse.sigma)
   pen = 1+[where (penal_pol[1:*] lt 1)]
   coefs_pol[pen] *= penal_pol[pen]^2 
endif

mpoly.poly = mpoly.leg_array[*,0:mpoly.lmdegree] # transpose(coefs_pol)

minmaxpol = minmax(mpoly.poly[goodPixels])
if minmaxpol[0] lt 0 then begin
   lmd = mpoly.lmdegree 
   while min(mpoly.poly) le 0 and lmd gt 0 do begin
      lmd -= 1
      coefp = mregress(cmul[goodPixels,0:lmd], galaxy, $
                       MEASURE_ERRORS=noise,       $
                       STATUS=status_mregress)
      if status_mregress eq 1 then status = 1
      if n_elements(polpen) ne 0 then begin
         penal_pol = abs(coefp) / (polpen*inv.sigma)
         pen = 1+[where (penal_pol[1:*] lt 1)]
         coefs_pol[pen] *= penal_pol[pen]^2 
      endif
      mpoly.poly = (mpoly.leg_array[*,0:lmd] # transpose(coefp))
   endwhile
   
   if not keyword_set(quiet) then $
      print, 'WARNING: the polynomial became negative, min:'+strtrim(minmaxpol[0],2),' its degree is reduced to:', lmd
   
   coefs_pol = 0*coefs_pol + coefp
   
   if min(mpoly.poly) le 0 then if not keyword_set(quiet) then begin
      message, /INFO, 'polynomial is negative'
      status = 1
      return, 0
   endif
endif

bestfitp *= mpoly.poly 
if keyword_set(lum_weight) and n_elements(cmp.weight) gt 0 then begin
   cmp.e_weight /= coefs_pol[0]
   cmp.l_weight /= coefs_pol[0]
endif
mpoly.mpolcoefs = coefs_pol

return, bestfitp
end

;------------------------------------------------------------------------------
; Determination of the weight of each component
function uly_fit_lin_weight, SignalLog, goodpixels, tmp, mpoly, $
                             ADEGREE=adegree,                   $
                             LUM_WEIGHT=lum_weight,             $
                             CMP=cmp,                           $
                             ADDCONT=addcont,                   $
                             STATE=bvls_state,                  $
                             QUIET=quiet,                       $
                             STATUS=status
; Input arguments:
;   SignalLog  : observation (spect structure)
;   goodpixels : list of pixels to use 
;   tmp        : templates
;   mpoly      : structure containing the multiplicative polynomial
; Keywords (input)
;   adegree    : degree of additive polynomial
;   lum_weight : keyword switch: if set, compute e_weight and l_weight
; Keyword (output)
;   cmp        : result returned as cmp.weight (and .e_weight .l_weight if /LUM_WEIGHT)
;   addcont    : additive continuum (array)
;   status     : error status: 0= no error
; RETURN
;   bestfit    : array with the best matching model

n_cmp = n_elements(cmp)                        ; number of components
npix = (size(*SignalLog.data, /DIM))[0]        ; number of wavelength bins
n_cmp = n_elements(cmp)                        ; number of components
galaxy = double((*SignalLog.data)[goodPixels]) ; signal for non-masked points
noise = double((*SignalLog.err)[goodPixels])   ; error for non-masked points
status = 0

; Fill the columns of the design matrix of the least-squares problem
c = rebin(double(mpoly.poly), npix, n_cmp) * tmp  ; Fill with models*mulpol

if adegree ge 0 then begin   ; If there is an addpol, prepend addpol*mulpol
    c1 = rebin(double(mpoly.poly), npix,(adegree+1)) 
    for j=0,adegree do begin 
        x = 2d*dindgen(npix)/npix - 1d 
        c1[0,j] *= legendre(x,j)
    endfor
    c = [[c1],[c]]
endif

a = c[goodpixels,*]                            ; Used to solve the system
for j=0,(adegree+1)+n_cmp-1 do a[*,j] /= noise ; Weight with errors

;-----------------------------------------------------------------------------
; Experiment with Tikhonov regularization (2011/02/05)
; We add to the system above some linear constraints , weighted with "lambda"
; See the chapters 18.4 and 18.5 of Numerical Recipes
lambda = 0 ; regularization weight : 0 => no regularization
;lambda = 1d-13 * double(npix)/( n_cmp-1)
;
; For the moment, we experimented with 1st order regularization, i.e.
; the additional constraints asymptotically lead to a solution where all
; the weights are equal.
; The solution that we tried should still be finalized: The regularization
; weight should be scaled more intelligently. lambda should be insensitive
; to: (i) the number of wavelength bins (or VELSCALE), (ii) the flux scale
; in the data and (iii) the flux scale in the models. 5This should not be
; an issue to implement this
; We shall also explore these other aspects:
;   - Other types of regularization, like 2nd order (linear dependence)
;     or regularization matrix passed by the user (we may provide a different
;     routine to help the creation of this matrix)
;   - We may want to regularize on the l_weight rather than on the weights.
; We have to find a compromize between the flexibility of the program and
; its ease of use. The regularization is a complex business: How can the
; user choose the proper value of lambda, and the type of regularization?
; If we just give a single regularization (say 1st order) and a default
; lambda ... the user may be happy to have something functional at low
; expense of brain. But it is not fully correct.
if lambda gt 0 then begin;
    ng = n_elements(goodpixels)
    nrel = ng + n_cmp - 1
    aa = dblarr(nrel, n_cmp)
    for k = 0, n_cmp-1 do aa[0,k] = a[*, k]
    a = temporary(aa)
    cr = dblarr((size(c, /DIM))[1])
    cr[adegree+1:adegree+2] = lambda * [-1, 1] ; 1st order regularization
    for k=ng, nrel-1 do begin
        a[k,*] = shift(cr, k-ng)
    endfor
    goodpixels = [goodpixels, lindgen(n_cmp-1)+npix]
    galaxy = [galaxy,fltarr(n_cmp-1)]
    noise = [noise,fltarr(n_cmp-1)+1]
endif
   
npoly = (adegree+1)   ; Number of additive polynomials in the fit

status_math = check_math(MASK=32)    ; galaxy/noise may cause an underflow: handle
if(adegree eq -1) then begin
    cmp.weight = uly_fit_fractsolve(a, galaxy/noise, npoly, BOUNDS=cmp.lim_weig, STATE=bvls_state)
    bestfit = c # cmp.weight
    addcont = dblarr(n_elements(bestfit)) ; needed to return arg addcont
endif else begin
    wght = uly_fit_fractsolve(a, galaxy/noise, npoly, BOUNDS=cmp.lim_weig, STATE=bvls_state)
    bestfit = c # wght
    cmp.weight = wght[adegree+1:*]
    addcont = (c[*,0:adegree] # wght[0:adegree])
endelse
if status_math eq 0 then status_math = check_math(MASK=32)

nan = where(finite(cmp.weight) eq 0 or abs(cmp.weight) gt 1d36, cnt)
if cnt gt 0 then begin
    message, /INFO, 'Could not determine the weight of cmp ' + $
             strjoin(strtrim(nan,2),',')+' ... assume 0'
    cmp[nan].weight = 0
endif
zeros = where(abs(cmp.weight) lt 1d-36, cnt)
if cnt gt 0 then cmp[zeros].weight = 0d
if abs(total(cmp.weight)) lt 1d-36 then begin
    if not keyword_set(quiet) then $
      message, /INFO, 'The total weight is null'
    status = 1
    common MPFIT_ERROR, error_code
    error_code = -2 ; tell MPFIT to stop, the best fit is a null spectrum
    return, 0
endif

if keyword_set(lum_weight) and n_elements(cmp.weight) gt 0 then begin
;   we call this only when we want to evaluate the errors on the
;   weights and the luminosity weights
    w_notnull = where(cmp.weight ne 0, nw_notnull)
    cmp.e_weight = dblarr(n_cmp)
    cmp.l_weight = dblarr(n_cmp)      
    if nw_notnull ne 0 then begin
        inv = 0
        if adegree gt -1 then $
           w_nn = [indgen(adegree+1), w_notnull+adegree+1] $
        else $
           w_nn =  w_notnull
        ww = mregress((c[goodPixels,*])[*,w_nn], galaxy, $
                      MEASURE_ERRORS=noise, INVERSE=inv, STATUS=status)
        cmp[w_notnull].e_weight = inv.sigma[adegree+1:*]
        cmp.l_weight = (total(tmp, 1) * cmp.weight) / npix
    endif
endif else begin
    cmp.e_weight = 0
    cmp.l_weight = 1
endelse

return, bestfit

end
;------------------------------------------------------------------------------
FUNCTION uly_fit_lin, ADEGREE=adegree,            $
                      PAR_LOSVD=par_losvd,        $
                      CMP=cmp,                    $
                      SPECTRUM=SignalLog,         $
                      GOODPIXELS=goodPixels,      $
                      VOFF=voff,                  $
                      MPOLY=mpoly,                $
                      ADDCONT=addcont,            $
                      LUM_WEIGHT=lum_weight,      $
                      POLPEN=polpen,              $
                      CACHE=cache,                $
                      STATE=bvls_state,           $
                      MODECVG=modecvg,            $
                      QUIET=quiet,                $
                      STATUS=status

compile_opt idl2

status = 0

; Clear the accumulated stack of math error, in order to diagnos errors
; from the current routine.
if !EXCEPT gt 0 then status_math = check_math(/PRINT) else $
  if check_math() then print, $
  'An arithmetic error occured before calling uly_fit_lin'

npix = (size(*SignalLog.data, /DIM))[0] ; number of wavelength bins
n_cmp = n_elements(cmp)                 ; number of components
npar_losvd = n_elements(par_losvd)      ; number of free params for LOSVD

galaxy = double((*SignalLog.data)[goodPixels]) ; signal for non-masked points
noise = double((*SignalLog.err)[goodPixels])   ; error for non-masked points

if n_elements(mpoly) eq 0 then begin ; default poly initialization
    mpoly = {lmdegree:0, mpolcoefs:dblarr(1), poly:dblarr(npix)}
    mpoly.poly = 1d
endif

; check that noise is strictly positive, at least on the place of goodpixels
if min(noise) le 0 then message,  'NOISE vector must be strictly positive'

;-----------------------------------------------------------------------------
; evaluate all the components
if size(cache,/TYPE) ne 8 then begin
   cache = {para:ptrarr(n_elements(cmp), /ALLOCATE_HEAP), spectra:dblarr(npix, n_cmp)}
   for i=0, n_elements(cmp)-1 do *(cache.para)[i] = -1D35 
endif

;init. models spectra
models = dblarr(npix,n_cmp,/NOZERO)

for i=0, n_cmp-1 do begin
    cached = 0
    if n_elements(*cmp[i].para) gt 0 then begin
        if min((*cmp[i].para).value eq *(cache.para)[i]) eq 1 then begin
            models[0,i] = cache.spectra[*,i] 
            cached = 1
        endif
    endif

    if cached eq 0 then begin
;       We treat separately SSPs and stars because it is (a bit) faster
;       than doing everything with call_function
        case cmp[i].eval_fun of
            'SSP' : mdl = $
              uly_ssp_interp(*cmp[i].eval_data, (*cmp[i].para).value)
            'STAR' : mdl = *cmp[i].eval_data
            else : mdl = $
              call_function(cmp[i].eval_fun, *cmp[i].eval_data, $
                            (*cmp[i].para).value)
        endcase
        nmd = min([npix, (size(mdl,/DIM))[0]]) - 1
        models[0:nmd,i] = mdl[0:nmd]
        if n_elements(*cmp[i].para) gt 0 then $
           *(cache.para)[i] = (*cmp[i].para).value
        cache.spectra[*,i] = models[*,i]
    endif

endfor

if n_elements(par_losvd) gt 0 then begin
;-----------------------------------------------------------------------------
; generate the LOSVD kernel and convolve
; Detailed tests have shown that this kernel is correct as long
; as sigma > 0.8 (it is still acceptable down to sig=0.6)
; For low values of sigma:
;   - The total flux of the kernel is too large (need a normalization)
;   - The function is too much picked (for sigma around 0.3)
;     when the kernel is centered and too flat when the kernel is
;     shifted by half a pixel
; This low-sigma bias of the kernel can be modelled as:
;     sigmac = sqrt(sigma^2 + 0.125*exp(-(sigma-0.28)^2/(0.04)))
; Where sigma is the value used in the program and sigmac the de-biased 
; value corresponding to the real broadening. For the lack of simplicity,
; this correction is not applied in the program, but the user should
; be aware of the bias. A healthy solution consist in rebinning the
; spectrum to a scale where sigma>0.6

    maxvel = abs(par_losvd[0])
    maxsig = par_losvd[1]
    dx = ceil(abs(voff)+abs(maxvel)+5d*maxsig) ; Sample the LOSVD at least to vel+5*sigma
    dx = min([dx, (npix-1.)/2.])
    n = 2*dx + 1
    x = dx - findgen(n) 
    losvd = dblarr(n,/NOZERO)
            
    vel = voff + par_losvd[0]
    sigma_pix = (par_losvd[1] gt 0.1) ? par_losvd[1] : 0.1d
    w = (x - vel)/sigma_pix
    w2 = w^2
;   Replace the values larger then 5*sigma with zeros
    wlarge = where(abs(w) gt 5., cnt, COMPLEMENT=wnorm)
    if cnt gt 0 then losvd[wlarge] = 0. 
    losvd[wnorm] = exp(-0.5d*w2[wnorm])/(sqrt(2d*!dpi)*sigma_pix) 
        
; Hermite polynomials normalized as in Appendix A of van der Marel &
; Franx (1993). They are given e.g. in Appendix C of Cappellari et al. (2002)
    nherm = (size(par_losvd,/DIM))[0] - 2
    if nherm gt 0 then begin
        poly = 1d + par_losvd[2]/sqrt(3d)*(w*(2d*w2-3d)) ; h3
        if nherm gt 1 then $
          poly += par_losvd[3]/sqrt(24d)*(w2*(4d*w2-12d)+3d) ; h4
        if nherm gt 2 then $
          poly += par_losvd[4]/sqrt(60d)*(w*(w2*(4d*w2-20d)+15d)) ; h5
        if nherm gt 3 then $    ; h6
          poly += par_losvd[5]/sqrt(720d)*(w2*(w2*(8d*w2-60d)+90d)-15d)
        
        losvd *= poly
    endif   
    losvd /= total(losvd) ; Normalized total(Gaussian)=1 (need for sigma<0.6)

;   convolution with the LOSVD 
    status_math = check_math(MASK=32) ; convolution may cause an underflow: handle
    tmp = convol(temporary(models),losvd,/EDGE_TRUNCATE) 
    if status_math eq 0 then status_math = check_math(MASK=32)

endif else tmp = temporary(models)

; Iterative loop to get the weight and mulcont
;   when n_cmp>1 and md>0 the determination of the weights
;   and of mulpol is not a linear process, and at present we
;   need iteration to converge... and unfortunately, our
;   current approach is slow (we tried to accelarete it but
;   we failed for the moment to obtain a stable solution).
;   However, in most cases, the LM iteration will converge
;   even if uly_fit_lin does not converge at each step.
; So, we are currently offering three options:
;   MODECVG=0 (the default one), this is the fastest, but 
;     it happens that the solution is missed
;   MODECVG=1 (the convergence is achieved only once for each
;     LM iteration, it is not reached for the derivatives
;   MODECVG=2, uly_fit_lin always converge, but this is slow
; We are actively working at solving the problem and make 
; something fast and reliable.
itermax = 0
if n_elements(modecvg) eq 0 then modecvg=0
if modecvg eq 1 or modecvg eq 2 then itermax=499

fstop = 0
niter = 0
while fstop eq 0 do begin
;  Determination of the weight of each component
   bestfitw = uly_fit_lin_weight(SignalLog, goodpixels, tmp, mpoly, $
                                 ADEGREE=adegree,                   $
                                 LUM_WEIGHT=lum_weight,             $
                                 CMP=cmp,                           $
                                 ADDCONT=addcont,                   $
                                 STATE=bvls_state,                  $
                                 QUIET=quiet,                       $
                                 STATUS=status)
   if status ne 0 then return, bestfitw

   posw = where(cmp.weight gt 0, cposw)
   if cposw eq 0 then begin
       if not keyword_set(quiet) then $
         message, /INFO, 'All the cmp have a null weight... cannot continue'
       common MPFIT_ERROR, error_code
       error_code = -2          ; tell MPFIT to stop
       status = 1
       return, 0
   endif

; Determination of the multiplicative polynomials
   bestfit = uly_fit_lin_mulpol(bestfitw, CMP=cmp, MPOLY=mpoly,    $
                                SPECTRUM=SignalLog,                $
                                GOODPIXELS=goodPixels,             $
                                QUIET=quiet,                       $
                                STATUS=status)
   if status ne 0 then return, float(bestfit)

   mpoly.poly /= mpoly.mpolcoefs[0]
   cmp.weight *= mpoly.mpolcoefs[0]
   mpoly.mpolcoefs /= mpoly.mpolcoefs[0]
   if cposw eq 1 then break ; no need to iterate when only on cmp is active

   niter += 1
   if abs(1d0-total(bestfit,/NAN)/total(bestfitw,/NAN)) lt 5d-9 or niter gt itermax then fstop = 1
;   if abs(mean(1d0-bestfit/bestfitw)) lt 5d-9 or niter gt itermax then fstop = 1
endwhile

bestfit = float(bestfit) 

; The test below is useful to identify that an error occured during
; the execution of uly_fit_lin.
; Note that we cannot use a CATCH handler, because it does not catch
; math errors !!
; Note that to find the actual statement causing an error, the
; variable !EXCEPT can be set to 3 in the IDL session.
if check_math(/PRINT) ne 0 then begin 
    message, /INFO, 'an arithmetic error occured'
    common MPFIT_ERROR, error_code
    error_code = -2 ; tell MPFIT to stop
    status = 1
    return, 0
endif


return, bestfit
end
