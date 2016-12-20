;+
; NAME:
;               ULY_FIT
;
; PURPOSE:
;               Fit a spectrum with a non-linear model
;
; USAGE:
;     solution = uly_fit(signalLog,                                           
;                        CMP=cmp,                                             
;                        KMOMENT=kmoment, KGUESS=kguess, KFIX=kfix, KPEN=kpen, KLIM=klim,
;                        ADEGREE=adegree,                                     
;                        MDEGREE=mdegree, POLPEN=polpen,                      
;                        /CLEAN, /QUIET,                                 
;                        MODECVG=modecvg,
;                        /ALL_SOLUTIONS)
;
; ARGUMENTS:
;   signalLog: structure containing the spectrum of the object to be
;       fitted. The models (see CMP) and the object have to be logarithmically
;       rebinned. The WCS of the spectrum is read in the structure, whose tags
;       are described in the documentation of ULY_SPECT_ALLOC.
;
; KEYWORDS:
;
;   CMP
;       Array describing the components of the model to fit. See the
;       description below.
;
;   KMOMENT: Order of the Gauss-Hermite moments to fit. If this keyword is 
;       set to 2 then only [cz,sigma] are fitted. 
;       Set this keyword to 4 to fit [h3,h4] and to 6 to fit [h3,h4,h5,h6]. 
;       Note that in all cases the Gauss-Hermite moments are fitted 
;       (nonlinearly) *together* with [cz,sigma]. 
;       If KMOMENT=0 then only the continuum polynomial, the weights of the components
;       and their coefficients are fitted.
;
;   KGUESS: vector (up to 6 elements)
;       [velStart,sigmaStart,h3,h4,h5,h6] with the initial estimate for the 
;       cz and the cz dispersion in km/s, and of the Gauss-Hermite coefficients.
;   
;   KLIM: 2D floating array. The first dimension are the low/high bounds. and the second
;       the parameter of the LOSVD (cz, sigma, h3 ,...)
;       The limits on cz and sigma are in km/s.
;
;   KFIX: vector of  zeros or ones (corresponds to "KGUESS"). 1 means
;       that the corresponding LOSVD parameter is fixed during the minimization
;       procedure. For example, to fix the velocity dispersion, specify
;       kfix=[0,1]
;
;   KPEN: This parameter biases the (h3,h4,...) measurements towards zero
;       (Gaussian LOSVD) unless their inclusion significantly decreases the
;       error in the fit. By default KPEN=0, meaning that the penalization is
;       disabled. With a strictly positive value, the solution (including cz 
;       and sigma]) will be less noisy. When using penalization, it is adivised
;       to test the choice of KPEN with Monte-Carlo simulations. A value in the
;       range 0.5 to 1.0 shall be a good starting guess. This keyword 
;       corresponds to the parameter lambda in Cappellari & Emsellem (2004, 
;       PASP 116, 138, Eqn. 9). 
;       Note that the penalty depends on the *relative* change of the fit 
;       residuals, so it is insensitive to proper scaling of the noise. A 
;       nonzero KPEN can be safely used even without a reliable error spectrum,
;       or with equal weighting for all pixels.
;       The implementation of penalization was taken from the ppxf program
;       by M. Cappellari.
;
;   ADEGREE: degree of the *additive* Legendre polynomial used to correct the
;       template continuum shape during the fit. Set ADEGREE = -1 not to 
;       include any additive polynomial. By default, no additive term is used.
;       
;
;   MDEGREE: degree of the *multiplicative* Legendre polynomial 
;       used to correct the continuum shape during the fit (default: 10). The
;       zero degree multiplicative polynomial is always included in the fit as
;       it corresponds to the weights assigned to the components.
;
;   POLPEN
;       Parameter to select only the significant terms of the multiplicative
;       polynomial. See ULY_FIT_LIN for details.
;
;   /CLEAN: set this keyword to iteratively detect and clip the outiers.
;       See below the description of the clipping algorithm.
;   
;   /QUIET
;       Set this keyword to suppress messages printed on the screen.
;
;   /ALL_SOLUTIONS: If CMP contains more than a single guess for each free 
;       parmeter (ie. a vector of guesses for at least parameter), all the 
;       'local' solutions are returned in the form of an array.
;       Otherwise, just the best solution is returned.
;   MODECVG: 0,1 or 2
;       when n_cmp>1 and md>0 the determination of the weights and of mulpol is
;       not a linear process, and at present we need iteration to converge... 
;       and unfortunately, our current approach is slow (we tried to accelerate
;       it but we failed for the moment to obtain a stable solution).
;       However, in most cases, the LM iteration will converge even if 
;       ULY_FIT_LIN does not converge at each step.
;       So, we are currently offering three options:
;        MODECVG=0 (the default one), this is the fastest, but 
;          it happens that the solution is missed
;        MODECVG=1 (the convergence is achieved only once for each
;          LM iteration, it is not reached for the derivatives
;        MODECVG=2, uly_fit_lin always converge, but this is slow
;       We are actively working at solving the problem and make something fast 
;       and reliable.
;
; OUTPUT:
;   solution:
;     Output structure containing the solution of the fit. It can
;     be printed/plotted by ULY_SOLUT_PRINT, ULY_SOLUT_TPLOT,
;     ULY_SOLUT_TWRITE, ... Its detailed content is described later
;     in this document.
;
;     The solution structure may be passed to ULY_SOLUT_TPRINT, 
;     ULY_SOLUT_PLOT, ULY_SOLUT_TWRITE, ...
;
;
; DESCRIPTION:
;   Fit a spectrum with a linear combination of non-linear components, 
;   convolved with a LOSVD. Multiplicative and additive Legendre
;   polynomials can also be included in the model.
;
;   The usage of ULY_FIT requires to specify the characteristics of the
;   LOSVD function (number of free parameters), the degrees of the polynomials,
;   and the 'description' of the model to fit.
;
;   The model to fit is described by the CMP argument. Since the model
;   is a linear combination of components, CMP is an array whose individual
;   element is a structure describing one of these components.
;   See the routine ULY_SSP for an example of how to define CMP.
;   Other types of components may be defined by the user (need some 
;   programming).
;
;   The CMP array may contains vectors of guesses for some parameters.
;   In this case, ULY_FIT will perform a local minimization for each
;   possible combination of these guesses and return the best solution
;   or the whole array of solutions if the keyword /ALL_SOLUTIONS is
;   given.
;
;   FITTING PROCEDURE:
;   The 'non-linear' parameters of the components of the model are determined
;   with the Levenberg-Marquart (LM) routine MPFIT implemented by C. Markwardt
;   and the other parameters (weight of each components, coefficients of 
;   the multiplicative and additive polynomials) are fitted linearly
;   at each iteration of the LM process.
;
;   LINE-OF-SIGHT VELOCITY DISTRIBUTION (LOSVD):
;   The line-of-sight velocity distribution (LOSVD) is described by a
;   Gaussian plus Gauss-Hermite expansion (V, sigma, h3, h4, h5, h6).
;   The fit is performed in pixels space using the Penalized Pixel
;   Fitting method described in Cappellari M., Emsellem E., 2004, 
;   PASP, 116, 138.
;
;   CLIPPING ALGORITHM:
;   The goal is to exclude the outliers of the fit when the /CLEAN keyword
;   is specified. This reduces the sensitivity of a fit to spikes or
;   undesirable features due to an unperfect observation, treatment or 
;   model.
;   The principle is to re-iterate the fit after the outliers are 
;   identified from the residuals and masked.
;   After a fit, the most descrepant pixels are identified using a
;   kappa x sigma detection limit (sigma is the standard deviation
;   of the residuals, and kappa a numerical factor). Kappa is chosen
;   as the lowest of 3,4,5 or 7 in order to find at maximum 3% of
;   outliers. An additional condition is that an outlier cannot be
;   explained by a small mismatch in wavelength (this is achieved by
;   comparing the residuals with the gradients of the model; this is 
;   particularly important if the model includes emission lines).
;   The immediate neighbours of the outliers are also masked if they
;   depart by 2 x sigma (this will for example detect some residuals of
;   sky subtraction). Finally, the features detected as outliers are
;   tracked until their absolute value decreases is above 1 sigma. This
;   achieve a clean removal of emission lines that are not in the model.
;   Up to 5 cleaning iterations are made if the algorithm does not converge.
;
; DESCRIPTION OF A COMPONENT OF THE MODEL TO FIT:
;   This section of the documentation will be mostly useful to users
;   intending to write their own model component. Some usual components
;   are included in the package, and to use them the following
;   information is not needed.
;
;   The ULY_FIT program fit a linear combination of non-linear components
;   convolved with a LOSVD and multiplied by a polynomial.
;   Each component is described by a CMP structure, and therefore the
;   complete fit is described by an array of CMP.
;
;   See the functions ULY_SSP or ULY_TGM for practical examples
;   of definition of a CMP.
;
;   The CMP structure contains all the information about a component,
;   like in particular the name of the functions to use to initialize
;   and to evaluate the corresponding spectrum and the list and
;   characteristics of its parameters. It is also used to return the
;   results of the fit...
;
;   The CMP array may be passed to ULYSS, ULY_FIT or other routines
; 
;   An individual component, element of CMP, is a structure made as:
;       .name     : Name of the component (string)
;                   This name must be unique for a given fit, the concatenation
;                   of .name and .para.name must identify a free parameter.
;       .descr    : Description of the component (string)
;                   Essentially for display and information purpose.
;       .init_fun : Name of the initialization routine (string)
;                   The initialization routine has the following template:
;                   cmp = <init_fun>(cmp, WAVERANGE=wave, VELSCALE=vels,
;                                    QUIET=quiet)
;                   where the input arguments are cmp, a single component
;                   WAVERANGE=dblarr(2), the range of wavelength of the
;                   generated model (must be the same as the observation),
;                   VELSCALE the velocity scale (size of the pixel) in km/s,
;                   and QUIET the verbosity control.
;                   <init_fun> is called by the routine ULY_FIT_INIT
;                   when the fit is performed.
;                   <init_fun> normally creates .eval_data described below.
;                   If no initialization is required this tag may be kept
;                   as a null string.
;       .init_data: data passed to init_fun (pointer)
;       .eval_fun : name of the routine to evaluate the component (string)
;                   The evaluation function has the template:
;                   spectrum = <eval_fun>(*cmp.eval_data, (*cmp.para).value)
;                   Where 'spectrum' the the array containing the spectrum 
;                   of the component of the model.
;                   *cmp.eval_data and (*cmp.para).value are described below.
;                   For SSP models (computed with ULY_SSP_INTERP,
;                   eval_fun='SSP' may be used instead of 
;                   eval_fun='uly_ssp_interp' (slightly faster).
;                   If the model can be evaluated as: spectrum=*cmp.eval_data
;                   i. e. if there is no free parameter, eval_fun='STAR'
;                   may be used.
;       .eval_data: data passed to eval_fun (pointer)
;                   Usually generated by the initialization routine
;                   like for example a grid of models convolved by the
;                   proper LSF that eval_fun will interpolate.
;       .para     : description of the parameters (array of pointers)
;                   Each free parameter of the model is described by the
;                   structure whose members are:
;            .name      : name of the parameter unique for a single component
;            .unit      : unit of the parameter (display purpose)
;            .guess     : (double) Starting value of the fit
;            .fixed     : 0 if the parameter is free, 1 if it is fixed
;            .limited   : 0 if no limit is set, 1 if limits are set
;            .limits    : dble array [low, high] of bounds on the parameters
;            .step      : (double) step for numerical derivation
;            .value     : (double) final value of the parameter
;            .error     : (double) uncertainty on value
;            .dispf     : (string) indicate how to display the value
;                         Empty if the parameters is directly printed
;                         or 'exp' to print exp(value) and compute error
;                         as err*exp(value)
;       .start   : WCS of the component, starting wavelength
;       .step    : WCS of the component, step 
;       .npix    : WCS of the component, number of pixels
;       .sampling: WCS of the component, sampling mode (0:linear, 1:log, 2:unused)
;       .mask    : mask of pixels (array) good=1, bad=0
;       .weight  : (double) Weight of the component, result of the fit
;       .e_weight: (double) Error on .weight
;       .l_weight: (double) light fraction of the weight
;       .lim_weig: dble array [low, high] of bounds on the weights
;
; DESCRIPTION OF THE SOLUTION STRUCTURE:
;   The solution structure returned by ULY_FIT has the following members:
;      .TITLE      : Description of the fit (name of the fitted file)
;      .HDR        : header as read from the analysed spectrum
;      .START      : WCS of the spectrum, starting wavelength
;      .STEP       : WCS of the spectrum, step in wavelength
;      .SAMPLING   : WCS 0:Linear/1:log/2:irregular
;      .DATA       : The input spectrum which is fitted
;      .ERR        : The error spectrum 
;      .BESTFIT    : spect structure containing the best fitting model
;                    including the LOSVD convolution and add/mult polynomials
;      .MASK       : Binary mask: 1: Fitted pixels, 0: Rejected pixels
;      .DOF_FACTOR : Ratio (number of pixels/number of independent pixels)
;      .CHISQ      : Reduced chi2
;      .ADDCONT    : Array with additive continuum
;      .MULCONT    : Array with multiplicative continuum
;      .ADDCOEF    : Coefficients of the additive polynomial (if used)
;      .MULCOEF    : Coefficients of the multiplicative polynomial (if used)
;      .LOSVD      : Fitted parameters of the LOSVD
;      .E_LOSVD    : Errors on the fitted LOSVD parameters
;      .CMP        : Array describing the components of the model,
;                    the fitted value of the parameters and their errors,
;                    see above the description of the CMP structure.
;
; REQUIRED ROUTINES:
;     MPFIT: by C.B. Markwardt from http://astrog.physics.wisc.edu/~craigm/idl/
;     ROBUST_SIGMA: by H. Freudenreich from http://idlastro.gfsc.nasa.gov/
;
; HISTORY AND CREDITS:
;    This program is part of the ULySS package, currently under development
;    by ULYSS team (Ph. Prugniel, M. Koleva, A. Bouchard & Y. Wu) 
;    ulyss@obs.univ-lyon1.fr
;
;    The developments started from the ppxf program from Cappellari
;    (http://www.strw.leidenuniv.nl/~mcappell/idl/ ; Cappellari &
;    Emsellem (2004, PASP, 116, 138) ), based
;    on  the Levenberg-Marquart routine, MPFIT from Craig Markwardt
;    (http://astrog.physics.wisc.edu/~craigm/idl/fitting.html)
; 
;    The parametric minimisation on age/metallicity was started by
;    I. Chilingarian, Obs. de Lyon & N. Bavouzet, Obs. de Paris.
;    
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
function uly_makeparinfo, cmp

compile_opt idl2, hidden

n_cmp = n_elements(cmp)

n_par = 0
for i=0,n_cmp-1 do n_par += n_elements(*cmp[i].para)
if n_par eq 0 then return, 0

pinf = REPLICATE({value:0D, step:1D-2,limits:[0D,0D],limited:[1B,1B],fixed:0B}, n_par)
n = 0
for i=0,n_cmp-1 do begin
    for j = 0, n_elements(*cmp[i].para)-1 do begin
        pinf[n].value = (*(*cmp[i].para)[j].guess)[0]
        pinf[n].step = (*cmp[i].para)[j].step
        pinf[n].limits = (*cmp[i].para)[j].limits
        pinf[n].limited = (*cmp[i].para)[j].limited
        pinf[n].fixed = (*cmp[i].para)[j].fixed
        n ++
    endfor
endfor

return, pinf

end


;----------------------------------------------------------------------------
; parse the content of 'pars' into par_losvd cmp
; cmp is initialized from its value in the common
pro uly_fit_pparse, pars, par_losvd, cmp,        $
                    KMOMENT=kmoment, QUIET=quiet, ERROR=error

compile_opt idl2, hidden

; The first kmoment elements of the pars array are the parameters of LOSVD
; The parameters of one LOSVD are: cz, sigma, h3, h4 (cz and sigma in pixels)
; If h3 and/or h4 are not determined, they are absent and the dimension
; of the array is smaller
if kmoment gt 0 then par_losvd = reform(pars[0:kmoment-1], kmoment)

; the components used in the fit are described in a cmp structure
; containing
;    - The name of the function to evaluate it
;    - Some frozen parameters (like the type of models, the id of the
;    star or whatever)
;    - The array of free parameters
; The evaluation function is called through CALL_FUNCTION

common uly_functargs_common, mpolyc, cmpc, cache, bvls_state, modecvg
cmp = cmpc  ; return in argument the cmp stored in the common

; load the parameters and data in the comp struct
i0 = kmoment
for i=0,n_elements(cmp)-1 do begin
    i1 = i0 + n_elements(*cmp[i].para) - 1
    if i1 ge i0 then (*cmp[i].para).value = pars[i0:i1]
    if n_elements(error) gt 0 then $
      if i1 ge i0 then (*cmp[i].para).error = error[i0:i1]
    i0 = i1 + 1
endfor

end

;------------------------------------------------------------------------------
pro uly_fitfunc_init, SPECTRUM=spec, CMP=cmp, MDEGREE=mdegree, MODECVG=modecvg

compile_opt idl2, hidden

; Declaration of commons used to cache some information used by uly_fit_lin.
; The purpose is to avoid to redo some previous calculation, and therefore 
; to improve performance.
common uly_functargs_common, mpoly, cmpc, cache, bvls_state, modecvgc

if (n_elements(cmp) eq 0) then message, 'Keyword CMP must be specified'
if (n_elements(spec) eq 0) then message, 'Keyword SPECTRUM must be specified'
if (n_elements(mdegree) eq 0) then message, 'Keyword MDEGREE must be specified'
npix = n_elements(*spec.data)

; initialize/re-nitialize mpoly
if n_elements(mpoly) eq 1 then $
   if mpoly.lmdegree ne mdegree or (size(mpoly.poly,/DIM))[0] ne npix then undefine, mpoly

if n_elements(mpoly) eq 0 then begin
   mpoly = {lmdegree:mdegree, mpolcoefs:dblarr(mdegree+1), poly:dblarr(npix),  $
            leg_array:dblarr(npix,mdegree+1)}
   x = 2d*dindgen(npix)/npix - 1d ; X in [-1,1] for Legendre Polynomials
   mpoly.leg_array = double(flegendre(X, mdegree+1))
endif

mpoly.mpolcoefs = 0d
mpoly.mpolcoefs[0] = 1d
mpoly.poly = 1d

cmpc = cmp  ; store cmp in the common

; initializing cache...
;   The cache allows to store the result of evaluation for the
;   case it is requested again.
;   Each component is cached separately.
;   The cached item is identified by the parameters of the component
;   cmp[i].para, stored in cache_para[i]
if size(cache,/TYPE) eq 8 then heap_free, cache
undefine, cache

; state of each cmp, kept for BVLS. undefined => cold start
undefine, bvls_state

if n_elements(modecvg) eq 0 then modecvg = 0
modecvgc = modecvg

END

;----------------------------------------------------------------------------
; MK 2006/08/10 
; ULY_FITFUNC is the user function used by mpfit.
; See the documentation of MPFIT for general information about user functions. 
; It returns the weighted deviation between the model and the data,
; and is used when derivatives are computed by finite differences.
; INPUT:
;   PARS:      array of all the parameters of the model. Positional
;              parameter. It's dimensions is [starting guesses+ ],
;              example: 2 moments of the gauss(v, sigma)+
;              4*1(nBursts) = 14
;   KPEN, POLPEN, ADEGREE, MDEGREE, KMOMENT, VOFF, SPECTRUM, GOODPIXELS
;              check the documentation in uly_fit
;   LUM_WEIGHT  used only when called from ULY_FIT
; OUTPUT:
;   OUTPIXELS, BESTFIT, MPOLY, ADDCONT, CMP 
;              These output arguments are used only when called from ULY_FIT,
;      
; RETURN:
;   The user function returns the wighted deviation between
;   the model and the data for the pixels in GOODPIXELS
;
FUNCTION uly_fitfunc, pars,                                      $
                      KPEN=kpen, POLPEN=polpen,                  $
                      ADEGREE=adegree, MDEGREE=mdegree,          $
                      KMOMENT=kmoment,                           $
                      VOFF=voff,                                 $
                      SPECTRUM=signalLog, GOODPIXELS=goodpixels, $
                      QUIET=quiet,                               $
                      LUM_WEIGHT=lum_weight,                     $
                      OUTPIXELS=outpixels, BESTFIT=bestfit,      $
                      MPOLY=mpoly, ADDCONT=addcont, CMP=cmp

compile_opt idl2, hidden

if n_elements(pars) gt 0 then if min(finite(pars)) eq 0 then message, '(pars) array contains NaN'

; OUTPIXELS keyword is used for the cleaning procedure in ULY_FIT
if n_elements(outpixels) eq 0 then outpixels=goodpixels

; MPOLY : multiplicative continuum polynomial
; structure {lmdegree, mpolcoefs, poly}
;   lmdegree (integer), maximum degree of legendre polynomials
;   mpolcoefs (double 1D array), LMDEGREE+1 terms.
;   poly (1D array)
; input: mpolcoefs are the coefficients determined at previous iterration
; output: mpolcoefs contains updated values of the coefficients
; The structure is initialized in the function uly_fit

; GOODPIXELS : combine the good pixels of model and observation

;Mina Koleva 2005/09/11
;Separate the evaluate procedure from the fitfunc
;call evaluate procedure

; The heavy data grids are passed through common, because the
; functargs mechanism used by MPFIT is inefficient
common uly_functargs_common, mpolyc, cmpc, cache, bvls_state, modecvg

uly_fit_pparse, pars, par_losvd, cmpc, KMOMENT=kmoment, QUIET=quiet

bestfit = uly_fit_lin(ADEGREE=adegree,                             $
                      POLPEN=polpen,                               $
                      VOFF=voff,                                   $
                      PAR_LOSVD=par_losvd, CMP=cmpc,               $
                      GOODPIXELS=goodPixels, SPECTRUM=signalLog,   $
                      MPOLY=mpolyc,                                $
                      ADDCONT=addcont,                             $
                      CACHE=cache,                                 $
                      STATE=bvls_state,                            $
                      QUIET=quiet,                                 $
                      MODECVG=modecvg,                             $
                      LUM_WEIGHT=lum_weight,                       $
                      STATUS=status)

if modecvg eq 1 then modecvg = 3  ; prevent iteration when computing the derivative, 

mpoly = mpolyc
cmp = cmpc
if status ne 0 then begin
   common MPFIT_ERROR, error_code
   error_code = -2              ; MPFIT fatal error
   return, 0
endif

; Compute 'err' the weighted deviates between the model and the observation.
err = double(((*signalLog.data)[outpixels]-bestfit[outpixels]) / (*signalLog.err)[outpixels])

; Penalize the solution towards (h3,h4,...)=0 if the inclusion of
; these additional terms does not significantly decrease the error.
if kmoment gt 2 and kpen ne 0 then $
  err = err + kpen*robust_sigma(err, /ZERO)*sqrt(total(pars[2:kmoment-1]^2,/DOUBLE))

return, err
END

;----------------------------------------------------------------------------
; This function is executed at each LM iteration.
PRO uly_iterproc, myfunct, p, iter, fnorm, FUNCTARGS=fcnargs, $
                  PARINFO=parinfo, DOF=dof,                   $
                  KPEN=kpen, POLPEN=polpen,                   $
                  ADEGREE=adegree, MDEGREE=mdegree,           $
                  KMOMENT=kmoment,                            $
                  VOFF=voff,                                  $
                  SPECTRUM=signalLog, GOODPIXELS=goodpixels,  $
                  QUIET=quiet

compile_opt idl2, hidden
  
common uly_functargs_common, mpolyc, cmpc, cache, bvls_state, modecvg

if modecvg eq 3 then modecvg = 1 ; make uly_fit_lin converge when calling at a new point

END


;----------------------------------------------------------------------------
FUNCTION uly_fit, signalLog,                                            $
                  CMP=cmp,                                              $
                  KMOMENT=kmoment, KGUESS=kguess, KFIX=kfix, KLIM=klim, $
                  KPEN=kpen, $
                  ADEGREE=adegree,                                      $
                  MDEGREE=mdegree, POLPEN=polpen,                       $
                  CLEAN=clean,                                          $
                  MODECVG=modecvg,                                      $
                  QUIET=quiet,                                          $
                  ALL_SOLUTIONS=all_solutions                          

compile_opt idl2
on_error, 2

c0 = 299792.458d    ;speed of the light in km/s

galaxy = *signalLog.data
npix = (size(galaxy))[1]

have_noise = 1
if n_elements(*SignalLog.err) eq 0 then begin
   have_noise = 0
   *SignalLog.err = 1 + 0 * galaxy
endif 

noise = *SignalLog.err 

; combine the goodpixels from observation and model
msk = uly_spect_get(SignalLog, /MASK)
for i=0,n_elements(cmp)-1 do begin
    if ptr_valid(cmp[i].mask) then $
      if n_elements(*cmp[i].mask) gt 0 then msk *= *cmp[i].mask
endfor

goodpix = uly_spect_get(SignalLog, /goodpix)
nan = where(finite(galaxy[goodpix]) eq 0, c)
if c gt 0 then begin
    if not keyword_set(quiet) then $
      message, strtrim(string(c),2)+' NaNs are masked in the spectrum', /INFO
    msk[goodpix[nan]] = 0
endif
nan = where(finite(noise[goodpix]) eq 0, c)
if c gt 0 then begin
    if not keyword_set(quiet) then $
      message, strtrim(string(c),2)+' NaNs are masked in the error spectrum', /INFO
    msk[goodpix[nan]] = 0
endif
goodPixels0 = where(msk eq 1, cnt)
if cnt eq 0 then begin
    if not keyword_set(quiet) then message, /INFO, 'No good pixel left'
    return, 0
endif

; Check that all the cmp have the same WCS
for k = 1, n_elements(cmp)-1 do begin
    if abs(cmp[k].start-cmp[0].start) gt 0.01*cmp[0].step then $
      message, 'Component '+strtrim(string(k),2)+' has a different start than cmp0'+ string(cmp[k].start) + string(cmp[0].start)
    if abs(cmp[k].step-cmp[0].step)*cmp[0].npix gt 0.01*cmp[0].step then $
      message, 'Component '+strtrim(string(k),2)+' has a different step than cmp0'+ string(cmp[k].step) + string(cmp[0].step)
endfor

; make a reduction for velocity, because the starting wl of
; the signal and the templates are not the same
voff = (cmp[0].start-signalLog.start)*c0 ; km/s

velScale = c0 * signalLog.step

; Do some input error checking
;
s2 = size(galaxy, /DIM)
s3 = size(noise, /DIM)

if (size(galaxy, /N_DIM) ne 1 or size(noise, /N_DIM) ne 1) then $
  message, 'Input data and noise should be 1D arrays'
if not array_equal(s2,s3) then begin
    message, 'Observation and noise spectra must have the same size ('+string(s2)+' vs. '+string(s3)+')'
endif
if (min(cmp.npix) lt s2[0]) then message, 'MODELS length cannot be smaller than observation' + $
  ' (models:'+strtrim(cmp.npix,2)+', data:'+strtrim(s2[0],2)+')'
if n_elements(adegree) le 0 then adegree = -1 else adegree = adegree
if n_elements(mdegree) le 0 then mdegree = 10 else mdegree = mdegree
if max(goodPixels0) gt s2[0]-1 then message, 'GOODPIXELS are outside the data range'
if n_elements(kpen) le 0 then kpen = 0.7d*sqrt(500d/n_elements(goodPixels0))
if n_elements(kmoment) eq 0 then kmoment = 2 else $
    if total(kmoment eq [0,2,4,6]) eq 0 then message, 'KMOMENT should be 0, 2, 4 or 6'
if kmoment ge 2 and n_elements(kguess) lt 2 then message, 'KGUESS must have two elements [V,sigma]'

nst=n_elements(kguess)
if n_elements(kfix) eq 0 then kfix=intarr(nst)*0 $ ;setting fixed flags to 0
else if n_elements(kfix) lt nst then $
   kfix = [kfix, intarr(nst-n_elements(kfix))] $
else if n_elements(kfix) gt nst then kfix = kfix[0:nst-1]

if kmoment gt 0 then begin
;  Set default limits for kinematics
;  The convolution kernel for the LOSVD is bad for sub-pixel sigma, for this 
;  reason we set the default low limit of sigma to 0.3 pix (FWHM=0.7 px). 
   klimits = dblarr(2, kmoment)
   klimits[0:1,0] = (kguess[0] + [-2d3,2d3]) ; +/-2000 km/s 
   if kmoment ge 2 then klimits[0:1,1] = [0.3*velScale,1d3]
   for k=2,kmoment-1 do klimits[0:1,k] =  [-0.3d,0.3d]
;  Override with the actual limits, if given
   if n_elements(klim) ne 0 then klimits[0] = klim
;  Diagnostic out of bounds
   for k=0,kmoment-1 do begin
      if kfix[k] eq 0 and kguess[k] lt klimits[0,k] then begin
         message, /CONT, 'Guess on kinematic moment'+string(k)+' :'+ $
                  string(kguess[k])+' is lower that the limit:'+ $
                  string(klimits[0,k])
         return, 0
      endif else if kfix[k] eq 0 and kguess[k] gt klimits[1,k] then begin
         message, /CONT, 'Guess on kinematic moment'+string(k)+' :'+ $
                  string(kguess[k])+' is lower that the limit:'+ $
                  string(klimits[1,k])
         return, 0
      endif
   endfor
;  Convert klimits into pixels
   klimits[0:1,0] = alog(1 + klimits[0:1,0]/c0) / signalLog.step
   if kmoment ge 2 then klimits[0:1,1] = alog(1 + klimits[0:1,1]/c0) / signalLog.step
endif

lmdegree = mdegree

n_nodes = 1
for n1=0,n_elements(cmp)-1 do for n2=0,n_elements(*(cmp[n1].para))-1 do $
  n_nodes *= n_elements(*(*cmp[n1].para)[n2].guess)

if not keyword_set(quiet) and n_nodes gt 1 then $
  print, 'Perform global optimization, number of nodes:', n_nodes

n_para = 0
for n1=0,n_elements(cmp)-1 do n_para += n_elements(*(cmp[n1].para))
if n_para gt 0 then indx = intarr(n_para, 3)
i = 0
for n1=0,n_elements(cmp)-1 do begin
    for n2 = 0, n_elements(*(cmp[n1].para))-1 do begin
        indx[i,0] = n1
        indx[i,1] = n2
        i++
    endfor
endfor

; In case n_nodes>1 and keyword_set(clean) we should do a first loop 
; on nodes with cleaning, and a second one with the best mask that was
; determined. *TO BE IMPLEMENTED*
; iter_global_clean = n_nodes gt 1 and keyword_set(clean) ? 2 : 1

for node=0,n_nodes-1 do begin   ; loop on the nodes of the guess grid

    goodPixels = goodPixels0

;   Test if the input parameters (guesses) are withing the permitted limits
    for i=0,n_elements(cmp)-1 do begin
       for k = 0, n_elements(*cmp[i].para)-1 do begin
          if (*(*cmp[i].para)[k].guess)[0] lt (*cmp[i].para)[k].limits[0] then begin
             message, /CONT, 'guess on '+(*cmp[i].para)[k].name + ' :' + $
                      string((*(*cmp[i].para)[k].guess)[0]) + $
                      ' is below the limit :' + $
                      string((*cmp[i].para)[k].limits[0])
             return, 0
          endif else if $
             (*(*cmp[i].para)[k].guess)[0] gt (*cmp[i].para)[k].limits[1] then begin
             message, /CONT, 'guess on '+(*cmp[i].para)[k].name + ' :' + $
                      string((*(*cmp[i].para)[k].guess)[0]) + $
                      ' is above the limit :' + $
                      string((*cmp[i].para)[k].limits[1])
             return, 0
          endif
       endfor
    endfor
    
;   make parinfok : the part of parinfo containing the kinematics 
    if kmoment eq 0 then undefine, parinfok
    if kmoment gt 0 then begin
        parinfok = REPLICATE({value:0D, step:1D-2, limits:[0D,0D], limited:[1B,1B], fixed:0B}, kmoment)
        parinfok[0].value = alog(1 + kguess[0]/c0) / signalLog.step ; Convert vel to pix
        parinfok[0].limits = klimits[*,0]
        if(kfix[0] eq 1) then begin ;; radial velocity
            parinfok[0].fixed = 1
            parinfok[0].step = 0
            parinfok[0].limits = [0.,0.] + parinfok[0].value
        endif
    endif
    if kmoment gt 1 then begin
        parinfok[1].value = alog(1 + kguess[1]/c0) / signalLog.step ; Convert vel to pix
        parinfok[1].limits = klimits[*,1]
        if(kfix[1] eq 1) then begin ;; velocity dispersion
            parinfok[1].fixed=1.0
            parinfok[1].step=0
            parinfok[1].limits=[0.,0.] + parinfok[1].value
        endif
        if min(cmp.npix) le (abs(voff)+abs(kguess[0])+5d*kguess[1])/Velscale then begin
            message, /CONT, $
              'Wavelength range is too small, or velocity shift too big'
            return, 0
        endif
    endif 
    if kmoment gt 2 then begin
       parinfok[2:kmoment-1].value=kguess[2:*]
       parinfok[2:kmoment-1].limits = [-0.3d,0.3d] ; -0.3<[h3,h4,...]<0.3
       parinfok[2:kmoment-1].step = 1d-3
       fixh = where(kfix eq 1, cnt)
       if cnt gt 0 then begin
          parinfok[fixh].fixed=1.0
          parinfok[fixh].step=0
          for ii = 0,cnt-1 do $
             parinfok[fixh[ii]].limits=[0.,0.] + kguess[fixh[ii]]
       endif
    endif
;   make parinfo : the cmp part
    parinfo = uly_makeparinfo(cmp)
    if size(parinfo, /TYPE) ne 8 then undefine, parinfo
;   combine parinfo with parinfok
    if n_elements(parinfok) gt 0 then $
       if size(parinfo, /TYPE) eq 8 then parinfo = [parinfok,parinfo] $
       else parinfo = [parinfok]

; Minimization for one set of guesses
;      If required by the /CLEAN kw, once the minimum is found, clean the 
;      outliers and repeat the minimization
    nclean = 0
    if keyword_set(clean) then nclean = 10 ; At most 11 cleaning iterations
    if nclean eq 0 then ftol=1d-5 else ftol=1d-2

    for j=0,nclean do begin
        time0=systime(1)
    
;       Parameters to be passed to the fitting function through functArgs
        if n_elements(quiet) eq 0 then qui = 0 else qui=1 
        if n_elements(polpen) ne 0 then $
          functArgs = {KMOMENT:fix(kmoment), KPEN:kpen,                $
                       ADEGREE:adegree,                                $
                       SPECTRUM:signalLog, GOODPIXELS:goodPixels,      $
                       VOFF:voff/velScale,                             $
                       POLPEN:polpen, QUIET:qui}                       $
        else $
          functArgs = {KMOMENT:fix(kmoment), KPEN:kpen,                $
                       ADEGREE:adegree,                                $
                       SPECTRUM:signalLog, GOODPIXELS:goodPixels,      $
                       VOFF:voff/velScale,                             $
                       QUIET:qui}

;       initialize the common variables
        uly_fitfunc_init, SPECTRUM=signalLog, CMP=cmp, MDEGREE=mdegree, MODECVG=modecvg

;       Keep the residuals before the minimization, to check later
;       the magnitude of the improvment.
        undefine, value
        if size(parinfo, /TYPE) eq 8 then value = parinfo.value
        resc0 = uly_fitfunc(value, KPEN=kpen, ADEGREE=adegree,      $
                            SPECTRUM=signalLog,                     $
                            GOODPIXELS=goodPixels,                  $
                            VOFF=voff/velScale,                     $
                            KMOMENT=kmoment,                        $
                            QUIET=quiet)    
        if n_elements(resc0) le 1 then return, 0

        nlin_free = 0  ; number of free non-linear parameters
        for k=0,n_elements(parinfo)-1 do begin
           if (parinfo[k]).fixed eq 0 then nlin_free += 1
        endfor

        if nlin_free gt 0 then begin
;        /NOCATCH : disable the error catching wich makes uneasy the debuging in
;        the user's function. The alternatine would be to have a traceback
;        message like: http://www.dfanning.com/widget_tips/widgettrace.html
;        implemented in the error handling of MPFIT (modif of MPFIT).
           res = mpfit('uly_fitfunc', PARINFO=parinfo, FUNCTARGS=functArgs, $ 
                       ITERPROC='uly_iterproc', ITERARGS=functArgs,         $ 
                       FTOL=ftol, XTOL=1d-10,                               $
                       PERROR=error, BESTNORM=bestnorm,                     $
                       STATUS=mpfstat, ERRMSG=errmsg, NFEV=ncalls,          $
                       /QUIET, /NOCATCH)

           uly_fit_pparse, res, par_losvd, cmp,                         $
                           ERROR=error, KMOMENT=kmoment, QUIET=quiet

           if not keyword_set(quiet) then begin
              if mpfstat eq 5  then $
                 print,'MPFIT STATUS= 5 (max number of iterations reached)' $
              else if mpfstat eq -1 then $ ; code set by ULY_FIT_LIN, indicates
                                ; that the best fit is the null spectrum
                 print,'ULY_FIT: bestfit is the NULL spectrum' $
              else if mpfstat lt 1 or mpfstat gt 3 then $
                 print,'MPFIT STATUS=', strtrim(mpfstat,2), ' ', errmsg
              if mpfstat eq 0 then begin
;                at least one of the input parameter is not within the assigned
;                limits, this should not happen, because it is tested before calling MPFIT
                 for k = 0, n_elements(parinfo)-1 do begin
                    print, 'param:', k, ' value:', parinfo[k].value, ' limits:', parinfo[k].limits
                 endfor
              endif
           endif
           
           if errmsg ne '' and mpfstat ne -1 then begin
              if not keyword_set(quiet) then $
                 message, 'MPFIT returned to ULY_FIT with an error: ' + errmsg, /CONT
              return, 0
           endif
           
           if min(finite(res)) eq 0 then begin
              if not keyword_set(quiet) then $
                 message, 'MPFIT did not return a valid result', /CONT
              return, 0
           endif

           if not keyword_set(quiet) and n_elements(ncalls) ne 0 then $
              print, 'number of model evaluations:', ncalls

        endif else begin
           bestnorm = total((resc0)^2)  
           if size(parinfo,/TYPE) eq 8 then begin
               res = parinfo.value ; the NL parameters have their input value
               error = dblarr(n_elements(parinfo))
           endif else error = 0
        endelse
                
        if j eq nclean then break
        
;;;;;  Next block is cleaning: recompute the goodpixel list
;;;;;  ;;;;;;;;;;;;;;;;;;;
        goodOld = goodPixels ; save the goopix list to see if clipping change
        rbst0=stddev(resc0)
        resc = fltarr(npix)

        resc[goodPixels0] = uly_fitfunc(res, KPEN=kpen, ADEGREE=adegree, $
                                        SPECTRUM=signalLog,              $
                                        GOODPIXELS=goodPixels,           $
                                        VOFF=voff/velScale,              $
                                        KMOMENT=kmoment,                 $
                                        OUTPIXELS=goodPixels0,           $
                                        BESTFIT=bestfit,                 $
                                        QUIET=quiet)
    
;       compute the gradients in the model.
;       this is used to prevent clipping of pixels that may be explained
;       by a shift of about 1/5 of pixel
        facsh = 1/5d
        if j eq 0 then facsh = 0.5d
        modelgrd = max([[abs(bestfit-shift(bestfit,-1))], $
                        [abs(bestfit-shift(bestfit,1))]], DIMENSION=2) $
                   *facsh / noise ; gradients in model

;       select errors larger than 3*sigma or 4, 5, 7 depending on
;       clipped fraction
        rbst_sig=stddev(resc[goodPixels])
        tmp = where(abs(resc[goodPixels])-modelgrd[goodPixels] gt 3*rbst_sig, m, COMPLEM=w) 
        clip_level = 3
        if (m gt 0.03*n_elements(goodPixels)) then begin ; dont remove too many
           tmp = where(abs(resc[goodPixels])-modelgrd[goodPixels] gt 4*rbst_sig, m, COMPLEM=w)
           clip_level = 4
           if (m gt 0.03*n_elements(goodPixels)) then begin
              tmp = where(abs(resc[goodPixels])-modelgrd[goodPixels] gt 5*rbst_sig, m, COMPLEM=w)
              clip_level = 5
              if (m gt 0.03*n_elements(goodPixels)) then begin
                 tmp = where(abs(resc[goodPixels])-modelgrd[goodPixels] gt 7*rbst_sig, m, COMPLEM=w)
                 clip_level = 7
              endif
           endif 
        endif 

        mask = fltarr(npix)
        mask[goodPixels0] = 1

        if m ne 0 then begin
           tmp = where(abs(resc[goodPixels0])-modelgrd[goodPixels0] gt clip_level*rbst_sig, m, COMPLEM=w)
           tmp = goodPixels0[tmp]
           mask[tmp] = 0b
        endif

        m2 = 0
        if (m ne 0) then begin ; clean the 1st neigbours at a 2*sig threshold
            near = [tmp-1, tmp+1]
            near = near[where(near ge 0 and near lt npix)]
            nnnn = where(abs(resc[near])-modelgrd[near] gt 2*rbst_sig and mask[near] eq 1b, m2)
            if m2 gt 0 then begin
               mask[near[nnnn]] = 0b
               tmp = [tmp, near[nnnn]]
            endif
        endif

        r_sig = 0
        if (m ne 0) then begin ;
            nclip = 0
            for k=1, 20 do begin
               r_sig=stddev(resc[where(mask eq 1)])
               ss = shift(resc, 1)
               nnn1 = where(mask eq 1 and shift(mask,1) eq 0 and ss*resc gt 0 and abs(resc)-modelgrd gt r_sig and abs(resc) le abs(ss), mmm1)
               ss = shift(resc, -1)
               nnn2 = where(mask eq 1 and shift(mask,-1) eq 0 and ss*resc gt 0 and abs(resc)-modelgrd gt r_sig and abs(resc) le abs(ss), mmm2)
               if mmm1 gt 0 then mask[nnn1] = 0b
               if mmm2 gt 0 then mask[nnn2] = 0b
               if mmm1 le 0 and mmm2 le 0 then break
               nclip += mmm1 + mmm2
            endfor
        endif
        
        if (m ne 0) then begin
            if not keyword_set(quiet) then $
               print, 'Number of clipped outliers:', strtrim(m,2), '+', $
                      strtrim(m2,2),'+',strtrim(nclip,2), $
              ' out of',n_elements(goodPixels0), rbst_sig
            goodPixels = where(mask eq 1)
        endif

        if array_equal(goodOld,goodPixels) then break
        ftol=1d-5  ; Normal precision (high) for the last iter
        if j+2 lt nclean then ftol=1d-3 $
        else if j+1 lt nclean then ftol=1d-4

; rbst0 : standard deviation of the residuals before the minimisation
; rbst_sig : standard deviation of the residuals after minimisation 
; rbst1 : standard deviation of the residuals after minimisation and clipping
; rbst_sig^2/(rbst0*rbst1) : ratio of the reduction due to the
;          optimization to the reduction due to the clipping
; If there is a major reduction of the redicuals due to the clipping,
; this ratio becomes large, and as the effect of optimizing
; the parameters have been minor, we restart from the previous
; position in the parameters space. Otherwise we continue from
; the current 'res' position.
        rbst1=stddev(resc[goodPixels])
        if rbst_sig^2/(rbst0*rbst1) lt 1.5 and size(parinfo,/TYPE) eq 8 $
          then  parinfo.value = res $
        else if j+1 lt nclean then ftol=1d-2
        
;;; End of the cleaning block ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    endfor

; get the output BESTFIT, WEIGHTS, ADDCONT, MPOLY
;

;   MPFIT assumes the bins are independent ... we have to correct the errors
    error *= sqrt(SignalLog.dof_factor) 
    uly_fit_pparse, res, par_losvd, cmp,                         $
      ERROR=error, KMOMENT=kmoment, QUIET=quiet


    deviates = uly_fitfunc(res, KPEN=kpen, ADEGREE=adegree, $
                           SPECTRUM=signalLog,              $
                           GOODPIXELS=goodPixels,           $
                           VOFF=voff/velScale,              $
                           KMOMENT=kmoment,                 $
                           OUTPIXELS=goodPixels,            $
                           BESTFIT=bestfit, CMP=cmp,        $
                           ADDCONT=addcont, MPOLY=mpoly,    $
                           /LUM_WEIGHT,                     $
                           QUIET=quiet)
    bestnorm = total(deviates^2)

; Computation of the reduced chi square: chi2
;   bestnorm = total(dev^2), determined by mpfit, where
;   dev = (observation-model)/err, computed by the user function uly_fitfunc.
;   Therefore, chi2 = bestnorm/nu, where nu is the number of degrees of freedom.
;   nu = (number of pixels) - (npara, number of free parameters) 
    para = kfix
    for k=0,n_elements(cmp)-1 do $
       if n_elements((*cmp[k].para)) gt 0 then para = [para, (*cmp[k].para).fixed]
    p = where(para eq 0, npara)
    npara += kmoment + n_elements(cmp) + mdegree + adegree 
    chi2 = bestnorm  / (n_elements(goodPixels)-npara)

; If no noise was provided, chi2 is meaningless. 
; We compute a S/N corresponding to the observed residuals (as if chi2=1), 
; and rescale the errors (the computing errors are hence upper estimates,
; under the assumption that the errors result only for the noise).
    snr_est = 0
    if have_noise eq 0 then begin
       snr_est = mean(galaxy[goodPixels])/sqrt(chi2)
       error *= sqrt(chi2)
       i0 = kmoment
       for i=0,n_elements(cmp)-1 do begin
          i1 = i0 + n_elements(*cmp[i].para) - 1
          if i1 ge i0 then (*cmp[i].para).error *= sqrt(chi2)
          i0 = i1 + 1
       endfor
    endif
    
; Load the kinematics in the cmp structure
    if kmoment gt 0 then begin
        losvd = c0 * (exp(signalLog.step * res[0:1]) - 1)
        e_losvd = error[0:1]*velScale
        if kmoment gt 2 then begin ; do not multiply h3,h4,.. by VelScale
           losvd = [losvd, res[2:kmoment-1]]
           e_losvd = [e_losvd, error[2:kmoment-1]]
        endif
        s_losvd = 0*kfix  ; to hold the status (0/1/2/3 for fixed/pegged)
        n = where(kfix eq 1, cnt)
        if cnt gt 0 then s_losvd[n] = 1
        n = where((res[0:kmoment-1] eq klimits[0,*]) eq 1, cnt)
        if cnt gt 0 then s_losvd[n] = 2
        n = where((res[0:kmoment-1] eq klimits[1,*]) eq 1, cnt)
        if cnt gt 0 then s_losvd[n] = 3
    endif 

; Create the output cmp array to attach to the output solution
;   We do not simply attach cmp into the output solut structure for
;   two reasons:
;     1) cmp contains the pointers  .init_data .eval_data and .mask
;     that are used by the input cmp, and if we copy the pointers
;     we will not know how to manage the cmp (we cannot free one without
;     damaging the others).
;     Therefore in the output copy of cmp we put these pointers to NULL.
;     2) The pointer para contains the data that we want to save.
;     We have to copy these data, not simply the pointer.
; The output solut structure can be safely freed with heap_free.
    cmpo = cmp

    for k=0,n_elements(cmpo)-1 do begin
       cmpo[k].init_data = ptr_new()
       cmpo[k].eval_data = ptr_new()
       cmpo[k].mask = ptr_new()
       cmpo[k].para = ptr_new(*cmp[k].para)

       if n_elements(*cmpo[k].para) gt 0 then begin
          nt = where(tag_names(*cmpo[k].para) eq 'STATUS')
          if nt[0] eq -1 then begin ; add the tag STATUS in the structure
             cmptmp = create_struct((*cmpo[k].para)[0], 'status', 0S)
             cmptmp = replicate(cmptmp, n_elements(*cmpo[k].para))
             struct_assign, *cmpo[k].para, cmptmp
             *cmpo[k].para = cmptmp
          endif

          for l=0,n_elements(*cmpo[k].para)-1 do begin
             (*cmpo[k].para)[l].guess = ptr_new((*(*cmpo[k].para)[l].guess)[0])
             
             (*cmpo[k].para)[l].status = 0
             if (*cmpo[k].para)[l].fixed eq 1 then (*cmpo[k].para)[l].status = 1 $
             else if (*cmpo[k].para)[l].value le (*cmpo[k].para)[l].limits[0] then $
                (*cmpo[k].para)[l].status = 2 $
             else if (*cmpo[k].para)[l].value ge (*cmpo[k].para)[l].limits[1] then $
                (*cmpo[k].para)[l].status = 3 
          endfor
       endif
    endfor

; Create the SOLUTION output structure:
    mask = bytarr(npix)
    mask[goodPixels] = 1
    h = n_elements(*signalLog.hdr) eq 0 ? ['END'] : *signalLog.hdr

    if kmoment gt 0 then $
      sol = {                                  $
              title:signalLog.title,           $
              hdr:ptr_new(h),                  $
              start:signalLog.start,           $
              step:signalLog.step,             $
              sampling:signalLog.sampling,     $
              refpix:1.,                       $
              data:ptr_new(*signalLog.data),   $
              err:ptr_new(*signalLog.err),     $
              bestfit:bestfit,                 $
              mask:mask,                       $
              dof_factor:signalLog.dof_factor, $
              chisq:chi2,                      $
              snr:snr_est,                     $
              addcont:addcont,                 $
              mulcont:mpoly.poly,              $
              model_id:cmp[0].name,            $
              mulcoef:mpoly.mpolcoefs,         $
              losvd:losvd,                     $
              e_losvd:e_losvd,                 $
              s_losvd:s_losvd,                 $
              cmp:cmpo                         $
            } $
    else $
      sol = {                                  $
              title:signalLog.title,           $
              hdr:ptr_new(h),                  $
              start:signalLog.start,           $
              step:signalLog.step,             $
              sampling:signalLog.sampling,     $
              refpix:1.,                       $
              data:ptr_new(*signalLog.data),   $
              err:ptr_new(*signalLog.err),     $
              bestfit:bestfit,                 $
              mask:mask,                       $
              dof_factor:signalLog.dof_factor, $
              chisq:chi2,                      $
              snr:snr_est,                     $
              addcont:addcont,                 $
              mulcont:mpoly.poly,              $
              model_id:cmp[0].name,            $
              mulcoef:mpoly.mpolcoefs,         $
              cmp:cmpo                         $
            } 

    if not keyword_set(quiet) and n_nodes gt 1 then $
      print, 'node:', node, ' chi2: ', sol.chisq

    if not keyword_set(all_solutions) then begin
       if n_elements(solution) gt 0 then if sol.chisq lt solution.chisq then begin
          heap_free, solution
          solution = sol
       endif else heap_free, sol
       if n_elements(solution) eq 0 then solution = sol    
    endif else begin
       if n_elements(solution) gt 0 then solution = [solution, sol] $
       else solution = sol
    endelse

; use another guess if required
    for n=0, n_para-1 do begin
        n1 = indx[n,0]
        n2 = indx[n,1]
        *(*cmp[n1].para)[n2].guess = shift(*(*cmp[n1].para)[n2].guess, -1)
        if indx[n,2] ne n_elements(*(*cmp[n1].para)[n2].guess)-1 then begin
            indx[n,2] ++
            for k=0,n-1 do indx[k,2] = 0
            break
        endif
    endfor
endfor                              ; end of the loop on the grid of guesses


if have_noise eq 0 then begin
   solution.chisq  = 0
;   undefine, *sol.err   ; solut_splot will fail if *sol.err is suppressed (shall debug)
endif

return, solution

end
;-- end -----------------------------------------------------------------------
