;+
; NAME:                 ULY_LSF
;
; PURPOSE:              Determine the line-spread function (LSF)
;
; USAGE:                uly_lsf, <spectrum>, <cmp>, <width>, <step>,   $
;                            FILE_OUT=file_out, /PLOT, /QUIET, /CLEAN, $
;                            KMOMENT=<kmoment>, DG=<dg>, DL=<dl>,      $
;                            KFIX=<kfix>, KPEN=<kpen>,                 $
;                            MDEGREE=<mdegree>, ADEGREE=<adegree>,     $
;                            VELSCALE=<velscale>
;
; ARGUMENTS:
;     <spectrum>       Filename or spectrum structure containing the spectrum
;                      to be analysed for the LSF.
;     <cmp>            Component structure containing the reference high
;                      resolution spectrum. The LSF of <spectrum> will be
;                      determined relatively to <cmp>
;                      In the simplest case, <cmp> can be a reference spectrum
;                      loaded with ULY_STAR, or it can be an interpolated
;                      stellar model (see ULY_TGM) or any other model, or even
;                      an array of several cmp.
;                      During the analysis, eventual non-linear parameters
;                      will be adjusted.
;     <width>, <step>  scalars, [Angstrom] 
;                      The LSF is measured in a succession of chunks of width 
;                      <width> separated by <step>.
;         
; KEYWORDS:
;     FILE_OUT         Name of the output file.
;                      Defaulted to 'output.lsf'
;     KMOMENT, DG, DL, KFIX, KPEN
;                      LOSVD parameters passed to ULySS.
;                      They are used to determine which terms of the LSF to fit.
;     MDEGREE, ADEGREE, VELSCALE
;                      Parameters passed to ULySS.
;                      Remember that MDEGREE (and ADEGREE) apply to the (small)
;                      chunks of spectra of width <width>, so, this parameters
;                      should be smaller than when the whole spectrum is
;                      analysed at once.
;     PLOT             Set this keyword if you want to display the
;                      losvd fit using the program ULY_SOLUT_SPLOT.
;     CLEAN
;                      Set this keyword to perform a kappa-sigma clipping
;                      of the outliers.
;     QUIET            Verbosity control
;
; EXAMPLE:
;   In the following example, we determine the relative LSF between a
;   spectrum to analyse and a model extracted from a grid.
;   Then we smooth this LSF and use it to analyse the spectrum.
;
;   filetest = uly_root+'/data/VazMiles_z-0.40t07.94.fits.fits'
;   model = uly_ssp_extr(filetest, uly_root+'/models/PHR_Elodie31.fits' ,[6000.,-0.4,0.], SIG=30)
;   cmp = uly_star(model)
;   uly_lsf, filetest, cmp, 200, 100, FIL='vaz.lsf', /QUIET
;   heap_free, cmp
;   uly_lsf_plot, 'vaz.lsf' 
;   uly_lsf_smooth, 'vaz.lsf', 'vazs.lsf'
;   uly_lsf_plot, 'vazs.lsf' 
;   cmp = uly_ssp(LSF='vazs.lsf')
;   ulyss, filetest, cmp, /PLOT
;
; LSF ANALYSIS AND INJECTION:
; 
;   A spectrum is characterized by its line spread function, LSF, ie. the 
;   instrumental broadening of the spectral lines. This is the analog of the 
;   PSF for images. The LSF can be modeled by a Gaussian or more generally 
;   by a Gauss-Hermite expansion of a Gaussian.
;   
;   For an observation, the LSF depends on the spectrograph and on the slit, 
;   it generally changes with the wavelength and, in the case of 2 or 3D 
;   spectra, with the position in the field. On the other side, the model
;   used to analyse this observation has also a finite resolution, i.e.
;   has a specific LSF hopefully finer than the one of the observation.
;   Therefore to analyse an observation we have to take into account the
;   relative LSF of the observation with respect to the model. This is of 
;   a primordial importance when studying the internal kinematics (physical 
;   broadening of the lines) of an astronomical source.  
;   
;   The process to properly take the LSF into consideration consists in 
;   two phases: (i) The determination of the relative LSF and (ii) its 
;   injection into the model, so that the result is comparable to an
;   observed spectrum.  
;   
;   Choice of the LSF standard
;   
;   In standard cases, the LSF can be considered as stable for a set of 
;   observations, for example for an observing run, or for a given setup of 
;   an instrument. In this case it is possible to determine it once using 
;   the calibration and standard stars observations. There are
;   different ways to determine the LSF.
;   
;   1) Using the arc spectrum used for the wavelength calibration. 
;   The profile of the lines can be determined at different wavelength 
;   and position in the field (sigma, h3, h4). 
;   This method is currently not implemented in ULySS and requires 
;   some precaution because of blends of the lines and limited sampling.  
;   
;   2) Using a twilight spectrum, and compare it with the solar spectrum 
;   having the same LSF as the model. The input aperture is uniformly
;   illuminated. This is the preferred method.
;   
;   3) Using a spectrum from a star for which a spectrum with the same LSF 
;   as the model is available. The critical aspect of this method is 
;   to be sure that the stellar light is uniformly distributed over
;   then entrance aperture or that the resultion is not determined by
;   the size and shape of the entrance aperture.
;   
;   4) Using a spectrum of a star of known atmospheric parameter,
;   and analyse it with the output of a library interpolator. Or use a
;   spectrum of a star with unknown atmospheric parameters and determine
;   first these parameters with ULySS (see ULY_TGM).
;   
;   Finally, if no stellar spectrum is available, the analysis can be done 
;   with an observation of a galaxy and the internal kinematics of all 
;   objects in the series will be determined relatively to this reference. 
;   (choose the galaxy with the smallest velocity dispersion or decrease 
;   the width of the measured LSF by a fixed amount (use the keyword SIGMA).
;
;   Caveat: To maximize the flux efficiency, some setups use a wide input 
;   aperture (slit) that may be eventually wider than the seeing and that 
;   the intrinsic resolution of the spectrograph. The result is that in 
;   case of good seeing the object may be resolved in the slit, and the 
;   measured broadening does not reflect the physical broadening. This is 
;   a poor practice, and such observations should not be used to study 
;   the internal kinematics.  
;   
;   Analysis of the LSF 
;
;   The analysis of the LSF is made with the program: ULY_LSF
;   which produce an ASCII file describing the LSF (see the
;   Description section in this article). 
;   This file may be viewed with ULY_LSF_PLOT, or smoothed with ULY_LSF_SMOOTH.
;
;   Injection of the LSF 
; 
;   The LSF injection is made by: ULY_SPECT_LSFCONVOL, normally called when
;   a component of the model is initialized or evaluated.
;   The principle of the injection of the LSF in the model is briefly 
;   introduced in Koleva et al. 2008.
;
; HISTORY:
;   Mina Koleva,       2008/03/11, created
;   Yue WU,            2009/03/05, updated
;   Philippe Prugniel, 2015/08/10, updated
;
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
pro  uly_lsf, spectrum, cmp, width, step, $
              FILE_OUT=file_out, PLOT=plot, CLEAN=clean, $
              KMOMENT=kmoment, DG=dg, DL=dl, KFIX=kfix, KPEN=kpen, $
              MDEGREE=mdegree, ADEGREE=adegree, VELSCALE=velscale, $
              QUIET=quiet   

if n_elements(spectrum) eq 0 or n_elements(cmp) eq 0 or $
   n_elements(width) eq 0 or n_elements(step) eq 0 then $
  message, 'Usage: ULY_LSF, <spectrum>, <cmp>, <width>, <step>, [optional keyw]'

if n_elements(file_out) ne 1 then lsf_file = 'output.lsf' $
else lsf_file = file_out

; ----------------------------------------------------------------------------
; Read the observed spectrum from which to determine the LSF
if uly_spect_get(spectrum, /VALID) eq 1 then begin
   SignalLog = uly_spect_alloc(SPECT=spectrum)
endif else begin
   SignalLog = uly_spect_read (spectrum, QUIET=quiet) 
endelse

if size(*SignalLog.data, /N_DIM) gt 1 then $
  message, /INFO, 'The file array to process has ' + $
  strtrim(size(*SignalLog.data, /N_DIM),2) + ' dimensions ... we use only the first spectrum'

SignalLog = uly_spect_extract(SignalLog, ONED=0, STATUS=status, /OVER)
if status ne 0 then message, 'Could not extract 1D spectrum'

SignalLog = uly_spect_logrebin(SignalLog, velscale, /OVER)

possig = where((*SignalLog.data) gt 0, cnt)
if cnt eq 0 then $
  message, /INFO, 'The signal is everywhere null or negative in the spectrum to analyse'

velscale = SignalLog.step * 299792.458d
wr = uly_spect_get(SignalLog, /WAVERANGE)


; ----------------------------------------------------------------------------
; Read the LSF reference cmp (used as an LSF standard)
status = uly_fit_init(cmp, WAVERANGE=wr, VELSCALE=velscale, QUIET=quiet)
if status ne 0 then begin
    message, 'Fit initialization failed, abort', /CONT
    return
endif

model_range = exp(cmp[0].start + [0.5,cmp[0].npix-1.5] * cmp[0].step)

; re-cut the input spectrum to the range of the model
SignalLog = uly_spect_extract(SignalLog, WAVERANGE=model_range, /OVER)
wr = uly_spect_get(SignalLog, /WAVERANGE)


; ----------------------------------------------------------------------------
npix = n_elements(*SignalLog.data)
lamrange = wr
nlsf = floor((lamrange[1]-lamrange[0]-width)/step) + 1

;loadct, 2, /SILENT 

if file_test(lsf_file) then file_delete, lsf_file

;openw, lun, lsf_file, BUFSIZE=0, /GET_LUN
;printf, lun, '# lsf_format=1', FORMAT='(a)'
;printf, lun, '# chi2, cz[km/s], e_cz, sigma[km/s], e_sigma, h3, e_h3, h4, e_h4, wl_c[Angstroms], Wstart, Wend', FORMAT='(a)'

undefine, sol_lsf
for nn = 0, nlsf-1 do begin
    minWL = lamrange[0]+step*nn
    maxWL = min([minWL+width, lamrange[1]])
;    centrWL= fix((maxWL-minWL)/2.+minWL)
    s = uly_spect_extract(SignalLog, WAVERANGE=[minWL,maxWL])
    if not keyword_set(quiet) then print, 'LSF chunk: ', nn
    ulyss, s, cmp, $
           KMOMENT=kmoment, DG=dg, DL=dl, KFIX=kfix, KPEN=kpen, $
           MDEGREE=mdegree, ADEGREE=adegree, QUIET=quiet, PLOT=plot, CLEAN=clean, SOL=sol
    uly_spect_free, s

;   write the results in the output file
    if size(sol,/TYPE) ne 8 then begin
; we have a pb there ... we expect a fixed set of wavelengths for each 
; LSF, it is in particular important when we measure a series of LSFs
; that we want to average (but maybe the solution is to change the way
; we read and average the lsfs
       print, 'hum ... what do we do when a fit failed?'
    endif else begin
       uly_solut_twrite, sol, lsf_file, /HEAD
       heap_free, sol
    endelse

;;   write the results in the output file
;    if size(sol,/TYPE) ne 8 then begin
;        if not keyword_set(quiet) then message, /INFO, 'The fit failed for LSF segment '+ strtrim(centrWL,2)
;        printf, lun, 0., 0., 0., 0., 0., 0., 0., 0., 0.,centrWL, minWL, maxWL, $
;          FORMAT='(12e15.7)'
;    endif else begin
;        printf, lun, sol.chisq, $
;          sol.losvd[0], sol.e_losvd[0], sol.losvd[1], sol.e_losvd[1],$
;          0.,0.,0.,0.,centrWL, minWL, maxWL, $
;          FORMAT='(12e15.7)'
;;        heap_free, sol
;    endelse
endfor


;close, lun
;free_lun, lun

if not keyword_set(quiet) then print, 'The results are saved in ', lsf_file

uly_spect_free, SignalLog 

end

;#############################################################################
