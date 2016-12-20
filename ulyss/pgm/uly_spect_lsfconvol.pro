;+
; NAME:           
;                 ULY_SPECT_LSFCONVOL
;
; PURPOSE:        
;                 Convolve a spectrum or a population model
;                 grid with a LSF.
;
; USAGE:          
;                 uly_spect_lsfconvol, <lsf_file>, <spec>
;
; DESCRIPTION:
;     Convolve the input <spect> with a line-spread function. <spec> can 
;     be 1 to 4D spectrum, like a population model grid.
;
;     The line-spread function, LSF, describes the intrinsic broadening 
;     of a spectrum and its variation with wavelength. It can be 
;     obtained with the program ULY_LSF.
;
;     The LSF is given as a series of Gauss or Gauss-Hermite expansions 
;     representing the instrumental broadening at different wavelengths.
;     To inject this LSF in <spec>, first <spec> is convolved
;     with each of these developments, and then the convolved spectra are
;     interpolated in wavelength to restitute the LSF at each wavelength.
;
; ARGUMENTS:
;   <lsf_file>:   Name of the file with the LSF it should contain 5
;                 columns with the parameters of the LOSVD(v, s, h3,
;                 h4) and the wavelength (usually the central wl were
;                 the fit for determining the LOSVD parameters was obtained.) 
;
;   <spec>:       Spectrum or a population model grid which *has to be logrebin*
;                 (can use ULY_SPECT_READ).
;
; EXAMPLE:   
;   1.Spectrum case:
;      First extract template at given Teff, Logg and [Fe/H], find the relative
;      change of the LSF between this template and the observed
;      spectrum (cflib_114642.fits), and finally inject it, which
;      means convolve this template with the derived LSF.
;
;          obs = uly_root+'/data/cflib_114642.fits'
;          template=uly_tgm_extr(obs,[6400,4.,0.])
;          cmp = uly_star(template)
;          uly_lsf, obs, cmp, 400, 200, FIL='lsfeg1.txt'
;
;          uly_spect_lsfconvol, 'lsfeg1.txt', template
;
;   2.Population model grid case:
;      Interpolate at given age/met; Find the LSF of Vazdekis/Miles
;      relative to PHR/Elodie3.1; read a model (PHR); and convolve
;      it with the derived LSF.
;
;         obs = uly_root+'/data/VazMiles_z-0.40t07.94.fits'
;         model = uly_ssp_extr(obs, uly_root+'/models/PHR_Elodie31.fits', $
;                              [6000., -0.4, 0], SIGMA=30)
;         cmp = uly_star(model)
;         uly_lsf, obs, cmp, 200, 100, FIL='lsfeg2.txt'
;
;         grid = uly_ssp_read(uly_root+'/models/PHR_Elodie31.fits',VELSCALE=30)
;         
;         uly_spect_lsfconvol, 'lsfeg2.txt', grid  
;         
;      or do the lsfconvol inside ULY_SSP_READ with giving LSF file           
;         newgrid = uly_ssp_read(uly_root+'/models/PHR_Elodie31.fits', $
;                                VELSCALE=30, LSF='lsfeg2.txt')
;    
; HISTORY: 
;                 Yue WU, 21/08/2008, 
;                 imitated from ULY_SPECT_LOSVDCONVOL made by Ph. Prugniel
;
;-
; CATEGORY:       ULY
; -----------------------------------------------------------------
function uly_lsfconvol_1D, lsf_file, data, start, step, QUIET=quiet

compile_opt idl2, hidden

velscale = 299792.458d * step

if file_test(lsf_file) eq 0 then begin
    message, /CONT, 'LSF file not found '+lsf_file
    return, 0
endif

rlsf = uly_lsf_read(lsf_file)


; Check if the LSF file was read sucessfully
nlsf = n_elements(rlsf.v)
if nlsf eq 0 then begin
    message, /CONT, 'Invalid LSF file '+lsf_file
    return, 0
endif

s = size(data)
npx = s[1]

alsf=replicate({wl:1d,conv:dblarr(npx)},nlsf) 
alsf.wl = rlsf.centrWL

for k=0,nlsf-1 do begin
    uly_slit, rlsf.v[k], rlsf.sig[k], rlsf.h3[k], rlsf.h4[k], velscale, lsf
    lsf = lsf/total(lsf)    ;normalising lsf    
    alsf[k].conv = convol(data, lsf, /EDGE_TRUNCATE) 
endfor

wl = exp([start + step*lindgen(npx)])

;check if the lsf.wl is in the wl range of the spec
limits = where(alsf.wl ge wl[0] and alsf.wl le wl[n_elements(wl)-1])
alsf = alsf[min(limits, MAX=mxl):mxl]
nlsf = n_elements(alsf.wl)

px = findgen(n_elements(wl))

alsf.wl = (alog(alsf.wl)-start) / step
minPX = min((alsf.wl)[0:n_elements(alsf.wl)-1], MAX=maxPX)

n1 = WHERE(px lt minPX)
n3 = WHERE(px ge maxPX)

(data)[n1] = alsf[0].conv[n1]
(data)[n3] = alsf[nlsf-1].conv[n3]

for k = 0, nlsf-2 do begin
    n2 = WHERE((px ge alsf[k].wl) and (px lt alsf[k+1].wl), nn2)
    ff = rebin((n2-alsf[k].wl) / (alsf[k+1].wl-alsf[k].wl ), nn2)
    (data)[n2] = alsf[k].conv[n2]  + $
      ff * (alsf[k+1].conv[n2]-alsf[k].conv[n2])
endfor

return, data
end

;------------------------------------------------------------------
function uly_lsfconvol_2D, lsf_file, data, start, step, QUIET=quiet

compile_opt idl2, hidden

velscale = 299792.458d * step

rlsf = uly_lsf_read(lsf_file) ; read the lsf file


; Check if the LSF file was read sucessfully
nlsf = n_elements(rlsf.v)
if nlsf eq 0 then begin
    message, /CONT, 'Invalid LSF file '+lsf_file
    return, 0
endif

s = size(data)
npx = s[1]
n2d = s[2]

alsf = replicate({wl:1d,conv:dblarr(npx,n2d)},nlsf) 
alsf.wl = rlsf.centrWL

for k=0,nlsf-1 do begin
    uly_slit, rlsf.v[k], rlsf.sig[k], rlsf.h3[k], rlsf.h4[k], velscale, lsf
    lsf = lsf/total(lsf)        ;normalizing lsf    
    alsf[k].conv = convol(data, lsf, /EDGE_TRUNCATE) 
endfor

wl = exp([start + step*lindgen(npx)])

;check if the lsf.wl is in the wl range of the spec
limits = where(alsf.wl ge wl[0] and alsf.wl le wl[n_elements(wl)-1])
alsf = alsf[min(limits, MAX=mxl):mxl]
nlsf = n_elements(alsf.wl)

px = findgen(n_elements(wl))

alsf.wl = (alog(alsf.wl)-start) / step
minPX = min((alsf.wl)[0:n_elements(alsf.wl)-1], MAX=maxPX)

n1 = WHERE(px lt minPX)
n3 = WHERE(px ge maxPX)

(data)[n1,*] = alsf[0].conv[n1,*]
(data)[n3,*] = alsf[nlsf-1].conv[n3,*]

for k = 0, nlsf-2 do begin
    n2 = WHERE((px ge alsf[k].wl) and (px lt alsf[k+1].wl), nn2)
    ff = rebin((n2-alsf[k].wl) / (alsf[k+1].wl-alsf[k].wl ), nn2, n2d)
    (data)[n2,*] = alsf[k].conv[n2,*]  + $
      ff * (alsf[k+1].conv[n2,*]-alsf[k].conv[n2,*])
endfor

return, data
end

;---------------------------------------------------
pro uly_spect_lsfconvol, lsf_file, spec, QUIET=quiet

;; input checks
if n_elements(lsf_file) eq 0 then begin 
    print, 'Usage: uly_spect_lsfconvol, <lsf_file>, <spec>'
    message, 'Missing argument <lsf_file>', /INFO
    return
endif

if file_test(lsf_file) ne 1 then begin
    f = ''
    if size(lsf_file, /TYPE) eq 7 then f = ' ('+lsf_file+')'
    message, 'No valid input LSF file: '+f, /INFO
    return
endif

if n_elements(spec) eq 0 then begin
    print, 'Usage: uly_spect_lsfconvol, <lsf_file>, <spec>'
    message, 'Missing argument <spec>'
    return
endif

if size(spec, /TYPE) ne 8 then begin
   spec = uly_spect_read(spec, /QUIET)
endif 

s = size(*spec.data)
type = ''
if size(spec, /TYPE) eq 8 then begin
    case s[0] of
        0: message, 'No valid input spectrum data' 
        1: type = '1D'
        2: type = '2D'
        else: type = '>2D'
    endcase
endif

if spec.sampling ne 1 then $
  message, 'The input spectrum should be sampled in log-wave'

; switch to the appropriate reading routine
if type eq '1D' then begin
    *spec.data = uly_lsfconvol_1D(lsf_file, *spec.data, spec.start, spec.step, QUIET=quiet)
    return
endif

if type eq '2D' then begin
    *spec.data = uly_lsfconvol_2D(lsf_file, *spec.data, spec.start ,spec.step, QUIET=quiet)
    return
endif

if type eq '>2D' then begin
    ndim = size(*spec.data, /N_DIM)
    data = *spec.data
    start = spec.start
    step = spec.step    
    dim = size(*spec.data, /DIM) 
    if ndim ge 3 then begin     ; > 2D spectrum  
        tot = 1
        for i = 1, ndim-1 do tot *= dim[i]
        data = reform(data, dim[0], tot, /OVER) 
    endif
    (*spec.data) = reform(uly_lsfconvol_2D(lsf_file, data, start, step, $
                                           /QUIET), dim, /OVER)
endif

end
;--- end ----------------------------------------------------------------------
