;+
; NAME:           
;                 ULY_LSF_SMOOTH
;
; PURPOSE:        
;                 Smooth a LSF file
;
; USAGE:       
;                 uly_lsf_smooth, <lsf_file>, <outfile>, [VDEG=<vdeg>], [SDEG=<sdeg>]
;
; ARGUMENTS:
;   <lsf_file>:   Name of the file containing the LSF
;
;   <outfile>:    Name of the file containing the smoothed LSF
;
; KEYWORDS:
;    VDEG:        degree of polynomial, used to smooth the velocity, default=3
;    SDEG:        degree of polynomial, used to smooth the dispersion, default=3
;    H3DEG:       degree of polynomial, used to smooth h3, default=3
;    H4DEG:       degree of polynomial, used to smooth h4, default=3
;
; DESCRIPTION:
;     Fit the LSF determined by ULY_LSF (cz vs. wave) and (sigma vs. wave)
;     respectively with degrees 3 and 2 polynomials and return the
;     fitted values as a 'smoothed' LSF.
;
;     As it is, this program is obviously a stub on which to build the 
;     smoothing. We may make a more elaborated version in the future.
;
; DEPENDENCE:     ULY_LSF_PLOT, READCOL (by W. Landsman, astrolib)
;
; AUTHOR:         Mina Koleva, Philippe Prugniel
;
;-
; CATEGORY:       ULY
;------------------------------------------------------------------------------
pro uly_lsf_smooth, lsf_file, outfile, VDEG=vdeg, SDEG=sdeg, H3DEG=h3deg, H4DEG=h4deg

if not keyword_set(vdeg) then vdeg=3
if not keyword_set(sdeg) then sdeg=3
if not keyword_set(h3deg) then h3deg=3
if not keyword_set(h4deg) then h4deg=3


lsf = uly_solut_tread(lsf_file, STATUS=status)
if status ne 0 then begin
   message, 'Could not read the LSF file (shall be in solut format)'
endif

if (size(lsf.losvd, /DIM))[0] ge 2 then begin
   v = reform(lsf.losvd[0,*])
   sig = reform(lsf.losvd[1,*])
endif
if (size(lsf.losvd, /DIM))[0] ge 3 then h3 = lsf.losvd[2,*] else h3 = 0*v
if (size(lsf.losvd, /DIM))[0] ge 4 then h4 = lsf.losvd[3,*] else h4 = 0*v
centrWL = 0*v
for k=0,n_elements(lsf)-1 do begin
   wave_min = sxpar(*(lsf[k]).hdr, 'WAVE_MIN')
   wave_max = sxpar(*(lsf[k]).hdr, 'WAVE_MAX')
   centrWL[k] = (wave_min + wave_max) / 2
endfor

; Check if the LSF file was read sucessfully
nlsf = n_elements(v)
if nlsf eq 0 then begin
    message, 'Invalid LSF file '+lsf_file
endif

fitfit = poly_fit(centrwl, v, fix(vdeg), YFIT=vfit)
fitfit = poly_fit(centrwl, sig, fix(sdeg), YFIT=sfit)
fitfit = poly_fit(centrwl, h3, fix(h3deg), YFIT=h3fit)
fitfit = poly_fit(centrwl, h4, fix(h4deg), YFIT=h4fit)

lsf.losvd[0,*] = reform(vfit,1,nlsf)
lsf.losvd[1,*] = reform(sfit,1,nlsf)
if (size(lsf.losvd, /DIM))[0] ge 3 then lsf.losvd[2,*] = reform(h3fit,1,nlsf)
if (size(lsf.losvd, /DIM))[0] ge 4 then lsf.losvd[3,*] = reform(h4fit,1,nlsf)

if file_test(outfile) then file_delete, outfile
for k=0,n_elements(lsf)-1 do uly_solut_twrite, lsf[k], outfile, HEAD=*(lsf[k]).hdr

;; note we can not use forprint which will substitute the next 7 raws
;; because "forprint" is adding by default an extension (.prt)  
;forprint, TEXTOUT=outfile, vfit, sfit, h3, h4, centrWL, /SILENT, /noCOMMENT

;openw, lun, outfile, BUFSIZE=0, /GET_LUN
;printf, lun, '# lsf_format=1', FORMAT='(a)'
;printf, lun, '# Smoothed LSF', FORMAT='(a)'
;printf, lun, '# chi2, cz[km/s], e_cz, sigma[km/s], e_sigma, h3, e_h3, h4, e_h4, wl_c[Angstroms], Wstart, Wend', FORMAT='(a)'
;for k = 0, n_elements(v)-1 do begin
;    printf, lun, chi2[k], vfit[k], e_v[k], sfit[k], e_sig[k], h3[k], e_h3[k], h4[k], e_h4[k], centrWL[k], wstart[k], wend[k], $
;      FORMAT='(12e14.7)'
;endfor
;close, lun
;free_lun, lun


uly_lsf_plot, lsf_file, PSYM=4, LINETHICK=3
uly_lsf_plot, outfile, LINECOLOR='RED', LINETHICK=3, /OVER 

end
;---end---------------------------------------------------
