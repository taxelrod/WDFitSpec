;+
; NAME:
;                 ULY_TEST
;
; PURPOSE:
;                 Test the ULySS package
;
; USAGE:
;		  uly_test [, /REGRESS][, CRASH=crash][, STATUS=status]
;
; KEYWORDS:
;   /REGRESS:     
;		  If this keyword is set, perform regression tests,
;                 otherwise crash tests are made.
;
;   CRASH:        
;		  0 or not set, perform all the tests
;                 n, perform test #n (there are 12 test available)
;
;   STATUS:       
;		  0 if the tests succeed, "n" if test #n failed
;                                                            
; DESCRIPTION:
;   Perform tests of the ULySS package.
;  
;   At present only crash tests are executed. They test all the
;   routines of the package in at minimum one configuration. The goal
;   of these tests is to detect bugs in the packages introduced during
;   the development.
;
; AUTHOR:
;                 Mina Koleva, Philippe Prugniel & Yue Wu
;
; HISTORY:
;                 14/02/2008
;-
; CATEGORY:    ULY
pro uly_test, REGRESS=regress, CRASH=crash, STATUS=status

common uly_path, uly_root

status = -1

if n_elements(uly_root) ne 0 then begin
    if file_test(uly_root+'/data/VazMiles_z-0.40t07.94.fits') eq 1 then $
       datadir = uly_root+'/data'
endif

if n_elements(datadir) ne 0 then begin
  galaxy = datadir + '/VazMiles_z-0.40t07.94.fits' 
  star = datadir + '/cflib_114642.fits'  ; sptype = 'F6V'
endif else $
  message, 'Could not find the test files'

modeldir = uly_root+'/models'

if keyword_set(regress) then begin
    print, 'No regression tests yet !!!!'
endif

if n_elements(crash) eq 0 and n_elements(datadir) ne 0 then crash=0

if crash eq 0 or crash eq 1 then begin
    status = 1
    plot_var = uly_plot_init(LINECOLOR='Black')
    uly_spect_plot, galaxy, PLOT_VAR=plot_var, /QUIET
    s1 = uly_spect_read(galaxy, /QUIET)
    s1 = uly_spect_losvdconvol(s1, 10000., 300., 0., 0., /OVERWRITE)
    uly_spect_plot, s1, /OVERPLOT, LINECOLOR='Blue'
    uly_spect_free, s1
    print, 'Crash test 1 succeeded ...'
endif

if crash eq 0 or crash eq 2 then begin
    status = 2
    model = modeldir+'/sun.fits'
    cmp = uly_star(model)
    ulyss, star, cmp, FILE_OUT='crashtest2', /PLOT, /QUIET
    heap_free, cmp
    file_delete, 'crashtest2.res', 'crashtest2.fits'
    print, 'Crash test 2 succeeded ...'
endif

if crash eq 0 or crash eq 3 then begin
    status = 3
    cmp = uly_ssp()
    ulyss, galaxy, cmp, /CLEAN, /QUIET, FILE_OUT='crashtest3'
    sp = uly_spect_read('crashtest3.fits', /QUIET)
    s = uly_solut_tread('crashtest3')
    heap_free, s
    heap_free, sp
    heap_free, cmp
    uly_solut_tprint, 'crashtest3.res'
    uly_solut_splot, 'crashtest3.fits', SMOOTH_RES=50
    file_delete, 'crashtest3.res', 'crashtest3.fits'
    print, 'Crash test 3 succeeded ...'
endif

if crash eq 0 or crash eq 4 then begin
    status = 4
    cmp = uly_ssp(MODEL_FILE=modeldir+'/Vaz_Miles.fits', AG=[1000d, 8000d])
    ulyss, galaxy, cmp,  KFIX=[0,1], DG=0., /PLOT, /QUIET, FILE_OUT='crashtest4'
    file_delete, 'crashtest4.res', 'crashtest4.fits'
    heap_free, cmp
    print, 'Crash test 4 succeeded ...'
endif

if crash eq 0 or crash eq 5 then begin
    status = 5
    if file_test('crashtest5.res') then file_delete, 'crashtest5.res'
    cmp = uly_ssp()
    spect = uly_spect_read(galaxy, /QUIET)
    if n_elements(*spect.err) eq 0 then $
       *spect.err = mean(*spect.data)/100. + *spect.data * 0
    ulyss, spect, cmp, NSIMUL=200, /PLOT, /QUIET, FILE_OUT='crashtest5'
    file_delete, /ALLOW_NO, 'crashtest5.res'
    heap_free, spect
    heap_free, cmp
    print, 'Crash test 5 succeeded ...'
endif

if crash eq 0 or crash eq 6 then begin
    status = 6
    cmp = uly_ssp()
    range = dindgen(20)*(12000.-500.)/19.+500
    map = uly_chimap(galaxy, cmp,[[0,0],[0,1]], /QUIET, $
                     FILE_OUT='crashtest6', $
                     XNODE=alog(range))
    uly_chimap_plot, map, /QUIET
    uly_chimap_plot, 'crashtest6', /QUIET, /XLOG
    file_delete, 'crashtest6'
    heap_free, cmp
    print, 'Crash test 6 succeeded ...'
endif

if crash eq 0 or crash eq 7 then begin
    status = 7
    model = uly_ssp_extr(galaxy, modeldir+'/PHR_Elodie31.fits' ,[6000.,-0.4,0.], SIG=30, /QUIET)
    cmp = uly_star(model)
    uly_lsf, galaxy, cmp, 400, 200, FIL='crashtest7', /QUIET
    heap_free, cmp
    uly_lsf_plot, 'crashtest7', YST=3, PSYM=4
    uly_lsf_smooth, 'crashtest7', 'crashtest7s'
    cmp = uly_ssp(LSF='crashtest7s')
    ulyss, galaxy, cmp, /PLOT, /QUIET, FILE_OUT='crashtest7'
    file_delete, /ALLOW_NO, 'crashtest7', 'crashtest7s', 'crashtest7.res', 'crashtest7.fits'
    print, 'Crash test 7 succeeded ...'
endif

if crash eq 0 or crash eq 8 then begin
    status = 8
    grid = uly_ssp_read(modeldir+'/PHR_Elodie31.fits', VELSCALE=30., /QUIET)
    uly_ssp_write, grid, 'crashtest8'
    file_delete, 'crashtest8'
    print, 'Crash test 8 succeeded ...'
endif

if crash eq 0 or crash eq 9 then begin
    status = 9
    cmp = uly_tgm()
    ; lim=3900,lmax=6800 according to elodie32 model wavelength range
    ulyss, star, cmp, LMIN=3900.0, LMAX=6800.0, /CLEAN, /PLOT, /QUIET, FILE_OUT='crashtest9'
    if file_test('crashtest9.res') then file_delete, 'crashtest9.res'
    if file_test('crashtest9.fits') then file_delete, 'crashtest9.fits'
    heap_free, cmp
    print, 'Crash test 9 succeeded ...'
endif

if crash eq 0 or crash eq 10 then begin
    status = 10
    cmp0 = uly_ssp()
    cmp1 = uly_star(modeldir+'/sun.fits')
    ulyss, galaxy, [cmp0,cmp1], /PLOT, /QUIET, FIL='crashtest10'
    heap_free, [cmp0, cmp1]
    file_delete, 'crashtest10.res', 'crashtest10.fits'
    print, 'Crash test 10 succeeded ...'
endif

if crash eq 0 or crash eq 11 then begin
    status = 11
    model = uly_tgm_extr(star,[6375.,3.94,-0.17])
    cmp = uly_star(model)
    uly_lsf, star, cmp, 800, 400, FIL='crashtest11', KMOMENT=4, /QUIET
    uly_lsf_plot, 'crashtest11', YST=3, PSYM=4
    file_delete, 'crashtest11'
    print, 'Crash test 11 succeeded ...'
endif


if crash eq 0 or crash eq 12 then begin    ; test AD
    status = 12
    cmp = uly_ssp()
    ulyss, galaxy, cmp, AD=2, FILE_OUT='crashtest12', /QUIET, MD=0
    if file_test('crashtest12.res') then file_delete, 'crashtest12.res'
    if file_test('crashtest12.fits') then file_delete, 'crashtest12.fits'
    print, 'Crash test 12 succeeded ...'
endif

if crash eq 0 or crash eq 13 then begin
;  print, 'Starting crash test 13'
;  add test for uly_solut_tplot
;  print, 'Crash test 13 succeeded ...'
endif

status = 0
print
print, 'All test succeeded!'
print 

return
end
