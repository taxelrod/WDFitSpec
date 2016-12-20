save_path = !path
save_root = uly_root
!path += ':' + uly_root + '/tuto'
uly_root = '.'

;------------------------------------------------------------------------------
; Basic examples
message, /INFO, 'TUTO BASE'
@ tuto_base_01
@ tuto_base_02
spawn, 'cat tuto/path.pro tuto/tuto_base_01.pro tuto/tuto_base_02.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_base_scr01.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

@ tuto_base_03
write_jpeg, "doc/tuto_base_img01.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
spawn, 'cat tuto/path.pro tuto/tuto_base_01.pro tuto/tuto_base_03.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_base_scr02.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

@ tuto_base_04
spawn, 'cat tuto/path.pro tuto/tuto_base_04.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_base_scr03.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

if file_test('tuto_base.res') then file_delete, 'tuto_base.res', 'tuto_base.fits'
message, /INFO, 'TUTO BASE 5'
@ tuto_base_05

@ tuto_base_06
write_jpeg, "doc/tuto_base_img02.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
spawn, 'cat tuto/path.pro tuto/tuto_base_06.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_base_scr04.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun
file_delete, 'tuto_base.res', 'tuto_base.fits'

@ tuto_base_07
write_jpeg, "doc/tuto_base_img03.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
spawn, 'cat tuto/path.pro tuto/tuto_base_07.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_base_scr05.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

@ tuto_base_08
write_jpeg, "doc/tuto_base_img04.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
spawn, 'cat tuto/path.pro tuto/tuto_base_08.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_base_scr06.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

@ tuto_base_09

@ tuto_base_10
file_delete, 'tuto_lsf.txt', 'tuto_lsf_s.txt'

;------------------------------------------------------------------------------
; tuto_chi
message, /INFO, 'TUTO CHI'
@ tuto_chi_01
@ tuto_chi_02
spawn, 'cat tuto/path.pro tuto/tuto_chi_01.pro tuto/tuto_chi_02.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_chi_scr01.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

if file_test('tuto_chi.res') then file_delete, 'tuto_chi.res'
if file_test('tuto_chi.fits') then file_delete, 'tuto_chi.fits'
@ tuto_chi_03   ; Monte-Carlo simulation (default S/N)
write_jpeg, "doc/tuto_chi_img01.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1

@ tuto_chi_04   ; plot the result of M-C simulations
write_jpeg, "doc/tuto_chi_img02.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1

if file_test('tuto_chi.res') then file_delete, 'tuto_chi.res'
if file_test('tuto_chi.fits') then file_delete, 'tuto_chi.fits'
@ tuto_chi_05   ; M-Csimulations S/N=15
write_jpeg, "doc/tuto_chi_img03.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1

@ tuto_chi_06   ; chi2 maps
wset, 0
write_jpeg, "doc/tuto_chi_img04.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
wset, 1
write_jpeg, "doc/tuto_chi_img05.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
wset, 0

@ tuto_chi_07  ; plot chi2 maps (overlays)
write_jpeg, "doc/tuto_chi_img06.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1

; file_delete, 'tuto_chi.res', 'tuto_chi.fits'

@ tuto_chi_08
write_jpeg, "doc/tuto_chi_img07.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1

file_delete, 'absolute.res', 'absolute.fits', 'tuto_chi.res'
;------------------------------------------------------------------------------
; tuto tgm
;@ tuto_tgm_01
;
;@ tuto_tgm_02
;write_jpeg, "doc/tuto_tgm_img01.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
;spawn, 'cat tuto/path.pro tuto/tuto_tgm_01.pro tuto/tuto_tgm_02.pro
;|idl -arg quiet', result, error
;openw, lun, 'doc/tuto_tgm_scr01.txt', BUFSIZE=0, /GET_LUN
;printf, lun, result, F='(A)'
;free_lun, lun 
;file_delete, 'tgm_test_1.res','tgm_test_1.fits'
;heap_free, cmp
;
;@ tuto_tgm_03
;write_jpeg, "doc/tuto_tgm_img02.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
;spawn, 'cat tuto/path.pro tuto/tuto_tgm_03.pro |idl -arg quiet', result, error
;
;@ tuto_tgm_04
;write_jpeg, "doc/tuto_tgm_img03.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
;spawn, 'cat tuto/path.pro tuto/tuto_tgm_01.pro tuto/tuto_tgm_04.pro
;|idl -arg quiet', result, error
;openw, lun, 'doc/tuto_tgm_scr02.txt', BUFSIZE=0, /GET_LUN
;printf, lun, result, F='(A)'
;free_lun, lun
;file_delete, 'tgm_test_2','tgm_test_2s','tgm_test_3.res','tgm_test_3.fits'
;heap_free, cmp
;
;@ tuto_tgm_05
;write_jpeg, "doc/tuto_tgm_img04.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
;spawn, 'cat tuto/path.pro tuto/tuto_tgm_01.pro tuto/tuto_tgm_05.pro
;|idl -arg quiet', result, error
;heap_free, cmp
;
;@ tuto_tgm_06
;write_jpeg, "doc/tuto_tgm_img05.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
;spawn, 'cat tuto/path.pro tuto/tuto_tgm_01.pro tuto/tuto_tgm_06.pro
;|idl -arg quiet', result, error
;heap_free, cmp
;
;@ tuto_tgm_07
;write_jpeg, "doc/tuto_tgm_img06.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
;spawn, 'cat tuto/path.pro tuto/tuto_tgm_01.pro tuto/tuto_tgm_07.pro
;|idl -arg quiet', result, error
;heap_free, cmp
;
;@ tuto_tgm_08
;write_jpeg, "doc/tuto_tgm_img07.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
;spawn, 'cat tuto/path.pro tuto/tuto_tgm_08.pro |idl -arg quiet', result, error
;heap_free, cmp
;file_delete, 'tgmcv*'
;------------------------------------------------------------------------------------
; tuto_ssp
message, /INFO, 'TUTO SSP'

if file_test('m67_phr.res') then file_delete, 'm67_phr.res', 'm67_phr.fits'
if file_test('m67_vaz.res') then file_delete, 'm67_vaz.res', 'm67_vaz.fits'
if file_test('lsf_m67_phr.txt') then file_delete, 'lsf_m67_phr.txt'
@ tuto_ssp_01   ; make a simple ssp with phr fit to prepare LSF analysis
if file_test('m67_phr.res') then file_delete, 'm67_phr.res', 'm67_phr.fits'
spawn, 'cat tuto/path.pro tuto/tuto_ssp_01.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_ssp_scr01.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

@ tuto_ssp_02   ; analyse the LSF, and make the fit with it
write_jpeg, "doc/tuto_ssp_img02.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
if file_test('m67_phr_2.res') then file_delete, 'm67_phr_2.res', 'm67_phr_2.fits'
spawn, 'cat tuto/path.pro tuto/tuto_ssp_02.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_ssp_scr02.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

@ tuto_ssp_03   ; make a simple ssp with vaz fit to prepare LSF analysis
if file_test('m67_vaz.res') then file_delete, 'm67_vaz.res', 'm67_vaz.fits'
spawn, 'cat tuto/path.pro tuto/tuto_ssp_03.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_ssp_scr03.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

@ tuto_ssp_04   ; inject lsf and analyze with vaz
write_jpeg, "doc/tuto_ssp_img04.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
if file_test('m67_vaz.res') then file_delete, 'm67_vaz.res', 'm67_vaz.fits'
spawn, 'cat tuto/path.pro tuto/tuto_ssp_04.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_ssp_scr04.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

@ tuto_ssp_05   ; present two plots
wset, 0
write_jpeg, "doc/tuto_ssp_img05a.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
wset, 1
write_jpeg, "doc/tuto_ssp_img05b.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1

message, /INFO, 'TUTO SSP_6'
@ tuto_ssp_06   ; quick fit to the L* SDSS gal
write_jpeg, "doc/tuto_ssp_img06.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
spawn, 'cat tuto/path.pro tuto/tuto_ssp_06.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_ssp_scr06.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

if file_test('m67_phr.res') then file_delete, 'm67_phr.res', 'm67_phr.fits'
if file_test('m67_phr_2.res') then file_delete, 'm67_phr_2.res', 'm67_phr_2.fits'
if file_test('m67_vaz.res') then file_delete, 'm67_vaz.res' 'm67_vaz.fits'
if file_test('m67_vaz_2.res') then file_delete, 'm67_vaz_2.res', 'm67_vaz_2.fits'
file_delete, 'lsf_m67_vaz.txt', 'lsfs_m67_phr.txt', /ALLOW_NON, /QUIET

message, /INFO, 'TUTO SSP_7'
@ tuto_ssp_07   ; TGM fit to the SDSS star
spawn, 'cat tuto/path.pro tuto/tuto_ssp_07.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_ssp_scr07.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun
@ tuto_ssp_08   ; LSF fit to SDSS, and plot
write_jpeg, "doc/tuto_ssp_img08.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1

message, /INFO, 'TUTO SSP_9'
@ tuto_ssp_09   ; inject in the SSP grid and save the grid

if file_test('sdss.res') then file_delete, 'sdss.res', 'sdss.fits'
@ tuto_ssp_10   ; SSP fit to the SDSS Lstar with the proper LSF
wset, 0
write_jpeg, "doc/tuto_ssp_img10a.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
wset, 1
write_jpeg, "doc/tuto_ssp_img10b.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
spawn, 'cat tuto/path.pro tuto/tuto_ssp_10.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_ssp_scr10.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

; Note: thehre is an error message in toto_ssp_11:
; Program caused arithmetic error: Floating illegal operand
; It comes from the fact that the fraction of the second burst is 0.
message, /INFO, 'TUTO SSP_11'
file_delete, 'sdss.res', 'sdss.fits'
@ tuto_ssp_11   ; SSP fit to the SDSS Lstar with the proper LSF
wset, 0
write_jpeg, "doc/tuto_ssp_img11a.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
wset, 1
write_jpeg, "doc/tuto_ssp_img11b.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
spawn, 'cat tuto/path.pro tuto/tuto_ssp_11.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_ssp_scr11.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun
file_delete, 'sdss.res', 'sdss.fits', 'sdss_lsf*.txt', 'PHR_Elodie31_SDSS.fits', /ALLOW_NON, /QUIET


;------------------------------------------------------------------------------
; creation of a new cmp
message, /INFO, 'TUTO CMP'

@ tuto_cmp_01

@ tuto_cmp_02
write_jpeg, "doc/tuto_cmp_img01.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
spawn, 'cat tuto/path.pro tuto/tuto_cmp_01.pro tuto/tuto_cmp_02.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_cmp_scr02.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

@ tuto_cmp_03
wset, 0
write_jpeg, "doc/tuto_cmp_img02.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
wset, 1
write_jpeg, "doc/tuto_cmp_img03.jpg", tvrd(TRUE=1), QUALITY=100, TRUE=1
spawn, 'cat tuto/path.pro tuto/tuto_cmp_01.pro tuto/tuto_cmp_03.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_cmp_scr03.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

;------------------------------------------------------------------------------
; Modify graphics

message, /INFO, 'TUTO PLOT'
@ tuto_plot_01
write_png, "doc/tuto_plot_img01.png", tvrd(TRUE=1)

@ tuto_plot_02

@ tuto_plot_03

@ tuto_plot_04

@ tuto_plot_05

@ tuto_plot_06
wset, 0
write_png, "doc/tuto_plot_img03.png", tvrd(TRUE=1)
wset, 1
write_png, "doc/tuto_plot_img04.png", tvrd(TRUE=1)
wset, 0

@ tuto_plot_07

@ tuto_plot_08
spawn, 'cat tuto/path.pro tuto/tuto_plot_08.pro |idl -arg quiet', result, error
openw, lun, 'doc/tuto_plot_scr01.txt', BUFSIZE=0, /GET_LUN
printf, lun, result, F='(A)'
free_lun, lun

@ tuto_plot_09
write_png, "doc/tuto_plot_img05.png", tvrd(TRUE=1)

@ tuto_plot_10

@ tuto_plot_11
write_png, "doc/tuto_plot_img06.png", tvrd(TRUE=1)

file_delete, 'filename.png', 'bw_settings.sav'

;------------------------------------------------------------------------------
; restore the state
!path = save_path
uly_root = save_root
;--- end ---------------------------------------------------------------------
