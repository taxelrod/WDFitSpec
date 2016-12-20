cmp = uly_ssp(MODEL='PHR_Elodie31_SDSS.fits')
ulyss, uly_root+'/data/sdss_lstar.fits', cmp, FILE='sdss', /QUIET
window, 0
uly_solut_splot, 'sdss', WAVERANGE=[4050, 4350]
window, 1
uly_solut_splot, 'sdss', WAVERANGE=[4800, 5200]
