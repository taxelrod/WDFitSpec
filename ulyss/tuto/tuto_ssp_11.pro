cmp1 = uly_ssp(MODEL='PHR_Elodie31_SDSS.fits', AG=[500.])
cmp2 = uly_ssp(MODEL='PHR_Elodie31_SDSS.fits', AG=[2000.])
cmp3 = uly_ssp(MODEL='PHR_Elodie31_SDSS.fits', AG=[12000.])
cmp = [cmp1, cmp2, cmp3]
ulyss, uly_root+'/data/sdss_lstar.fits', cmp, FILE='sdss'
window, 0
uly_solut_splot, 'sdss', WAVERANGE=[4050, 4350]
window, 1
uly_solut_splot, 'sdss', WAVERANGE=[4800, 5200]
