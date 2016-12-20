grid=uly_ssp_read(uly_root+'/models/PHR_Elodie31.fits', SIGMA=50., VELSCALE=50.)
young = uly_ssp_interp(grid, [alog(500.), 0.])
old = uly_ssp_interp(grid, [alog(10000.), -1.])

spectrum = uly_spect_alloc(DATA=0.1*young+0.9*old, START=grid.start, STEP=grid.step, SAMPLING=1)
