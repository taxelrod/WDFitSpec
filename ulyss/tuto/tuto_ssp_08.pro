star = uly_root+'/data/sdss_star.fits'
model = uly_tgm_extr(star,[4700.,3.,0.])
cmp = uly_tgm(TG=4700, LG=3., ZG=0)
uly_lsf, star, cmp, 600, 300, FIL='sdss_lsf.txt', /QUIET
uly_lsf_smooth, 'sdss_lsf.txt', 'sdss_lsfs.txt'
uly_lsf_plot, 'sdss_lsfs.txt'
