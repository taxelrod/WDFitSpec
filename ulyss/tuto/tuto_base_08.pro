star = uly_root+'/data/cflib_114642.fits'
model = uly_tgm_extr(star,[6400.,4.,0.]) ;  Teff=6400K log(g)=4 [Fe/H]=0
cmp = uly_star(model)
uly_lsf, star, cmp, 600, 300, FILE='tuto_lsf.txt', /QUIET
plotsym,0,1,/FILL
uly_lsf_plot, 'tuto_lsf.txt', YST=3, PSYM=8, CONNECT='Red'
