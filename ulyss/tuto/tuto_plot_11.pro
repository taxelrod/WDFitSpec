cmp = uly_ssp() 
chimap = uly_chimap(star, cmp, [[0,0],[0,1]], /QUIET)
plot_var = uly_plot_init(HALFTONECT=5)
uly_chimap_plot, chimap, /XLOG, XRANGE=[300., 18000.], PLOT_VAR=plot_var
uly_plot_save, 'filename.png', /PNG
