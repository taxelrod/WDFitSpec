chi1 = uly_chimap(spectrum, cmp, [[0,0],[0,1]])
chi2 = uly_chimap(spectrum, cmp, [[1,0],[1,1]])
window, 0
uly_chimap_plot, chi1, /XLOG, YRANGE=[-1.2, 0.3]
window, 1
uly_chimap_plot, chi2, /XLOG, YRANGE=[-1.2, 0.3]
