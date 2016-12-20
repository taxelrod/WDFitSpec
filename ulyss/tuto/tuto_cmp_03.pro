g = uly_spect_read(uly_root+'/data/VazMiles_z-0.40t07.94.fits')
g = uly_spect_logrebin(g, /OVERWRITE)
tau = -3d
p = exp(findgen(n_elements(*g.data))*g.step*tau)
*g.data += 0.5 * mean(*g.data)/mean(p) * p 

window, 0
uly_spect_plot, galaxy, LINECOLOR='Black'
uly_spect_plot, g, /OVERPLOT

window, 1
cmp = [power, uly_ssp(MODEL_FILE=uly_root+'/models/Vaz_Miles.fits')]
ulyss, g, cmp, KMOMENT=0, MD=0, /PLOT
