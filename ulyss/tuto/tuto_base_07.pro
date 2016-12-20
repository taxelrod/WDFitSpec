galaxy = uly_root+'/data/VazMiles_z-0.40t07.94.fits'
cmp1 = uly_ssp(AG=[1000.], ZG=[0], AL=[200.,3000.])
cmp2 = uly_ssp(AG=[10000.], ZG=[0], AL=[3000., 12000.])
cmp = [cmp1, cmp2]
ulyss, galaxy, cmp, /PLOT
