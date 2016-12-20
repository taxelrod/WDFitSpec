galaxy = uly_root+'/data/VazMiles_z-0.40t07.94.fits'
spectrum = uly_spect_read(galaxy)
specConv = uly_spect_losvdconvol(spectrum, 0., 30., 0, 0)
uly_spect_free, spectrum
ulyss, specConv, MODEL='models/PHR_Elodie31.fits', /QUIET, FILE_OUT='tuto_base'
