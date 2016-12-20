gal=uly_root+'/data/VazMiles_z-0.40t07.94.fits'
cmp = uly_ssp(AG=[500d,1000d,2000d,4000d,8000d,15000d,18000d], ZG=[-2d,-1.5d,-1d,-0.5d,0d,0.3d])
ulyss, gal, cmp, /ALL_SOL, FIL='absolute'
uly_solut_tplot, 'absolute.res', XA=[0,0], YA=[0,1], XR=[300,25000], YR=[-1.6,0.5], /CV, /XLOG
