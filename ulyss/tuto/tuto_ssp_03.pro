cmp_vaz = uly_ssp(MODEL=uly_root+'/models/Vaz_Miles.fits') ; Vazdekis/Miles
ulyss, uly_root+'/data/m67.fits', cmp_vaz, FILE='m67_vaz', /CLEAN, /QUIET, /PLO
uly_solut_tprint, 'm67_vaz'
