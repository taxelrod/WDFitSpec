cmp_peg = uly_ssp()           ; Pegase_HR / Elodie.3.1 model
ulyss, uly_root+'/data/m67.fits', cmp_peg, FILE='m67_phr', /CLEAN, /QUIET
uly_solut_tprint, 'm67_phr'
