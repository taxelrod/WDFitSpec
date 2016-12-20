m67 = uly_root+'/data/m67.fits'
model = uly_ssp_extr(m67, uly_root+'/models/Vaz_Miles.fits', [3920.,-0.1])
cmp = uly_star(model)
uly_lsf, m67, cmp, 400, 200, FILE='lsf_m67_vaz.txt', /QUIET
uly_lsf_smooth, 'lsf_m67_vaz.txt', 'lsfs_m67_vaz.txt'
cmp_vaz =  uly_ssp(MODEL=uly_root+'/models/Vaz_Miles.fits', LSF='lsfs_m67_vaz.txt')
ulyss, m67, cmp_vaz, KMOMENT=0, /CLEAN, FILE='m67_vaz_2', /PLOT
