m67 = uly_root+'/data/m67.fits'
model = uly_ssp_extr(m67, uly_root+'/models/PHR_Elodie31.fits', [3920.,-0.1])
cmp = uly_star(model)
uly_lsf, m67, cmp, 400, 200, FILE='lsf_m67_phr.txt', /QUIET
uly_lsf_smooth, 'lsf_m67_phr.txt', 'lsfs_m67_phr.txt'
cmp_phr =  uly_ssp(LSF='lsfs_m67_phr.txt')
ulyss, m67, cmp_phr, KMOMENT=0, /CLEAN, FILE='m67_phr_2', /PLOT
