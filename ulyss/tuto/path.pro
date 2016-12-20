; alter path and uly_root for the tutorial
; uly_root is altered in order to hide the full path of the file
; (useless and inelegant)
save_path = !path
save_root = uly_root
!path += ':' + uly_root + '/tuto'
uly_root = '.'
