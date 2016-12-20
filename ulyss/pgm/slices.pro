;+
; NAME:
;	SLICES
;
; PURPOSE:
;	Slice an image in horizontal/vertical direction
;
; USAGE:
;	slices, image, xx, yy
;
; INPUTS:
;	IMAGE:	Input image, which should be already displayed 
;       XX:     The x coordinates of the image
;       YY:     The y corrdinates of the image     
;
; PROCEDURE:
;       Plot interactively 2d slices of a 3d image in a seperate
;       window. Use the left mouse button to switch between row and
;       column profiles, while if the right mouse button is pressed
;       the program exits.
;
;       This routine is a modified version of the original IDL routine
;       profiles.pro. The differences are:
;       - The analysed image does not need to be scaled to the window
;         size. 
;       - It takes as input the x/y coordinates of the image
;       - It plots the profiles in DATA coordinates.
;
; EXAMPLE:
;	Create an chi2 map:
;       rangex = dindgen(10)*(10000.-3000.)/9.+3000
;       rangey = dindgen(10)/9.-1
;       map = uly_chimap(uly_root+'/data/VazMiles_z-0.40t07.94.fits', uly_ssp(),[[0,0],[0,1]], /QUIET, FILE_OUT='chimap.fits', XNODE=alog(rangex), YNODE=rangey)
;       Plot the resulting chi2 map:
;       uly_chimap_plot, map, /QUIET
;       Invoke slices:
;       slices, map.map, map.xnode, map.ynode
;
; HISTORY: This routine is a modification of the original IDL Library routine
;    Profiles, M. Koleva, 2012/11/08
;-
; CATEGORY:
;	ULY_UTIL
pro slices, image, xx, yy

  on_error,2                       ;Return to caller if an error occurs

  s = size(image)
  nx = s[1]                     ;Cols in image
  ny = s[2]                     ;Rows in image
  maxv = max(image, min=minv)             ;Get extrema
  orig_w = !d.window            ;Get main window def
  orig_xs = !x.s
  orig_ys = !y.s
  orig_zs = !z.s
  xtype = !X.TYPE
  ytype = !Y.TYPE
  ztype = !Z.TYPE
  xsize = 0.5*!D.X_VSIZE
  ysize = 0.5*!D.Y_VSIZE
  xstart = !D.X_SIZE + 0.5 * !D.X_SIZE
  ystart = !D.Y_SIZE - 0.15 * !D.Y_SIZE

  tvcrs, 1
  print, 'Left mouse button to toggle between rows and columns.'
  print, 'Right mouse button to Exit.'
  window, /FREE, xs=xsize, ys=ysize, xpos=xstart, ypos=ystart, title='Profiles' ;Make new window
  new_w = !d.window
  old_mode = -1                 ;Mode = 0 for rows, 1 for cols
  old_font = !p.font            ;Use hdw font
  !p.font = 0
  mode = 0
  !mouse.button = 0

  while (!mouse.button ne 4) do begin 
     wset, orig_w		;Image window
     !x.s = orig_xs
     !y.s = orig_ys 
     !z.s = orig_zs 
     !X.TYPE =  xtype
     !Y.TYPE = ytype
     !Z.TYPE = ztype
     cursor, a, b, 2, /DATA         ;Read position
     x =  where((abs(xx-a)) eq min(abs(xx-a)))
     y =  where((abs(yy-b)) eq min(abs(yy-b)))

     if !mouse.button eq 1 then begin
        mode = 1-mode           ;Toggle mode
        repeat begin
           cursor, a, b, 0, /DATA 
           x =  where((abs(xx-a)) eq min(abs(xx-a)))
           y =  where((abs(yy-b)) eq min(abs(yy-b)))
        endrep until !mouse.button eq 0
     endif

     device, window_state = ws
     wins = where(ws ne 0)
     if ((where(wins eq new_w))[0] ge 0) then begin
        wset, new_w
        endif else $             ;Graph window
           return

     if mode ne old_mode then begin
        old_mode = mode
        first = 1
        if mode then begin	;Columns?
           plot,[minv,maxv],minmax(yy),/nodata,title='Column Profile'
           vecy = yy
        end else begin
           plot,minmax(xx),[minv,maxv],/nodata,title='Row Profile'
           vecx = xx
        endelse
     endif

     if (x lt nx) and (y lt ny) and $
        (x ge 0) and (y ge 0) then begin           ;Draw it
        if first eq 0 then begin                   ;Erase?
           erase
        endif else first = 0
        value = string(xx[x],FORMAT='(i8)')+string(yy[y], Format='(+f8.2)')

        if mode then begin                ;Columns?
           vecx = image[x,*]              ;get column
           crossx = image[x,y]
           crossy = yy[y]
           title='Column Profile'
        endif else begin
           vecy = image[*,y]	;get row
           crossx = xx[x]
           crossy = image[x,y]
           title = 'Row Profile'
        endelse
        
;        print, !d.window
;        print, y, minmax(vecy), minmax(vecx), format='(i5, 2f+8.2, 2i10)'
        plot, vecx, vecy, TITLE=title            ;Graph
        xyouts, .1, 0, /NORM, value    ;Text of locn
        oplot, [crossx], [crossy], PSYM=1, SYMSIZE=3
     endif
  endwhile

;; Quit
  wset,orig_w
  !x.s = orig_xs
  !y.s = orig_ys 
  !z.s = orig_zs 
  tvcrs, 0                      ;Invisible cursor to the old window
  wdelete, new_w
  !p.font = old_font

end
; ------------- end ----------------------------
