; Substitutes for intrinsic IDL routines missing in GDL

pro heap_free, var, OBJ=obj, PTR=ptr, VERBOSE=verbose

if size(var, /TYPE) eq 10 then begin
    for i = 0, n_elements(var)-1 do begin
        if ptr_valid(var[i]) then begin
            heap_free, *var[i] 
            if keyword_set(verbose) then help, *var[i]
            ptr_free, var[i]
        endif
    endfor
    return
endif
if size(var, /TYPE) ne 8 then return

if keyword_set(obj) then message, 'Object heap variables are not freed', /INFO

for i = 0, n_elements(var)-1 do begin
    tags = tag_names(var[i])
    for j=0,n_elements(tags)-1 do begin ; scan all the tags in the structure
        ok = Execute('heap_free, OBJ=obj, PTR=ptr, VERBOSE=verbose, var[i].' + tags[j])
    endfor
endfor

end

pro file_delete, f1, f2, f3, f4, f5, f6, f7
; minimum implementation (does not support any keyword!)

if n_elements(f1) eq 0 then return
if n_elements(f1) gt 1 then begin
    spawn, 'rm ' + strjoin(f1, ' ')
    return
endif

cmd = 'rm '+ f1
if n_elements(f2) eq 1 then cmd += ' ' + f2
if n_elements(f3) eq 1 then cmd += ' ' + f3
if n_elements(f4) eq 1 then cmd += ' ' + f4
if n_elements(f5) eq 1 then cmd += ' ' + f5
if n_elements(f6) eq 1 then cmd += ' ' + f6
if n_elements(f7) eq 1 then cmd += ' ' + f7
spawn, 'rm ' + strjoin(f1, ' ')

end
;------------------------------------------------------------------------------
