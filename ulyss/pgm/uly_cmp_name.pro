;----------------------------------------------------------------------------
;+
; NAME:           
;		  ULY_CMP_NAME
;
; PURPOSE:        
;		  Create a unique name for a fit component
;              
; USAGE:
;                 name_unique = uly_cmp_name([name])
;
; ARGUMENTS:
;     [NAME]:  
;		  Proposed name
;
; DESCRIPTION:
;     ULY_CMP_NAME returns a unique name for a fit component.
;     If the proposed name is already used, a new
;     name is searched by appending a number. If no name is passed
;     in input, the routine creates one like 'cmp<n>', where <n> is
;     a number (incremented after each call to this routine).
;
;     The name is also trimmed from embedded space characters, and
;     the characters '_' are replaced with '-'.
;
;     ULY_CMP_NAME is used by ULY_SSP, ULY_STAR, ULY_TGM when a fit
;     component is being defined.
;
; HISTORY:
;     CREATION  2008/08/10 Philippe Prugniel
;-
; CATEGORY:    ULY
;
function uly_cmp_name, name
common uly_cmp_namelist, cmp_names

nnames = n_elements(cmp_names) + 1

if n_elements(name) gt 0 then begin
    if strlen(strtrim(name)) gt 0 then name_out = strcompress(name,/REM) $
    else name_out = 'cmp'+strtrim(string(nnames),2)
endif else name_out = 'cmp'+strtrim(string(nnames),2)

if n_elements(cmp_names) gt 0 then begin
    if max(cmp_names eq name_out) eq 1 then begin ; name already exists
        name_try = name_out + strtrim(string(nnames),2)
        k = nnames+1
        while max(cmp_names eq name_try) eq 1 do begin
            name_try = name_out + strtrim(string(k),2)
            k++
        endwhile
        name_out = name_try
    endif 
    cmp_names = [cmp_names, name_out]
endif else cmp_names = name_out

p = strpos(name_out, '_')
while p gt -1 do begin
    strput, name_out, '-', p
    p = strpos(name_out, '_')
endwhile

return, name_out
end
