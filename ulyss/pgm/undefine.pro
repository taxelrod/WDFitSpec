;+
; NAME:
;                       UNDEFINE
;
; PURPOSE:
;                       Destroy a variable
;
; USAGE:
;                       undefine, varname
;
; ARGUMENT:
;   varname : Named variable to delete
;
; DESCRIPTION:
;   Deletes a variable so that after n_elements(varname) will be 0.
;   Note that if the variable is a pointer or structure containing
;   pointers, this routine does not free the memory.
;
; HISTORY:
;   Philippe Prugniel, 2008, created
;-
; CATEGORY:    ULY_UTIL
;------------------------------------------------------------------------------
PRO undefine, varname
if n_elements(varname) eq 0 then return
tempvar = size(temporary(varname))
END
