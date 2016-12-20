; power is the component structure describing the power-low component
power = {name:'power', $
         descr:'Power-law component -- demonstration', $
         init_fun:'', $
         init_data:ptr_new(), $
         eval_fun:'tuto_power', $
         eval_data:ptr_new(/ALLO), $ ; will be filled by ULY_FIT_INIT
         para:ptr_new(/ALLO), $
         start:0d, $              ; this will be filled by ULY_FIT_INIT
         step:0d, $               ; this will be filled by ULY_FIT_INIT
         npix: 0l, $              ; this will be filled by ULY_FIT_INIT
         sampling:-1s, $          ; this will be filled by ULY_FIT_INIT
         mask:ptr_new(), $
         weight:0d, $
         e_weight:0d, $
         l_weight:0d, $
         lim_weig:[!VALUES.D_NaN,!VALUES.D_NaN] $
        }

; We have also to define the para structure, describing the parameter tau:
*power.para = {name:'tau', unit:'', guess:ptr_new(0D), step:1D-2, $
                       limits:[-10d,10d], limited:[1,1], fixed:0S, $
                       value:0D, error:0D, dispf:''}

; Since init_func is empty, ULY_FIT_INIT will fill the WCS information
; with the values corresponding to the spectrum to fit, and will
; pass the cmp structure as eval_data.
