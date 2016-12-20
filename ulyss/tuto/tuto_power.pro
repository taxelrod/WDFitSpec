function tuto_power, eval_data, pars
; pars is the array of parameters, actually only one value: the exponent
step = eval_data.step
npix = eval_data.npix

; simply return the array:
return, exp(findgen(npix)*step*pars[0])

end
