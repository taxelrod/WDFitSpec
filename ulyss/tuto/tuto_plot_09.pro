black_and_white = uly_plot_init()
black_and_white.linecolor=FSC_Color('Black')
black_and_white.modelcolor='Grey'
black_and_white.polycolor='Black'
black_and_white.polystyle=2
black_and_white.showsigma=0

uly_solut_splot, solution, PLOT_VAR=black_and_white

save, black_and_white, filename='bw_settings.sav'
