from mechtest_ufmg.Properties import *
import pandas as pd 

data = pd.read_csv('data/astm1055.tsv', sep='\t', decimal=',')

strain = data['Strain']
stress = data['Stress']

plot_eng_SSC(strain, stress, fig_label = 'ASTM 1055', show_plot = False, save = True, name = 'eng_SSC_1055')

hooke = young_modulus(strain, stress, fig_label = 'ASTM 1055', show_plot = False, save = True, name = 'elasticity_1055')

sig_y = sigma_y(strain, stress, hooke[0], hooke[1], fig_label = 'ASTM 1055', show_plot =False, save = True, name = 'yield_strength_1055')

E = hooke[0]

uts = UTS(strain, stress, show = True)

uniform_elong(strain, stress, show = True)

aprox_resilience(E, sig_y, True)

resilience(strain, stress, sig_y)

aprox_toughness(strain, sig_y, uts, show = True)

toughness(strain, stress, show = True)

plot_flow_curve(strain, stress, sig_y, uts, fig_label = 'ASTM 1055', show_plot = False, save = True, name = 'flow_curve_1055')

plot_true_SSC(strain, stress, sig_y, uts, fig_label = 'ASTM 1055', show_plot = False, save = True, name = 'true_curve_1055')

flow_model(strain, stress, sig_y, uts, show_plot = False, save = True, name = 'Hollomon_1055')