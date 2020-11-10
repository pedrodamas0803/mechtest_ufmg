import tensile_test as tt
import pandas as pd 

data = pd.read_csv('data/astm1016.tsv', sep='\t', decimal=',')

strain = data['Strain']
stress = data['Stress']

tt.plot_eng_SSC(strain, stress, fig_label = 'ASTM 1016', show_plot = False, save = True, name='eng_SSC_1016')

hooke = tt.young_modulus(strain, stress, fig_label = 'ASTM 1016', show_plot = False, save = True, name='elasticity_1016')

sig_y = tt.sigma_y(strain, stress, hooke[0], hooke[1], fig_label = 'ASTM 1016', show_plot = False, save = True, name='yield_strength_1016')

E = hooke[0]

uts = tt.UTS(strain, stress, show = True)

tt.uniform_elong(strain, stress, show=True)

tt.aprox_resilience(E, sig_y, True)

tt.resilience(strain, stress, sig_y)

tt.aprox_toughness(strain, sig_y, uts, show = True)

tt.toughness(strain, stress, show = True)

tt.plot_flow_curve(strain, stress, sig_y, uts, fig_label = 'ASTM 1016', show_plot = False, save = True, name='flow_curve_1016')

tt.plot_true_SSC(strain, stress, sig_y, uts, fig_label = 'ASTM 1016', show_plot = False, save = True, name='true_curve_1016')

tt.flow_model(strain, stress, sig_y, uts, show_plot = False, save = True, name='Hollomon_1016')