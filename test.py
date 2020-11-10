import tensile_test as tt
import pandas as pd 

data = pd.read_csv('data/am1.csv', sep = ';',
                    header = 1, usecols = ['(s)', '(mm)', '(kN)', '(mm/mm)', '(MPa)'],
                    skiprows = 4, encoding = 'latin_1')

strain = data['(mm/mm)']
stress = data['(MPa)']

tt.plot_eng_SSC(strain, stress, fig_label = 'ASTM 1016', show_plot = False, save = True)

hooke = tt.young_modulus(strain, stress, fig_label = 'ASTM 1016', show_plot = False, save = True)

E = hooke[0]

sig_y = tt.sigma_y(strain, stress, hooke[0], hooke[1], fig_label = 'ASTM 1016', show_plot = False, save = True)


uts = tt.UTS(strain, stress, show = True)

tt.uniform_elong(strain, stress, show=True)

tt.aprox_resilience(E, sig_y, True)

tt.resilience(strain, stress, sig_y)

tt.aprox_toughness(strain, sig_y, uts, show = True)

tt.toughness(strain, stress, show = True)

tt.plot_flow_curve(strain, stress, sig_y, uts, fig_label = 'ASTM 1016', show_plot = False, save = True)

tt.plot_true_SSC(strain, stress, sig_y, uts, fig_label = 'ASTM 1016', show_plot = False, save = True)

tt.flow_model(strain, stress, sig_y, uts, show_plot = False, save = True)

# print(eps_u)