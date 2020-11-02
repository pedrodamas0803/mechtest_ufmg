import tensile_test as tt
import pandas as pd 

data = pd.read_csv('data/am1.csv', sep = ';',
                    header = 1, usecols = ['(s)', '(mm)', '(kN)', '(mm/mm)', '(MPa)'],
                    skiprows = 4, encoding = 'latin_1')

strain = data['(mm/mm)']
stress = data['(MPa)']

tt.plot_eng_SSC(strain, stress, save = True)

hooke = tt.young_modulus(strain, stress, save = True)

tt.sigma_y(strain, stress, hooke[0], hooke[1], save = True)
# print(hooke)

uts = tt.UTS(strain, stress)
print(uts)