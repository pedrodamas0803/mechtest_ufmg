
import pandas as pd 
import matplotlib.pyplot as plt 
from mechtest_ufmg.Tensile_properties import *
from mechtest_ufmg.Utils import *
import numpy as np
import os

data = pd.read_csv('data/astm1016.tsv', sep='\t', decimal=',')

strain = data['Strain']
stress = data['Stress']

astm1055 = Tensile_test(strain, stress, specimen_name='ASTM 1016')


# astm1055.plot_young_modulus(save=True)

# astm1055.plot_yielding(save = True)

# astm1055.plot_true_curve(save = True)

# astm1055.plot_conventional_curve(save = True)

# astm1055.plot_flow_model(model = 'Datsko', save = True)

# astm1055.plot_flow_model(model = 'Ludwik', save = True)

# astm1055.plot_flow_model(save = True)

# astm1055.summary()
