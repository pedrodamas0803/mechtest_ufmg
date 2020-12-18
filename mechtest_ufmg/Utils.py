
import numpy as np 
import matplotlib.pyplot as plt 
import os
'''
Auxiliary functions module.
'''



def Hooke(strain, E = 150000, b = 0):
        '''
        Defines the linear equation for the linear portion of the stress-strain curve.

        Inputs:

        strain - the vector containing the strain values
        E - default = 150,000; initial guess for the modulus of elasticity in MPa.
        b - default = 0; intercept of the liear equation.
        '''

        return E * strain + b

def r_squared(x, y, model, modelParams):

        '''
        Calculates the R² statistics to evaluate the quality of the data regression.

        Inputs:

        x - the vector containing the x values of the independent variable of the function.
        y - the vector containing the dependent (measured) values as a function of x.
        model - the function that receives x as variable and calculates a y vector to model measured y.
        modelParams - the parameters of the model function passed as a list.

        Output:

        r_squared - the R² statistics. 
        '''

        res = y - model(x, *modelParams)
        ss_res = np.sum(res ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r_squared = 1 - (ss_res/ss_tot)
        
        return r_squared

def find_index(array, value):
    
        i = 0
        while array[i] < value:
                i = i + 1
        return i
  


def Hollomon(x, K = 300, n = 0.25):

        '''
        Defines Hollomon's equation.

        '''
        return K * x ** n

def Ludwik(x, sig_o = 300, K = 600, n = 0.24):

        return sig_o + K * x ** n

def Datsko(x, K = 300, x0 = 300, n = 0.20):

        return K * (x0 + x) ** n 
      

# def plot_mech(strain, stress, fig_label = 'Sample', stress_unit = 'MPa',
#                 show_plot = True, save = False, name = None, plot_type = 'ssc'):

#         if name.isnone():
#                 raise ValueError("Name must be provided")        
        
#         if plot_type not in ['ssc', 'young', 'yield', 'flow', 'true', 'model']:
#                 raise ValueError('Please provide the correct type of plot.')

#         elif plot_type == 'ssc':
#                 title = f'Engineering stress/strain curve'
#         elif plot_type == 'young':
#                 title = f'Apparent modulus of elasticity determination'
#         elif plot_type == 'yield':
#                 title = f'Yield strength determination'
#         elif plot_type == 'flow':
#                 title = f'Flow stress curve'
#         elif plot_type == 'true':
#                 title = f'True stress/strain curve'
#         elif plot_type == 'model':
#                 title = f'Plastic flow model'
       
#         plt.figure(figsize=(8, 4.5), facecolor = 'white')
#         plt.plot(strain, stress, 'b-', label = fig_label)
#         plt.xlabel('strain [mm/mm]')
#         plt.ylabel(f'stress [{stress_unit}]')
#         plt.xlim(0, 1.05 * max(strain))
#         plt.ylim(0, 1.05 * max(stress))
#         plt.title(title)
#         plt.legend(fontsize = 12, loc = 'lower right', frameon = False)

#         if save == True:
#             save_path = os.path.abspath(os.path.join('output', name))
#             plt.savefig(save_path, dpi = 300, bbox_inches = 'tight',transparent = False)
        
#         if show_plot == True:
#             plt.show()