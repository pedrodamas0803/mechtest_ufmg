from mechtest_ufmg.Utils import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps, trapz
from scipy.optimize import curve_fit
import os
import csv
# from mechtest_ufmg.Utils import *
# Conversion factor from psi to MPa
psi_to_mpa = 0.00689476

class Tensile_test:
    '''
    Tensile_test

    This function will receive you data and create an object with your data entries.
    UNITS: The expected measurements units of the data are mm/mm² for length/area sizes, and MPa/Newton for stress/load.
            All math will be performed according to international standards units.

    Inputs:
    x - an array or series containing the strain, or displacement of your sample test.
    y - an array or series containing the stress, of force readings of your sample test.
    is_strain - default = True; it must be set to False if displacement data is provided.
    is_stress - default = True; it must be set to False if force data is provided.
    stress_unit - default = 'MPa'; it must be set to 'psi' if the data is psi.
    specimen_name - default = 'Sample'; the name that will appear in the plots.
    specimen_length - the specimen length in mm, it must be provided if is_strain != True.
    crossec_area - the cross section of the specimen in mm², it must be provided if is_stress != True.
    '''
    

    def __init__(self, x, y, is_strain = True, is_stress = True, stress_unit = 'MPa', specimen_name = 'Sample',
                specimen_length = None, crossec_area = None):

        if is_strain not in [True, False]:
            raise ValueError('is_strain must be True or False.')

        if is_stress not in [True, False]:
            raise ValueError('is_stress must be True or False.')

        if stress_unit not in ['MPa', 'psi']:
            raise ValueError("stress unit must be 'MPa' or 'psi'." )
        elif stress_unit == 'psi' and is_stress == True:
            y = y * psi_to_mpa
            
        if is_strain != True and specimen_length.isnumeric() != True:
            raise ValueError('Specimen length must be provided if is_strain == False')

        if is_stress != True and crossec_area.isnumeric() != True:
            raise ValueError('Cross-section area must be provided if is_stress != True')


        self.name = specimen_name
        self.crossec_area = crossec_area
        self.specimen_length = specimen_length
        self.stress_unit = stress_unit

        if is_strain == True:
            self.strain = x
        else:
            self.strain = x/specimen_length

        if is_stress == True:
           
            self.stress = y
        else:
            self.strain = x/crossec_area

        create_output_folder()


    @property
    def young_modulus(self):
        # taking only the elastic portion of a curve
        x = self.strain[self.strain < 0.002]
        y = self.stress[0:len(x)]

        # performing the linear regression on the elastic part of the data
        init_guess = [100000, 0]
        model = curve_fit(Hooke, x, y, p0 = init_guess)
        ans, *_ = model
        E_mpa, intercept = ans
        fit_curve = E_mpa * x + intercept
        E_gpa = round(E_mpa / 1000)

        # calculating the R_squared statistic
        r2 = r_squared(x, y, Hooke, ans) 
        
        return E_mpa, E_gpa, intercept, r2
        

    @property
    def yield_strength(self):

        '''
        Calculates the yielding stress as the stress related to a permanent deformation of 0.002.

        Outputs:
        sigma_yield - the yielding stress in the same unit of stress input.

        '''
        # taking part of the data
        x = self.strain[self.strain < 0.05]
        k = x - 0.002
        z = Hooke(k, self.young_modulus[0])
        # finding the index of sig_y
        i = 0
        while self.stress[i] > z[i]:
            i = i + 1

        sig_y = self.stress[i]  
        
        return sig_y

    @property
    def UTS(self):

        '''
        Calculates the ultimate tensile stress taking the higher stress from the stress vector input.

        Output:
        uts - the ultimate tensile strength in the same unit of the stress input.

        '''

        return max(self.stress)

    @property
    def uniform_elongation(self):

        '''
        Calculates the uniforme elongation from the input data taking the strain point related to the uts value.

        Output:
        u_elong - the uniform elongation of the sample.
        '''
        # transform data to numpy arrays
        eps = self.strain.to_numpy()
        sig = self.stress.to_numpy()

    # finds the index of the maximum stress
        index = np.argmax(sig)
        u_elong = eps[index]

        return u_elong

    @property
    def approx_resilience(self):

        '''
        Calculates the resilience by the approximation formula U_r = sigma_yield^2/E.

        Output:

        U_r - the resilience calculated by the presented formula.

        '''

        U_r = (self.yield_strength ** 2)/(2 * self.young_modulus[0])

        return U_r

    @property
    def resilience(self):

        '''
        Calculates the resilience by numerical integration using the trapezoidal method up to the yield stregth point.

        Output:

        uts - the ultimate tensile strength in the same unit of the stress input.

        '''

        i = find_index(self.stress, self.yield_strength)
        eps = self.strain[0:i]
        sig = self.stress[0:i]
        U_r = trapz(sig, eps)

        return U_r

    @property
    def approx_toughness(self):

        '''
        Calculates an approximation of the toughness using the formula U_t = eps_f * ((sig_y + uts)/2).

        Output:

        U_t - the approximate toughness of the material
        '''

        eps_f = max(self.strain)

        U_t = eps_f * ((self.yield_strength + self.UTS) / 2)
        

        return U_t 

    @property
    def toughness(self):

        '''
        Calculates the toughness of the material performing the numerical integration of the stress-strain curve.

        Output:

        mat_toughness - the material toughness calculated by numerical integration, in MJ/m³.
        '''

        mat_toughness = trapz(self.stress, self.strain)
        
        return mat_toughness

    def uniform_deformation(self):
        # TODO: deal with discontinuous yielding; how can I determine where it ends?

        '''
        Plots the flow stress curve from the material stress-strain curve with values between the yield strength and the ultimate tensile strength.

        Outputs:    

        Figure that can be displayed in matplotlib interface and/or saved to the output folder.
        '''
        eps = self.strain.to_numpy()
        sig = self.stress.to_numpy()

        yield_index = find_index(sig, self.yield_strength)
        uts_index = find_index(sig, self.UTS)

        eps_c = eps[yield_index:uts_index]
        sig_c = sig[yield_index:uts_index]

        return eps_c, sig_c


    def real_values(self):
        
        '''
        Plots the true stress-strain conversion of the input data.

        Output:

        eps_t, sig_t - the true strain/stress values calculated between the yield strength and ultimate tensil strength.
        '''
        strain, stress = self.uniform_deformation()
             
        eps_t = np.log(1 + strain)
        sig_t = stress * (1 + strain)

        return eps_t, sig_t
       
    def flow_model(self, model = 'Hollomon'):

        '''
        Calculates the regression coefficients for a model of plasticity, i.e. Hollomon's equation.

        Input:

        model - default = 'Hollomon' - the model to describe the plastic flow of the material.
                                    Available: Hollomon, Ludwik and Datsko.
        Output:

        The coefficients of each model and the R² statistics.        
        '''

        if model not in ['Hollomon', 'Ludwik', 'Datsko']:
            raise ValueError('The model must be Hollomon, Ludwik, or Datsko')

        x, y = self.real_values()

        if model == 'Hollomon':

            init_guess = [300, 0.25]
            
            model_fit = curve_fit(Hollomon, x, y, p0 = init_guess)

            ans, *_ = model_fit
            Kexp, nexp = ans

            sig_h = Hollomon(x, K = Kexp, n = nexp)

            R2 = r_squared(x, y, Hollomon, ans)

            return nexp, Kexp, R2

        elif model == 'Ludwik':

            init_guess = [600, 600, 0.25]
            
            model_fit = curve_fit(Ludwik, x, y, p0 = init_guess)

            ans, *_ = model_fit
            sig_0, K, n = ans

            sig_h = Ludwik(x, sig_o = sig_0, K = K, n = n)

            R2 = r_squared(x, y, Ludwik, ans)

            return sig_0, K, n, R2

        elif model == 'Datsko':

            init_guess = [300, 300, 0.25]
            
            model_fit = curve_fit(Datsko, x, y, p0 = init_guess)

            ans, *_ = model_fit
            K, sig_0, n = ans

            sig_h = Datsko(x, K = K, x0 = sig_0, n = n)

            R2 = r_squared(x, y, Datsko, ans)

            return sig_0, K, n, R2


    def plot_conventional_curve(self, save = False):

        '''
        Plots and displays the engineering stress/strain curve for the provided data.

        Input:

        save - default = False, if True, the figure will be saved to the output folder.

        Output:

        The engineering stress/strain curve. 
        '''

        plt.figure(figsize=(8, 4.5))
        plt.plot(self.strain, self.stress, color = 'blue', label = self.name)
        plt.xlabel('strain [mm/mm]')
        plt.ylabel(f'stress[{self.stress_unit}]')
        plt.xlim(0, 1.05 * max(self.strain))
        plt.ylim(0, 1.05 * max(self.stress))
        plt.title('Engineering stress/strain curve')
        plt.legend(fontsize=12, loc = 'lower right', frameon = False)

        if save == True:

            save_path = os.path.abspath(os.path.join('output', self.name + '_engineering'))
            plt.savefig(save_path, dpi = 300, bbox_inches = 'tight', transparent = False)

        plt.show()

    def plot_young_modulus(self, save = False):

        '''
        Plots and displays the engineering stress/strain curve for the provided data highlighting
        the data used to calculate the apparent young modulus.

        Input:

        save - default = False, if True, the figure will be saved to the output folder.

        Output:

        Figure. 
        '''
        

        E_mpa, E_gpa, intercept, r2 = self.young_modulus

        x = self.strain[self.strain < 0.002]

        plt.figure(figsize=(8, 4.5), edgecolor='white', facecolor='white')
        plt.plot(self.strain, self.stress, color = 'blue', label = self.name)
        plt.plot(x, Hooke(x, E=E_mpa, b= intercept), color = 'red', label = 'Fitted E',
                linewidth = 2)
        plt.xlabel('strain [mm/mm]')
        plt.ylabel(f'stress [{self.stress_unit}]')
        plt.xlim(0, 1.05 * max(self.strain))
        plt.ylim(0, 1.05 * max(self.stress))
        plt.title(f'Modulus of elasticity determination')
        plt.legend(fontsize = 12, loc = 'lower right', frameon = False)
        plt.text(0.1 * max(self.strain), 0.1 * max(self.stress), f'The elasticity modulus is {E_gpa} GPa, R² = {round(r2, 4)}', fontsize = 12)

        if save:
            save_path = os.path.abspath(os.path.join('output', self.name + '_elasticity'))
            plt.savefig(save_path, dpi = 300, bbox_inches = 'tight',transparent = False)
        
        plt.show()

    def plot_yielding(self, save = False):

        '''
        Plots and displays the engineering stress/strain curve for the provided data highlighting the
        calculation of the yield stress.

        Input:

        save - default = False, if True, the figure will be saved to the output folder.

        Output:

        Figure. 
        '''

        E_mpa, E_gpa, intercept, r2 = self.young_modulus

        x = self.strain[self.strain < 0.05]
        k = x - 0.002
        z = Hooke(k, E_mpa)

        plt.figure(figsize = (8,4.5))
        plt.plot(self.strain, self.stress, 'b-', label = self.name)
        plt.plot(x, Hooke(x, E_mpa, intercept), 'r--', linewidth = 1, label = 'Hooke\'s law' )
        plt.plot(x, Hooke(k, E_mpa), 'r:', linewidth = 1 )
        plt.xlabel('strain [mm/mm]')
        plt.ylabel(f'stress [{self.stress_unit}]')
        plt.xlim(0, 1.05 * max(self.strain))
        plt.ylim(0, 1.05 * max(self.stress))
        plt.title(f'Yield strength determination')
        plt.legend(fontsize = 12, loc = 'lower right', frameon = False)
        plt.text(0.1 * max(self.strain), 0.1 * max(self.stress), 
                f'The yield strength is {round(self.yield_strength)} MPa.', fontsize = 12)

        if save == True:

            save_path = os.path.abspath(os.path.join('output', self.name + '_yielding'))
            plt.savefig(save_path, dpi = 300, bbox_inches = 'tight',transparent = False)

        plt.show()

    def plot_flow_curve(self, save = False):

        '''
        Plots and displays the flow curve for the provided data.

        Input:

        save - default = False, if True, the figure will be saved to the output folder.

        Output:

        Figure.
        '''

        yield_index = find_index(self.stress, self.yield_strength)

        plt.figure(figsize = (8,4.5))
        plt.plot(self.strain[yield_index: ], self.stress[yield_index: ],
                 color = 'blue', label = self.name)
        plt.xlabel('strain [mm/mm]')
        plt.ylabel(f'stress [{self.stress_unit}]')
        plt.xlim(0, 1.05 * max(self.strain))
        plt.ylim(0, 1.05 * max(self.stress))
        plt.title(f'Flow stress curve')
        plt.legend(fontsize = 12, loc = 'lower right', frameon = False)

        if save == True:

            save_path = os.path.abspath(os.path.join('output', self.name + '_flow'))
            plt.savefig(save_path, dpi = 300, bbox_inches = 'tight',transparent = False)

        plt.show()

    def plot_true_curve(self, save = False):

        '''
        Plots and displays the true stress/strain curve for the provided data.

        Input:

        save - default = False, if True, the figure will be saved to the output folder.

        Output:

        Figure.
        '''

        x, y = self.real_values()

        plt.figure(figsize=(8, 4.5))
        plt.plot(x, y, color = 'blue', label = self.name)
        plt.xlabel('true strain [mm/mm]')
        plt.ylabel(f'true stress [{self.stress_unit}]')
        plt.xlim(0, 1.05 * max(x))
        plt.ylim(0, 1.05 * max(y))
        plt.title(f'True stress/strain curve')
        plt.legend(fontsize = 12, loc = 'lower right', frameon = False)

        if save == True:

            save_path = os.path.abspath(os.path.join('output', self.name + '_true'))
            plt.savefig(save_path, dpi = 300, bbox_inches = 'tight', transparent = False)

        plt.show()

    def plot_flow_model(self, model = 'Hollomon', save = False):

        '''
        Plots and displays the flow stress/strain curve for the provided data and the model
        derived from it by linear regression.

        Input:

        save - default = False, if True, the figure will be saved to the output folder.

        Output:

        The engineering stress/strain curve. 
        '''
        

        if model not in ['Hollomon', 'Ludwik', 'Datsko']:
            raise ValueError('Please provide a valid model.')

        x, y = self.real_values()

        plt.figure(figsize = (8, 4.5))
        plt.plot(x, y, color = 'blue', label = self.name)
        
        if model == 'Hollomon':

           nexp, Kexp, R2 = self.flow_model(model = 'Hollomon')

           fit_y = Hollomon(x, K = Kexp, n = nexp)

           plt.plot(x, fit_y, 'r--', label = r'$\sigma = K \epsilon ^ n$')
           plt.plot([], [], ' ', label = f'K = {round(Kexp)}\nn = {round(nexp, 3)}\nR² = {round(R2, 3)}')

        if model == 'Ludwik':

           sig_0, Kexp, nexp, R2 = self.flow_model(model = 'Ludwik')

           fit_y = Ludwik(x, sig_o = sig_0, K = Kexp, n = nexp)

           plt.plot(x, fit_y, 'r--', label = r'$\sigma = \sigma_0 + K \epsilon ^ n$')
           plt.plot([], [], ' ', 
                    label = f'K = {round(Kexp)}\nn = {round(nexp,3)}\n' + r'$\sigma_0 =$' + f'{round(sig_0)}\nR² = {round(R2, 3)}')
        
        if model == 'Datsko':

           sig_0, Kexp, nexp, R2 = self.flow_model(model = 'Datsko')

           fit_y = Datsko(x, x0 = sig_0, K = Kexp, n = nexp)

           plt.plot(x, fit_y, 'r--', label = r'$\sigma = K (\epsilon_0 + \epsilon) ^ n$')
           plt.plot([], [], ' ', 
                    label = f'K = {round(Kexp)}\nn = {round(nexp,3)}\n' + r'$\epsilon_0 =$' + f'{round(sig_0)}\nR² = {round(R2, 3)}')



        plt.xlabel('strain [mm/mm]')
        plt.ylabel(f'stress [{self.stress_unit}]')
        plt.xlim(0, 1.05 * max(x))
        plt.ylim(0, 1.05 * max(y))
        plt.title(f'{model} flow model fitting')
        plt.legend(fontsize = 12, loc = 'lower right', frameon = False)
        
        if save == True:

            save_path = os.path.abspath(os.path.join('output', self.name + f'{model} fit'))
            plt.savefig(save_path, dpi = 300, bbox_inches = 'tight', transparent = False)
        
        plt.show()


    def summary(self):

        '''
        Saves a csv file with the properties calculated from the input data. 

        '''

        save_path = os.path.join('output', f'{self.name}_summary.csv')

        summ = [('Property [unit]', 'Value'),
                ('Elasticity mod [GPa]', round(self.young_modulus[1])), 
                ('Yield strength [MPa]', round(self.yield_strength)),
                ('Ultimate tens strength [MPa]', round(self.UTS)),
                ('Uniform elongation [mm/mm]', round(self.uniform_elongation, 4)),
                ('Approx. resilience [MJ/m³]', round(self.approx_resilience, 3)),
                ('Resilience [MJ/m³]', round(self.resilience, 3)),
                ('Approx. toughness [MJ/m³]', round(self.toughness, 3)),
                ('Toughness [MJ/m³]', round(self.toughness, 3))]

        if os.path.exists(save_path):

            with open(save_path, 'w') as csvfile:
                
                writer = csv.writer(csvfile)
                writer.writerows(summ)         

        else:

            with open(save_path, 'x') as csvfile:
                
                writer = csv.writer(csvfile)
                writer.writerows(summ)

