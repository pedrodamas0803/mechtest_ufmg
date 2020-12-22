import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps, trapz
from scipy.optimize import curve_fit
import mechtest_ufmg.Utils 
import os
from mechtest_ufmg.Utils import *

class Tensile_test:
    '''
    Tensile_test

    This function will receive you data and create an object with your data entries.
    UNITS: The expected measurements units of the data are mm/mm² for length/area sizes, and MPa/Newton for stress/force.
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
    # Conversion factor from psi to MPa
    psi_to_mpa = 0.00689476

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
        self.young_modulus = 0
        # self.intercept = 0
        # self.young_gpa = 0
        # self.young_r2 = 0

        if is_strain == True:
            self.strain = x
        else:
            self.strain = x/specimen_length

        if is_stress == True:
           
            self.stress = y
        else:
            self.strain = x/crossec_area
    
    @classmethod
    def young_modulus(cls):
        # taking only the elastic portion of a curve
        x = cls.strain[cls.strain < 0.002]
        y = cls.stress[0:len(x)]

        # performing the linear regression on the elastic part of the data
        init_guess = [100000, 0]
        model = curve_fit(Hooke, x, y, p0 = init_guess)
        ans, cov = model
        E, intercept = ans
        fit_curve = E * x + intercept
        E_gpa = round(E / 1000)

        # calculating the R_squared statistic
        r2 = r_squared(x, y, Hooke, ans) 
        
        cls.young_modulus = E

    


    # def ultimate_tens_stren(strain, stress, show = True):

    #     '''
    #     Calculates the ultimate tensile stress taking the higher stress from the stress vector input.

    #     Inputs:
    #     strain - vector containing the strain data in admensional units, if the data is provided in %, divide it by 100; e.g.: [mm/mm].
    #     stress - vector containing the stress data relative to the strain vector.
    #     show - default = True; if true, prints the value of UTS in the command line.

    #     Output:
    #     uts - the ultimate tensile strength in the same unit of the stress input.

    #     '''
            
    #     uts = max(stress)
    #     if show == True:
    #         print(f'The ultimate tensile strength is {round(uts)} MPa. ')

    #     return uts

    # def uniform_elongation(strain, stress, show = True):

    #     '''
    #     Calculates the uniforme elongation from the input data taking the strain point related to the uts value.

    #     Inputs:
    #     strain - vector containing the strain data in admensional units, if the data is provided in %, divide it by 100; e.g.: [mm/mm].
    #     stress - vector containing the stress data relative to the strain vector.
    #     show - default = True; if true, prints the value of UTS in the command line.

    #     Output:
    #     u_elong - the uniform elongation of the sample.
    #     '''
    #     # transform data to numpy arrays
    #     eps = strain.to_numpy()
    #     sig = stress.to_numpy()

    # # finds the index of the maximum stress
    #     index = np.argmax(sig)
    #     u_elong = eps[index]

    # # prints depending on the user's choice
    #     if show == True:
    #         print(f'The uniform elongation is {round(u_elong, 4)}.')

    #     return u_elong

    # def aprox_resilience(E, sig_y, show = True):

    #     '''
    #     Calculates the resilience by the approximation formula U_r = sigma_yield^2/E.

    #     Inputs:
    #     E - the apparent elasticity modulus, in MPa.
    #     sig_y - the yield strength, in MPa
    #     show - default = True; if True prints a string containing the calculated value.

    #     Output:

    #     U_r - the resilience calculated by the presented formula.

    #     '''

    #     U_r = (sig_y ** 2)/(2*E)

    #     if show == True:
    #         print(f'The resilience calculated by the formula yield²/2E is {round(U_r,2)} MJ/m³. ')

    #     return U_r

    # def resilience(strain, stress, sig_y, dx = 1.0, show = True):

    #     '''
    #     Calculates the resilience by numerical integration using the trapezoidal method up to the yield stregth point.

    #     Inputs:
    #     strain - vector containing the strain data in admensional units, if the data is provided in %, divide it by 100; e.g.: [mm/mm].
    #     stress - vector containing the stress data relative to the strain vector.
    #     sig_y - the yield strength in MPa.
    #     dx - default = 1.0; the step size for the integration.
    #     show - default = True; if true, prints the value of UTS in the command line.

    #     Output:

    #     uts - the ultimate tensile strength in the same unit of the stress input.

    #     '''
        
    #     i = find_index(stress, sig_y)
    #     eps = strain[0:i]
    #     sig = stress[0:i]
    #     U_r = trapz(eps, sig, dx)

    #     if show == True:
    #         print(f'The resilience calculated by trapezoidal integration with dx = {dx} is {round(U_r,2)} MJ/m³. ')

    #     return U_r


    # def aprox_toughness(strain, sig_y, uts, show = True):

    #     '''
    #     Calculates an approximation of the toughness using the formula U_t = eps_f * ((sig_y + uts)/2).

    #     Inputs:

    #     strain - the vector containing the strain values from the test.
    #     sig_y - the yield strength, in MPa
    #     uts - the ultimate tensile strength, in MPa
    #     show - default = True; if True prints the string containing the material toughness. 

    #     Output:

    #     U_t - the approximate toughness of the material
    #     '''

    #     eps_f = max(strain)

    #     U_t = eps_f * ((sig_y + uts) / 2)
        
    #     if show == True:
    #         print(f'The material toughness is approximately {round(U_t)} MJ/m³.')

    #     return U_t 


    # def toughness(strain, stress, show = True):

    #     '''
    #     Calculates the toughness of the material performing the numerical integration of the stress-strain curve.

    #     Inputs:

    #     strain - the vector containing the strain values obtained from the tests.
    #     stress - the vector containing the stress values that refer to the strain vector.
    #     show - default = True; if True prints a string containing the toughness in the command line.

    #     Output:

    #     mat_toughness - the material toughness calculated by numerical integration, in MJ/m³.
    #     '''

    #     mat_toughness = trapz(stress, strain)

    #     if show == True:
    #         print(f'The material toughness computed by Trapezoidal method is {round(mat_toughness)} MJ/m³. ')
        
    #     return mat_toughness

 

    
    # def flow_model(strain, stress, sig_y, uts, func = 'Hollomon', show = True, show_plot = True, save = False, name = 'flow_model'):

        # TODO: implement functions other than Hollomon.
        '''
        Calculates the regression coefficients for a model of plasticity, i.e. Hollomon's equation.

        Inputs:

        strain - the vector containing the strain values obtained from the tests.
        stress - the vector containing the stress values that refer to the strain vector.
        sig_y - the yield strength, in MPa.
        uts - the ultimate tensile strength, in MPa.
        fig_label - default = Sample; the string that identify the sample name in the output figure.
        show_plot - default = True; if True, the plot is shown in a matplotlib interface.
        save - default = False, if True, saves the figure in the folder output;
        name - default = flow_model; the name of the file saved in the output folder.

        Output:

        shex - strain hardening coefficient, n.
        Koeff - resistance modulus K, in MPa.
        
        Figure that can be displayed in matplotlib interface and/or saved to the output folder.
        # '''

        # eps_c, sig_c = uniform_plast(strain, stress, sig_y, uts)

        # eps_t, sig_t = true_values(eps_c, sig_c)

        # eps_tl = np.log(eps_t)
        # sig_tl = np.log(sig_t)

        # init_guess = [300, 0.25]

        # ans, cov = curve_fit(log_Hollomon, eps_tl, sig_tl, p0 = init_guess)

        # Koeff = ans[0]
        # shex = ans[1]

        # sig_h = log_Hollomon(eps_tl, K = Koeff, n = shex)

        # R2 = r_squared(eps_tl, sig_tl, log_Hollomon, ans)

        # plt.figure(figsize = (8,4.5))
        # plt.plot(eps_tl, sig_tl, 'b-', label = 'Stress-Strain log values')
        # plt.plot(eps_tl, sig_h, 'r:', label = f'Linear Hollomon, K = {round(Koeff)} MPa, n = {round(shex, 2)}, R²={round(R2, 4)}')
        # plt.xlabel('log(strain) [mm/mm]')
        # plt.ylabel('log(stress) [MPa]')
        # plt.title(f'{func} fitted to the data')
        # plt.legend(fontsize = 12, loc = 'lower right', frameon = False)        

        # if show_plot == True:

        #     plt.show()

        # if save == True:
            
        #     save_path = os.path.abspath(os.path.join('output', name))
        #     plt.savefig(save_path, dpi = 300, bbox_inches = 'tight',transparent = False)
            
        # if show == True:

        #     print(f'The resistance modulus is {round(Koeff)} MPa and the strain-hardening exponent is {round(shex, 2)}.')

        # return shex, Koeff


if __name__ == "__main__":
    
    import pandas as pd 

    data = pd.read_csv('tests/data/astm1055.tsv', sep='\t', decimal=',')

    strain = data['Strain']
    stress = data['Stress']

    astm1055 = Tensile_test(strain, stress, specimen_name='ASTM 1055')

    # astm1055.plot_eng_SSC(show_plot=True)
    astm1055.young_modulus()
    print(astm1055.young_modulus)
    # print(astm1055.young_gpa)
        