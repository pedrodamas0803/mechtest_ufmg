import matplotlib.pyplot as plt 
import numpy as np 
import scipy as sp 
from scipy.optimize import curve_fit
from utils import Hooke


def plot_eng_SSC(strain, stress, save = False):
    plt.figure(figsize=(8, 4.5), facecolor='white')
    plt.plot(strain, stress, 'b-')
    plt.xlabel('strain [mm/mm]')
    plt.ylabel('stress [MPa]')
    plt.xlim(0, 1.05* max(strain))
    plt.ylim(0, 1.05*max(stress))  
    
    if save != False:       
        plt.savefig('output/eng_SSC', dpi=300, bbox_inches='tight',transparent=False)
    else:
        plt.show()


    
def young_modulus(strain, stress, save = False):
    
# taking only the elastic portion of a curve
    x = strain[strain < 0.002]
    y = stress[0:len(x)]    
# performing the linear regression on the elastic part of the data
    init_guess = [100000, 0]
    model = curve_fit(Hooke, x, y, p0 = init_guess)
    ans, cov = model
    E, b = ans    
    fit_curve = E * x + b
    E_gpa = round(E / 1000)
    
# calculating the R_squared statistic
    res = y - Hooke(x, *ans)
    ss_res = np.sum(res ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res/ss_tot)   

# plotting the figure and showing or saving the figure
    plt.figure(figsize = (8,4.5))
    plt.plot(strain, stress, 'b-')
    plt.plot(x, fit_curve, 'r-', linewidth = 2 )
    plt.xlabel('strain [mm/mm]')
    plt.ylabel('stress [MPa]')
    plt.xlim(0, 1.05* max(strain))
    plt.ylim(0, 1.05*max(stress))
    plt.text(0.1*max(strain), 0.1*max(stress), f'The elasticity modulus is {E_gpa} GPa, R² = {round(r_squared, 4)}', fontsize=12)  

    if save != False:     
        plt.savefig('output/elasticity', dpi=300, bbox_inches='tight',transparent=False)
    else:
        plt.show()
    
    return E, int(b), r_squared
    
