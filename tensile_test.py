import matplotlib.pyplot as plt 
import numpy as np 
import scipy as sp 
from scipy.optimize import curve_fit
from utils import Hooke


def plot_eng_SSC(strain, stress, save = False):
    
    if save == True:
        plt.figure(figsize=(8, 4.5), facecolor='white')
        plt.plot(strain, stress, 'b-')
        plt.xlabel('strain [mm/mm]')
        plt.ylabel('stress [MPa]')
        plt.xlim(0, 1.05* max(strain))
        plt.ylim(0, 1.05*max(stress))  
        plt.savefig('output/eng_SSC', dpi=300, bbox_inches='tight',transparent=False)
    else:
        plt.figure(figsize=(8, 4.5), facecolor='white')
        plt.plot(strain, stress, 'b-')
        plt.xlabel('strain [mm/mm]')
        plt.ylabel('stress [MPa]')
        plt.xlim(0, 1.05* max(strain))
        plt.ylim(0, 1.05*max(stress))
        plt.show()


    
def young_modulus(strain, stress, save = False):

    x = strain[strain < 0.002]
    y = stress[0:len(x)]    

    init_guess = [100000, 0]
    model = curve_fit(Hooke, x, y, p0 = init_guess)

    ans, cov = model
    E, b = ans    
    fit_curve = E * x + b
    E_gpa = round(E / 1000)

    plt.figure(figsize = (8,4.5))
    plt.plot(strain, stress, 'b-')
    plt.plot(x, fit_curve, 'r-', linewidth = 2 )
    plt.xlabel('strain [mm/mm]')
    plt.ylabel('stress [MPa]')
    plt.xlim(0, 1.05* max(strain))
    plt.ylim(0, 1.05*max(stress))
    plt.text(0,0, f'The elasticity modulus is {E_gpa} GPa', fontsize=12)  
    if save == False:     
        plt.show()
    else:
        plt.savefig('output/elasticity', dpi=300, bbox_inches='tight',transparent=False)

        
    return E, int(b)
    
