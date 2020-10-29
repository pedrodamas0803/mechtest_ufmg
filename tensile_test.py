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
    
def young_modulus(strain, stress):

    x = strain
    y = stress
    strain = strain[0:2500]
    stress = stress[0:2500]
    
    E, b = [150000, 0]

    hooke = Hooke(strain, E, b)

    model = curve_fit(hooke, strain, stress, p0 = (E, b))

    ans, cov = model
    E, b = ans    
    fit_curve = E * strain + b

    plt.figure(figsize = (16,9))
    plt.plot(x, y, 'b.')
    plt.plot(strain, fit_curve, 'r-' )
    plt.show()
    return E, int(b)
    
