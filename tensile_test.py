import matplotlib.pyplot as plt 
import numpy as np 
import scipy as sp 
from scipy.optimize import curve_fit
from utils import Hooke, r_squared


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
    r2 = r_squared(x, y, Hooke, ans)
# plotting the figure and showing or saving the figure
    plt.figure(figsize = (8,4.5))
    plt.plot(strain, stress, 'b-')
    plt.plot(x, fit_curve, 'r-', linewidth = 2 )
    plt.xlabel('strain [mm/mm]')
    plt.ylabel('stress [MPa]')
    plt.xlim(0, 1.05* max(strain))
    plt.ylim(0, 1.05*max(stress))
    plt.text(0.1*max(strain), 0.1*max(stress), f'The elasticity modulus is {E_gpa} GPa, R² = {round(r2, 4)}', fontsize=12)  

    if save != False:     
        plt.savefig('output/elasticity', dpi=300, bbox_inches='tight',transparent=False)
    else:
        plt.show()
    
    return E, int(b), r2

def sigma_y(strain, stress, E_mpa, b, save = False):
    # taking part of the data
    x = strain[strain < 0.005]
    k = x - 0.002
    z = Hooke(k, E_mpa)
    # finding the index of sig_y
    i = 0
    while stress[i] > z[i]:
        i = i + 1 

    sig_y = stress[i]


    plt.figure(figsize = (8,4.5))
    plt.plot(strain, stress, 'b-')
    plt.plot(x, Hooke(x, E_mpa, b), 'r--', linewidth = 1 )
    plt.plot(x, Hooke(k, E_mpa), 'r:', linewidth = 1 )
    plt.xlabel('strain [mm/mm]')
    plt.ylabel('stress [MPa]')
    plt.xlim(0, 1.05* max(strain))
    plt.ylim(0, 1.05*max(stress))
    plt.text(0.1*max(strain), 0.1*max(stress), f'The yield strength is {round(sig_y)} MPa.', fontsize=12)

    if save != False:     
        plt.savefig('output/sigma_yield', dpi=300, bbox_inches='tight',transparent=False)
    else:
        plt.show()
    
    print(f'The yield strength is {round(sig_y)} MPa.')
    return sig_y

def UTS(strain, stress):
    return max(stress)
    
