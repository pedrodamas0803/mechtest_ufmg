import matplotlib.pyplot as plt 
import numpy as np 
import scipy as sp 
# from utils import *


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
    
# def young_modulus(strain, stress):
#     x = strain.to_numpy()
#     y = stress.to_numpy()
    
