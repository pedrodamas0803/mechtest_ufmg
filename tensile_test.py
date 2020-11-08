import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps, trapz
from scipy.optimize import curve_fit
from utils import Hooke, r_squared, find_index


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



def young_modulus(strain, stress, show_plot = False, save = False):

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

    if show_plot == True:
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

def sigma_y(strain, stress, E_mpa, b = 0, show_plot = True, save = False):
    # taking part of the data
    x = strain[strain < 0.005]
    k = x - 0.002
    z = Hooke(k, E_mpa)
    # finding the index of sig_y
                                        # TODO: replace find_index
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

    if save == True:
        plt.savefig('output/sigma_yield', dpi=300, bbox_inches='tight',transparent=False)
    
    if show_plot== True:
        plt.show()

    print(f'The yield strength is {round(sig_y)} MPa.')
    return sig_y

def UTS(strain, stress, show = True):
    uts = max(stress)
    if show == True:
        print(f'The ultimate tensile strength is {round(uts)} MPa. ')
    return uts

def uniform_elong(strain, stress, show = True):
# transform data to numpy arrays
    eps = strain.to_numpy()
    sig = stress.to_numpy()

# finds the index of the maximum stress
    index = np.argmax(sig)
    u_elong = eps[index]

# prints depending on the user's choice
    if show == True:
        print(f'The uniform elongation is {round(u_elong, 4)}.')

    return u_elong

def aprox_resilience(E, sig_y, show = True):

    U_r = (sig_y ** 2)/(2*E)

    if show == True:
        print(f'The resilience calculated by the formula yield²/2E is {round(U_r,2)} MJ/m³. ')

    return U_r

def resilience(strain, stress, sig_y, dx = 1.0, show = True):
    
    i = find_index(stress, sig_y)
    eps = strain[0:i]
    sig = stress[0:i]
    U_r = trapz(eps, sig, dx)

    if show == True:
        print(f'The resilience calculated by trapezoidal integration with dx = {dx} is {round(U_r,2)} MJ/m³. ')

    return U_r


def aprox_toughness(strain, sig_y, uts, show = True):

    eps_f = max(strain)

    U_t = eps_f * ((sig_y + uts) / 2)
    
    if show == True:
        print(f'The material toughness is approximately {round(U_t)} MJ/m³.')
    return U_t 


def toughness(strain, stress, show = True):

    mat_toughness = trapz(stress, strain)

    if show == True:
        print(f'The material toughness computed by Trapezoidal method is {round(mat_toughness)} MJ/m³. ')
    return mat_toughness

def plot_flow_curve(strain, stress, sig_y, uts, show_plot = True, save = False):
                                                    # TODO: deal with discontinuous yielding; how can I determine where it ends?
    eps = strain.to_numpy()
    sig = stress.to_numpy()

    yield_index = find_index(sig, sig_y)
    uts_index = find_index(sig, uts)

    eps_c = eps[yield_index:uts_index]
    sig_c = sig[yield_index:uts_index]

    plt.figure(figsize = (8,4.5))
    plt.plot(eps_c, sig_c, 'b-')
    plt.xlabel('strain [mm/mm]')
    plt.ylabel('stress [MPa]')
    plt.xlim(0, 1.05* max(eps_c))
    plt.ylim(0.95 * sig_c[0], 1.05* max(sig_c))

    if show_plot == True:
        plt.show()

    if save == True:
        plt.savefig('output/flow_curve', dpi=300, bbox_inches='tight',transparent=False)

def flow_model(strain, stress, func, *func_params):
    pass
    