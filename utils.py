
import numpy as np 



def Hooke(strain, E = 150000, b = 0):
    return E * strain + b

def r_squared(x, y, model, modelParams):
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

def uniform_plast(strain, stress, sig_y, uts):
    eps = strain.to_numpy()
    sig = stress.to_numpy()

    yield_index = find_index(sig, sig_y)
    uts_index = find_index(sig, uts)

    eps_c = eps[yield_index:uts_index]
    sig_c = sig[yield_index:uts_index]

    return eps_c, sig_c

def true_values(strain, stress, sig_y, uts):
    eps, sig = uniform_plast(strain, stress, sig_y, uts)
    
    eps_t = np.log(1 + eps)
    sig_t = sig * (1 + eps)

    return eps_t, sig_t


def log_Hollomon(strain, stress, K = 0.25, N = 300):
    
    return np.log(K) + N * np.log(strain)
    
