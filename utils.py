
import numpy as np 



def Hooke(strain, E = 150000, b = 0):
    return E * strain + b

def r_squared(x, y, model, modelParams):
    res = y - model(x, *modelParams)
    ss_res = np.sum(res ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res/ss_tot)
    return r_squared