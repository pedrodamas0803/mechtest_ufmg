
import numpy as np 
import matplotlib.pyplot as plt 
import os
'''
Auxiliary functions module.
'''



def Hooke(strain, E = 150000, b = 0):
        '''
        Defines the linear equation for the linear portion of the stress-strain curve.

        Inputs:

        strain - the vector containing the strain values
        E - default = 150,000; initial guess for the modulus of elasticity in MPa.
        b - default = 0; intercept of the liear equation.

        Output:

        Array with the calculated values of E * strain + b.
        '''

        return E * strain + b

def r_squared(x, y, model, modelParams):

        '''
        Calculates the R² statistics to evaluate the quality of the data regression.

        Inputs:

        x - the vector containing the x values of the independent variable of the function.
        y - the vector containing the dependent (measured) values as a function of x.
        model - the function that receives x as variable and calculates a y vector to model measured y.
        modelParams - the parameters of the model function passed as a list.

        Output:

        r_squared - the R² statistics. 
        '''

        res = y - model(x, *modelParams)
        ss_res = np.sum(res ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r_squared = 1 - (ss_res/ss_tot)
        
        return r_squared

def find_index(array, value):

        '''
        Finds the index of a value in an array.
        '''
    
        i = 0
        while array[i] < value:
                i = i + 1
        return i
  


def Hollomon(x, K = 300, n = 0.25):

        '''
        Defines the function to perform flow model fitting.
        '''

        return K * x ** n

def Ludwik(x, sig_o = 600, K = 600, n = 0.24):

        '''
        Defines the function to perform flow model fitting.
        '''

        return sig_o + K * x ** n

def Datsko(x, K = 300, x0 = 300, n = 0.20):

        '''
        Defines the function to perform flow model fitting.
        '''

        return K * (x0 + x) ** n 
      
def create_output_folder(path = 'output'):
        '''
        Creates the output folder where the results will be saved.
        '''

        if not os.path.exists(path):
            os.mkdir(path)
        else:
            print('The path already exists.')