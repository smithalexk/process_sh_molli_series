# python module to fit (sh)MOLLI data

import math
from scipy.optimize import curve_fit
    
# Define fit function
# This will fit the T1 curve to the data
def fit_func2(x, a, b, t):
    return ( a - b*pow(math.e, -x/t) )

def my_fit(x,y):
    guess = [1000, 1000, 1000]
    params = curve_fit(fit_func2, x, y , p0=guess, method='trf', bounds=(0, float('Inf')), maxfev=10000) #
    [a,b,t] = params[0]
    sse = (sum([fit_func2(xi,a,b,t) for xi in x]))
    return [a,b,t,sse]