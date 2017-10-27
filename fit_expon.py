# --- Fit exponential curve
import numpy as np
from scipy.optimize import curve_fit

def main(xvalues,yvalues):
    popt, pcov = curve_fit(exponenial_func, xvalues, yvalues, p0=(1, 1, yvalues[-1]))
    yexp = []
    for value in xvalues:
        yexp.append(popt[0]*np.exp(-popt[1]*value)+popt[2])
    fitted = {}
    fitted['params'] = popt
    fitted['values'] = yexp
    return fitted
    
def exponenial_func(x, a, b, c):
    return a*np.exp(-b*x)+c
