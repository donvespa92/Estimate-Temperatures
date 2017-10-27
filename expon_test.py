# Scientific libraries
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


x = np.array([399.75, 989.25, 1578.75, 2168.25, 2757.75, 3347.25, 3936.75, 4526.25, 5115.75, 5705.25])
y = np.array([109,62,39,13,10,4,2,0,1,2])

def exponenial_func(x, a, b, c):
    return a*np.exp(-b*x)+c


popt, pcov = curve_fit(exponenial_func, x, y, p0=(1, 1e-6, 1))


yexp = np.array
for value in x:
    np.append(yexp,popt[0]*np.exp(-popt[1]*value)+popt[2])