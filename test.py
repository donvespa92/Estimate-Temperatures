import fit_expon
import numpy as np
import matplotlib.pyplot as plt

x = np.array([399.75, 989.25, 1578.75, 2168.25, 2757.75, 3347.25, 3936.75, 4526.25, 5115.75, 5705.25])
y = np.array([109,62,39,13,10,4,2,0,1,2])

fitted = {}
fitted = fit_expon.main(x,y)

plt.plot(x,fitted['values'])
