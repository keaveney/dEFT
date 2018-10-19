from tools import splash
from tools import scanner
from asciimatics.screen import Screen
import numpy as np
import itertools
import matplotlib as plt


#intro graphics
Screen.wrapper(splash.deft_splash)

data = [1.0, 2.0, 3.0]
invcov = [[1.5, -0.49, -0.5], [-0.49, 1.5,  -0.5], [-0.50, -0.5,1.5]]

max_coeff_power = 1

#configure scan
coefficients = ["ctg", "ctw"]

preds = {'SM': np.array([2.0,3.0,4.0]),
              'ctg-': np.array([1.5,2.5,3.5]),
              'ctg+': np.array([2.5,3.5,4.5]),
              'ctw-': np.array([1.2,2.0,2.4]),
              'ctw+': np.array([2.1,3.2,4.1])}

sc  = scanner.scanner()
sc.init(data,invcov,preds,coefficients,max_coeff_power)

scan_params = {
        'ctg': np.linspace(-1.0, 1.0, 5),
        'ctw': np.linspace(-1.0, 1.0, 3),
}

ranges = [ ]

for coefficient in coefficients:
    ranges.append(scan_params[coefficient].tolist())

chi_sq_array = np.array()

# build parameter space
for scenario in itertools.product(*ranges):
    chi_sq_array.append(sc.eval(scenario))
    #print "scenario = " + str(scenario) + "chi sqared equals = " + str(chi_sq)


chi_sq_array.reshape(len(coefficients),len(coefficients) )

plt.contourf(scan_params['ctg'], scan_params['ctw'], chi_sq_array, 8, alpha=.75, cmap='jet')


#C = plt.contour(X, Y, f(X, Y), 8, colors='black', linewidth=.5)



