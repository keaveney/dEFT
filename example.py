from tools import splash
from tools import scanner
from asciimatics.screen import Screen
import numpy as np

#intro graphics
Screen.wrapper(splash.deft_splash)

data = [1.0, 2.0, 3.0]
invcov = [[1.5, -0.49, -0.5], [-0.49, 1.5,  -0.5], [-0.50, -0.5,1.5]]

#configure scan
coefficients = ["ctg"]

pred_basis = {'SM': np.array([2.0,3.0,4.0]), 'ctg-': np.array([1.5,2.5,3.5]), 'ctg+': np.array([2.5,3.5,4.5]) }

sc  = scanner.scanner()
sc.init(data,invcov,pred_basis)

scan_params = [-1.0, 1.0, 21]

paramspace = np.linspace(scan_params[0],scan_params[1],scan_params[2])

for scenario in paramspace:
    chi_sq = sc.eval(scenario)
    print "scenario = " + str(scenario) + "chi sqared equals = " + str(chi_sq)


