import numpy as np
import matplotlib.pyplot as plt

# the n values of c_i that define the 'input' distributions: c_sm, c_i1...c_in
# there are Nops columns.
#first column is the SM coupling squared values for each of
# the input samples,
# second is the squared value of the second term of the cross section expansion
# and so on for the third and subsequent solumns
# The convention here is that after the SM term, the squared term for the
# c_i comes next, then the ci-SM interference term, before continuing on
# for c_i+1 etc.

cInputAr = np.array([[1.0, 1.0, 1.0], [0.0, 4.0, 4.0], [0.0, 2.0, -2.0]])
inWeightsAr = np.linalg.inv(cInputAr)

Tin = np.array([ [20.76], [34.48], [12.71]])
TinT = Tin.transpose()

# 'sanity check' example
y = np.array([])
x = np.array([])

for target in range(-10,10):
    ciTargetAr = np.array([[1], [target**2.0], [target]])
    morphWeightsAr = inWeightsAr.dot(ciTargetAr)
    Tout = TinT.dot(morphWeightsAr)
    print("target = " + str(target) + " tout " + str(Tout))
    x = np.append(x, target)
    y = np.append(y, Tout)


valPtsX = [-3.31, 5.72]
valPtsY = [11.17, 77.03]


inPtsX = [0.0, 2.0, -2.0]
inPtsY = [20.76, 34.38, 12.71]

fig, ax = plt.subplots()
ax.set_xlabel('ctg')
ax.set_ylabel('cross section [pb]')

plt.plot(x,y, label='morphing model')
plt.plot(inPtsX, inPtsY, "b+", markersize=7, label='input predictions')
plt.plot(valPtsX, valPtsY, "r+", markersize=7,label='validation predictions')

legend = plt.legend(loc='upper center', shadow=True, fontsize='x-large')

plt.show()
