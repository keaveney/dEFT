from tools.predBuilder import predBuilder
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

testSamples = [[1.0, -49.0, -46.0], [1.0, -41.0, 45.0], [1.0, -43.0, 40.0], [1.0, 49.0, -50.0], [1.0, 10.0, -46.0], [1.0, -49.0, 50.0], [1.0, 10.0, -11.0], [1.0, 20.0, -25.0], [1.0, 20.0, -20.0], [1.0, -10.0, -10.0], [1.0, 49.0, -48.0], [1.0, -200.0, -200.0]]
testPreds = np.asarray([[278], [169.2], [156.2], [192.2], [105.4],[211.9], [26.48],[52.67], [46.06], [29.72], [183.2], [4738]])


prior_limits ={
    "cQq83": [-50.0, 25.0],
    "cQq81": [-40.0, 40.0]
}

#preds = np.asarray( [[189.4], [182.0], [302], [258.6], [184.2], [290.1]])
preds = np.asarray( [[189.7], [180.5], [299.8], [258.7], [183.4], [20.46]])

samples = [[1.0, -47.0, 46.0], [1.0, -50.0, 40.0], [1.0, 47.0, 47.0], [1.0, -48.0, -44.0], [1.0, -47.0, 44.0], [1.0, 0.0, 0.0]]

pb = predBuilder()
pb.init(len(prior_limits), samples, preds)

vMakePred = np.vectorize(pb.makePred)

ciAr = np.arange(-100, 100, 5.0)
cjAr = np.arange(-100, 100, 5.0)


ciP = ([])
cjP = ([])
xs = ([])

#z = np.array([pb.makePred(([ci, cj])) for ci in ciAr for cj in cjAr])

z = np.array([])

for ci in ciAr:
    for cj in cjAr:
        z = np.append(z, pb.makePred([ci, cj]))
        print("ci = " + str(ci) + " cj = " + str(cj) + " z = " + str(pb.makePred([ci, cj])))


z = z.reshape(len(cjAr),len(ciAr))

print("z shape = " + str(z.shape))

print("z max = " + str(np.amax(z)))
print("z min = " + str(np.amin(z)))

cols = np.arange(np.amin(z), np.amax(z), 1)
contour_filled = plt.contourf(cjAr, ciAr, z, levels=cols)

contour = plt.contour(cjAr, ciAr, z, levels=[0.0])

cbar = plt.colorbar(contour_filled)
cbar.set_label(r"$\sigma_{t\bar{t}}$ [pb]", rotation=90)

plt.xlabel(r"$cQq83$")
plt.ylabel(r"$cQq81$")

plt.show()


#Z = z.reshape(len(ciAr), len(ciAr))
#plt.imshow(Z, interpolation='bilinear')
#plt.savefig("test.png")


#x = data[:, 0]
#y = data[:, 1]
#z = data[:, 2]

#for ci in ciAr:
#    for cj in cjAr:
#        cPoint = ([ci, cj])
##        ciP = np.append(ciP,ci)
#        cjP = np.append(cjP,cj)
#        xs = np.append(xs,pb.makePred(cPoint) )
        #print(str(ci) + ", " + str(cj) + ", " + str(pb.makePred(cPoint)) )

#fig = plt.figure(figsize=(8, 6))
#ax = fig.add_subplot(111, projection="3d")
#ax.plot_trisurf(ciP,cjP,xs)

#cols = np.arange(0, 85, 1)
#contour_filled = plt.contourf(ciP, cjP, xs, levels=cols)
#cbar = plt.colorbar(contour_filled)
#cbar.set_label(r"$\sigma_{t\bar{t}}$ [pb]", rotation=90)

#plt.xlabel(r"$c_tG$")
#plt.ylabel(r"$c^{3,8}_Qq$")

#plt.show()
