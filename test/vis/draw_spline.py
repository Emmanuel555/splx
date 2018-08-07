import json
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('splinefile', help='json file that contains the spline')
args = parser.parse_args()

f = open(args.splinefile, "r")
l = [x[:-1] for x in f.readlines()]
f.close()
print(l)

i = 0
x = []
y = []

xv = []
yv = []
dvx = []
dvy = []
while(i<len(l)):
    type = l[i]
    xx = float(l[i+1].strip())
    yy = float(l[i+2].strip())
    i+=3
    if(type == 'c'):
        plt.plot([xx], [yy], marker='o')
    elif(type == 'v'):
        xv.append(x[-1])
        yv.append(y[-1])
        dvx.append(xx)
        dvy.append(yy)
    elif(type == 'd'):
        x.append(xx)
        y.append(yy)
    elif(type == 'd2'):
        pass
    else:
        raise "no"
plt.plot(x, y)
plt.plot(x, y, marker='o', ms=2)
plt.quiver(xv, yv, dvx, dvy, width = 0.0025, angles = 'xy')

#ymin, ymax = plt.ylim()
#xmin = -ymin/2
#xmax = -ymax/2
#plt.plot([xmin, xmax], [ymin, ymax])


plt.show()
