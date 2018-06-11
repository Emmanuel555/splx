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
while(i<len(l)):
    type = l[i]
    xx = float(l[i+1].strip())
    yy = float(l[i+2].strip())
    i+=3
    if(type == 'c'):
        plt.plot([xx], [yy], marker='o')
    else:
        x.append(xx)
        y.append(yy)
plt.plot(x, y)
plt.show()
