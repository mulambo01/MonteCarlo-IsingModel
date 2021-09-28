from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math
import random
import collections
import os


txtfile=open("data/ene.dat", "r")
data=txtfile.readlines()
X=list()
Y=list()
for line in data:
	line=line.replace("\n","")
	linesplit=line.split(" ")
	X.append(float(linesplit[0]))
	Y.append(float(linesplit[1]))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(X, Y);
plt.show()
plt.clf()


