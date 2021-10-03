from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math
import random
import collections
import os

def loadtxt(filename):
	txtfile=open(filename, "r")
	data=txtfile.readlines() #[70:160]
	X=list()
	Y=list()
	for line in data:
		line=line.replace("\n","")
		linesplit=line.split(" ")
		# T=float(linesplit[0])
		x=float(linesplit[0])
		# beta=1./T
		# beta=float(linesplit[0])
		# p=1.0-math.e**(-2*beta)

		# X.append(beta)
		y=float(linesplit[1])
		if(y!=0.0 or x!=0.0):
			X.append(x)
			Y.append(y)
	return(X,Y)

# x0,y0=loadtxt("data/mag-100x100.dat")
# x0,y0=loadtxt("data/mag.dat")
x1,y1=loadtxt("data/mag-10x10-T1.3.dat")
x2,y2=loadtxt("data/mag-30x30-T1.3.dat")
x3,y3=loadtxt("data/mag-100x100-T1.3.dat")

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(x1, y1, color="blue", label="N=10x10, bound=0");
ax.plot(x2, y2, color="red", label="N=30x30, bound=0");
ax.plot(x3, y3, color="green", label="N=100x100, bound=0");


plt.grid()
plt.legend()
plt.title("Magnetização em função de h (T=1.3)")
plt.xlabel("h")
plt.ylabel("Magnetização")
plt.show()
plt.clf()


