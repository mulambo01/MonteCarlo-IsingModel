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
	data=txtfile.readlines()
	X=list()
	Y=list()
	for line in data:
		line=line.replace("\n","")
		linesplit=line.split(" ")
		beta=float(linesplit[0])
		p=1.0-math.e**(-2*beta)
		X.append(p)
		# Y.append(abs(float(linesplit[1])))
		Y.append(abs(float(linesplit[1])))
	return(X,Y)

x0,y0=loadtxt("data/mag.dat")
x1,y1=loadtxt("data/mag-b010.dat")
x2,y2=loadtxt("data/mag-b050.dat")
x3,y3=loadtxt("data/mag-b0100.dat")
x4,y4=loadtxt("data/mag-b+10.dat")
x5,y5=loadtxt("data/mag-b+50.dat")
x6,y6=loadtxt("data/mag-b+100.dat")
x7,y7=loadtxt("data/mag-b-10.dat")
x8,y8=loadtxt("data/mag-b-50.dat")
x9,y9=loadtxt("data/mag-b-100.dat")

fig = plt.figure()
ax = fig.add_subplot(111)

# ax.plot(x1, y1, color="blue", label="N=10, bound=0");
# ax.plot(x4, y4, color="orange", label="N=10, bound=+");
# ax.plot(x7, y7, color="green", label="N=10, bound=-");

# ax.plot(x2, y2, color="blue", label="N=50, bound=0");
ax.plot(x5, y5, color="purple", label="N=50, bound=+");
# ax.plot(x8, y8, color="green", label="N=50, bound=-");

# ax.plot(x3, y3, color="blue", label="N=100, bound=0")
ax.plot(x6, y6, color="red", label="N=100, bound=+")
# ax.plot(x9, y9, color="green", label="N=100, bound=-")

# ax.plot(x0, y0, color="blue", label="1d N=1000, bound=+")

plt.legend()
plt.title("Magnetização absoluta em função de p")
plt.xlabel("p=1-e^-2β")
plt.ylabel("Magnetização absoluta")
plt.show()
plt.clf()


