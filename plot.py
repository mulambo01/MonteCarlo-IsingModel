from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math
import random
import collections
import os

k=0
size=51
spins=size*[list()]
# path="data/"
# files=os.listdir(path)
# for name in files:
# 	if(name.count("spins-h")):



while(k<size):
	txtfile=open("data/spins-T-"+str(k)+".dat", "r")
	data=txtfile.readlines()
	spins[k]=list()
	for line in data:
		line=line.replace("\n","")
		linesplit=line.split(" ")
		i=0
		for data in linesplit:
			linesplit[i]=float(data)
			i=i+1
		spins[k].append(linesplit)
	k=k+1

# print(spins[2])
qt=len(spins[0])
print(qt)
x=np.linspace(0,qt-1,qt)
y=np.linspace(0,qt-1,qt)
X,Y=np.meshgrid(x, y)
k=0
fig = plt.figure()
while(k<=size):
	ax = fig.add_subplot(111)
	data2=np.array(spins[k])
	plt.xlim(-.5,len(data2)-1+.5)
	plt.ylim(-.5,len(data2)-1+.5)
	plt.ion()
	plt.title("T="+str(k))
	ax.scatter(X, Y, data2+2,c=data2+2, cmap='coolwarm', linewidth=5, marker='o');
	# ax.set_xticklabels([])
	# ax.set_yticklabels([])
	# fig.colorbar(a)
 # ax.set_zticklabels([])
 # ax.scatter(X,Y,spins, cmap='warm')
	plt.show()
	plt.pause(0.5)
	k=k+1
	plt.clf()


