from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math
import random
import collections


#fig = plt.figure()



Tstart=.1
Tstop=10.1

hstart=0
hstop=0
dh=.1


T=Tstart
h=hstart
N=10
qtsteps=1000
stepTermalize=900
dT=0.5


spin=np.zeros((N+1))
meanMag=list()
AbsMag=list()
ene=list()
eneByStep=list()


i=1
spin[0]=0
spin[N]=0
while(i<N):
#+-+-
 # if(i%2==0):
 #  spin[i]=1
 # else:
 #  spin[i]=-1
#todos spins+ -
# spin[i]=1
#random
 rand=random.randint(0,10)
 if(rand<6):
  spin[i]=-1
 else:
  spin[i]=1

 i=i+1



def montecarlo():
 global spin
 i=1
 while(i<N):
#os spins 0 e N nao podem flipar
  posi=random.randint(1,N-1)
  eneBef=eneSys()
  spin[posi]=spin[posi]*(-1)
  eneAft=eneSys()
  deltaE=eneAft-eneBef
  rand=random.random()
  # if(rand<math.e**(-deltaE)):
  if(rand<math.e**(-1/T*deltaE)):
   ene.append(eneSys())
   pass
  else:
   eneAft=eneBef
   spin[posi]=spin[posi]*(-1)
  i=i+1



def eneSys():
 energy=0
 i=1
 while(i<N+1):
  # energy=energy-1/T*(spin[i]*spin[i-1])
  energy=energy-(spin[i]*spin[i-1])
  energy=energy-h*spin[i]
  i=i+1
 return(energy)

T=Tstart
h=hstart
xT=list()







data=list()
mE=0
mM=0
mAM=0

while(h<=hstop):
 j=0
 while(T<=Tstop):
  AbsMag.clear()
  meanMag.clear()
  eneByStep.clear()
  j=0
  while(j<qtsteps):
   AbsMag.append(abs(sum(spin[1:N-1])/(N-2)))
   meanMag.append(sum(spin[1:N-1])/(N-2))
   if(j==stepTermalize): ene.clear()
   eneByStep.append(eneSys())
   montecarlo()
   j=j+1
  mE=sum(eneByStep[stepTermalize:])/len(eneByStep[stepTermalize:])
  mM=sum(meanMag[stepTermalize:])/len(meanMag[stepTermalize:])
  mAM=sum(AbsMag[stepTermalize:])/len(AbsMag[stepTermalize:])
  data.append((h,T,mE,mM,mAM))
  T=T+dT
 h=h+dh


namedir="data/grafs/"
namefile=namedir+"data-N"+str(N)+".dat"

np.savetxt(namefile,data)


# distrib=list(collections.Counter(ene).items())

# xEne=list()
# DistEne=list()

# for value in distrib:
#  xEne.append(value[0])
#  DistEne.append(value[1])

# plt.scatter(xEne,DistEne)
# plt.show()

# xStep=np.linspace(0,len(eneByStep),num=len(eneByStep))

#passo onde a termalizacao ja ocorreu
# print("mean mag:",sum(meanAbsMag[stepTermalize:])/len(meanAbsMag[stepTermalize:]))
# print("mean abs mag:",sum(meanMag[stepTermalize:])/len(meanMag[stepTermalize:]))

# plt.plot(xStep, eneByStep)
# plt.plot(xStep, meanAbsMag, color="blue")
# plt.plot(xStep, meanMag, color="red")

# plt.title("valor medio da energia")
# plt.xlabel("Temperatura")
# plt.ylabel("energia")

# plt.plot(xT, meanEneByT)


#plt.show()

