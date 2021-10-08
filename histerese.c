#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include <math.h>
#include <string.h>


mtx genSpins(int nrow, int ncol, float bound){
 int Nrow=nrow+2;
 int Ncol=ncol+2;
 int value;
 float random;
 mtx Spins = crystalmatrix(Nrow,Ncol,bound);
 for(int i=1; i<Nrow-1; i++){
  for(int j=1; j<Ncol-1; j++){
   random=rand()%10000/10000.0;
   if(random<0.5){
    value=1;
   }
   else{
    value=-1;
   }
   Spins.data[i][j]=value;
  }
 }
 return(Spins);
}


float eneSys(mtx spins, float T, float h){
 float energy=0.0;
 int Nrow=spins.nrows;
 int Ncol=spins.ncols;
 int i,j;
 for(i=1;i<Nrow-1;i++){
  for(j=1;j<Ncol-1;j++){
   energy=energy-1.0/T*spins.data[i][j]*(spins.data[i-1][j]+spins.data[i][j-1]);
   // energy=energy-spins.data[i][j]*(spins.data[i-1][j]+spins.data[i][j-1]);
   energy=energy-h*spins.data[i][j];
  }
 }
 for(i=1;i<Ncol-1;i++){
  energy=energy-1.0/T*spins.data[Nrow-1][i]*spins.data[Nrow-2][i];
  // energy=energy-spins.data[Nrow-1][i]*spins.data[Nrow-2][i];
 }
 for(i=1;i<Nrow-1;i++){
  energy=energy-1.0/T*spins.data[i][Ncol-1]*spins.data[i][Ncol-2];
  // energy=energy-spins.data[i][Ncol-1]*spins.data[i][Ncol-2];
 }
 return(energy);
}

void montecarlo(mtx spins,float T, float h){
 int Nrow=spins.nrows-2;
 int Ncol=spins.ncols-2;
 int row, col, posi;
 float eneAft, eneBef, random, deltaEne, value;
 eneBef=eneSys(spins,T,h);
 for(int i=0;i<Nrow*Ncol;i++){
  posi=rand()%(Nrow*Ncol);
  row=posi%(Nrow)+1;
  col=(int) (posi-row)/(Nrow)+1.1;

  // row=rand()%Nrow+1;
  // col=rand()%Ncol+1;

  eneBef=eneSys(spins, T,h);
  spins.data[row][col]=-spins.data[row][col];
  eneAft=eneSys(spins,T,h);
  deltaEne=eneAft-eneBef;
  random=rand()%1000/1000.0;
  if(random<expf(-deltaEne)){
  // if(random<expf(-1.0/T*deltaEne)){
   // eneBef=eneAft;
   // printf("aceito\n");
  }
  else{
   // eneAft=eneBef;
   spins.data[row][col]=-spins.data[row][col];
  }
 }
}

float magnet(mtx spins){
 int Nrow=spins.nrows-2;
 int Ncol=spins.ncols-2;
 int i,j;
 // float* data=(float *)malloc(2*sizeof(float));
 float mag=0.0, absMag=0.0;
 int qt=0;
 for(i=1;i<Nrow+1;i++){
  for(j=1;j<Ncol+1;j++){
   mag=mag+spins.data[i][j];
  }
 }
 return(mag/(Nrow*Ncol));
}
void setSpins(mtx spins, float value){
 int nrow=spins.nrows;
 int ncol=spins.ncols;
 int i, j;
 for(i=1; i<nrow-1; i++){
  for(j=1; j<ncol-1; j++){
   spins.data[i][j]=value;
  }
 }
}

void setBound(mtx spins, float* bound){
 float left=bound[0], right=bound[1], down=bound[2], up=bound[3];
 int i;
 int ncol=spins.ncols, nrow=spins.nrows;
 for(i=1; i<nrow-1; i++){
  spins.data[i][0]=left;
  spins.data[i][ncol-1]=right;
 }
 for(i=1; i<ncol-1; i++){
  spins.data[0][i]=up;
  spins.data[nrow-1][i]=down;
 }
}

int main(){
 srand(91);
 int Nrow=100, Ncol=100;
 int qtsteps=10;
 float T,h, hmax, hmin, dh;
 float meanMag;
 char filename[100];
 int qtdata, qtMean;

 mtx Spins=genSpins(Nrow, Ncol,0.0), oldSpins;
 // setSpins(Spins, 1.0);
 float bound[4]={0.0, 0.0,0.0,0.0};
 setBound(Spins, bound);
 int i,count;
 h=0.0;

 hmin=-1.0;
 hmax=1.0;
 dh=0.05;
 qtdata=1000;

 mtx mag=nullmatrix(qtdata, 2);
 count=0;
 qtMean=5;
 T=5.3;
 for(h=0.0; h<hmax; h=h+dh){
  for(i=0;i<qtsteps;i++){
   montecarlo(Spins, T,h);
  }
  meanMag=0.0;
  for(i=0;i<qtMean;i++){
   meanMag=meanMag+magnet(Spins);
   montecarlo(Spins, T,h);
  }
  meanMag=meanMag/qtMean;
  printf("%d %f %f \n", count,h,meanMag);
  mag.data[count][1]=meanMag;
  mag.data[count][0]=h;
  count++;
 }

 for(; h>hmin; h=h-dh){
  for(i=0;i<qtsteps;i++){
   montecarlo(Spins, T,h);
  }
  meanMag=0.0;
  for(i=0;i<qtMean;i++){
   meanMag=meanMag+magnet(Spins);
   montecarlo(Spins, T,h);
  }
  meanMag=meanMag/qtMean;
  printf("%d %f %f \n", count,h,meanMag);
  mag.data[count][1]=meanMag;
  mag.data[count][0]=h;
  count++;
 }

 for(; h<hmax; h=h+dh){
  for(i=0;i<qtsteps;i++){
   montecarlo(Spins, T,h);
  }
  meanMag=0.0;
  for(i=0;i<qtMean;i++){
   meanMag=meanMag+magnet(Spins);
   montecarlo(Spins, T,h);
  }
  meanMag=meanMag/qtMean;
  printf("%d %f %f \n", count,h,meanMag);
  mag.data[count][1]=meanMag;
  mag.data[count][0]=h;
  count++;
 }

 sprintf(filename, "data/mag.dat");
 mtxsave(filename, mag);

}
