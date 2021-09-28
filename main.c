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
   	value=1;
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
	for(i=1;i<Nrow;i++){
		for(j=1;j<Ncol;j++){
			energy=energy-1/T*spins.data[i][j]*(spins.data[i-1][j]+spins.data[i][j-1]);
			energy=energy-h*spins.data[i][j];
		}
	}
	for(i=1;i<Ncol;i++){
		energy=energy-1/T*spins.data[Nrow-1][i]*spins.data[Nrow-2][i];
	}
	for(i=1;i<Nrow;i++){
		energy=energy-1/T*spins.data[i][Ncol-1]*spins.data[i][Ncol-2];
	}
	return(energy);
}

void montecarlo(mtx spins,float T, float h){
 int Nrow=spins.nrows-2;
 int Ncol=spins.ncols-2;
 int row, col;
 float eneAft, eneBef, random, deltaEne, value;
 eneBef=eneSys(spins,T,h);
 for(int i=0;i<Nrow*Ncol;i++){
 	row=rand()%Nrow+1;
 	col=rand()%Ncol+1;

 	// col=row;
  eneBef=eneSys(spins, T,h);
  spins.data[row][col]=-spins.data[row][col];
  eneAft=eneSys(spins,T,h);
  deltaEne=eneAft-eneBef;
  random=rand()%100000/100000.0;
  if(random<expf(-deltaEne)){
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
void reUpSpins(mtx spins){
 int nrow=spins.nrows;
	int ncol=spins.ncols;
 int i, j;
	for(i=1; i<nrow; i++){
		for(j=1; j<ncol; j++){
			spins.data[i][j]=1.0;
		}
	}
}

int main(){
srand(4);
int Nrow=100, Ncol=100, qtsteps=40, stepTerm=90;
float Tstart=0.01, Tstop=10.1, dT=0.2, T;
float hstart=-5.0, hstop=5.0, dh=1.0,h=0.0;
mtx Spins=genSpins(Nrow, Ncol,1.0), oldSpins;
int i;
int j=0;
char filename[100];
int jlim=50;
mtx ene=nullmatrix(qtsteps, 2);
mtx mag=nullmatrix(jlim, 2);

float beta, betaStart=0.01, dbeta=0.02;
beta=betaStart;
for(j=0;j<jlim;j++){
	reUpSpins(Spins);
 T=1.0/beta;
	for(i=0;i<qtsteps;i++){
		// printf("%d\n", i);
		montecarlo(Spins, T,h);
		// ene.data[i][1]=eneSys(Spins,T,h);
	}
		mag.data[j][1]=magnet(Spins);
		printf("%f\n", mag.data[j][1]);
		mag.data[j][0]=beta;
		beta=beta+dbeta;
}	// sprintf(filename, "data/ene.dat");
	sprintf(filename, "data/mag.dat");
 // mtxsave(filename, ene);
 mtxsave(filename, mag);

 // printf("nome %.0f\n", T);
	// sprintf(filename, "data/spins-h-%.0f.dat",h);
	// printf("%s\n", filename);
 // mtxsave(filename, Spins);
 // memset(filename,0,sizeof(filename));
// }





}