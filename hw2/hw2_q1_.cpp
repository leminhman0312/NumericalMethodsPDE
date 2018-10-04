//Question 1 

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

//DECLARE STRUCTURE

struct Hyperbolic
{
	double dx;
	double dt;
	int imax;
	double tmax; //max time
	double nmax; //max time steps
	double alpha; //speed of sound
	double pi;
	vector <double> u0; //initial u //USED IN LOOPS
	vector <double> u; //u at n+1 level, USED IN LOOPS
	vector <double> ufull; //end solution, for print
	vector <double> uInitial; //initial sol, for print 
	vector <double> uAnalytical;
 
};

void InitializeU(struct Hyperbolic *wave){

	//BASIC PROPERTIES	
	wave->dx = 0.01;
	wave->dt = 0.5*wave->dx;
	wave->imax = 101;
	wave->tmax = 4.0;
	wave->nmax = (wave->tmax/wave->dt);
	wave->pi = 4.0*atan(1.0);
	
	//SETTING U AND UN
	wave->u0.resize(wave->imax);
	wave->u.resize(wave->imax);
	wave->ufull.resize(wave->imax);
	wave->uInitial.resize(wave->imax);
	wave->uAnalytical.resize(wave->imax);
	

	//SETTING UP INITAL CONDITIONS, U0
	for (int i = 0; i<wave->imax;i++){
		
		wave->u[i] = sin(4.0*wave->pi*(i)*wave->dx);		
	}	
}


void LaxWendroff(struct Hyperbolic *wave){	

	int imax = wave->imax; 
	int nmax = wave->nmax;
	//Finite Difference Loop 

	//Time loop, FULL TIME
	for (int n = 0; n<=nmax;n++){
		//copy loop 
	
		
		//FDE
		for (int i = 1; i<imax-1;i++){
			
			wave->u0[i] = 0.25*(wave->u[i+1]-wave->u[i-1])+ 0.5*(wave->u[i+1]+wave->u[i-1]);		

		}

		//BC 

		//Periodic BC

		wave->u0[0] = 0.25*(wave->u[1]-wave->u[imax-2])+ 0.5*(wave->u[1]+wave->u[imax-2]);
		wave->u0[imax-1] = wave->u0[0];

		for (int i = 0; i<imax;i++){
			wave->u[i] = wave->u0[i];			
		}	


	}


	for (int i = 0; i<imax;i++){
		printf("%1.5f\n",wave->u0[i]);
	}





}


void runLaxWendroff(){
	Hyperbolic wave;
	InitializeU(&wave);
	LaxWendroff(&wave);
}

int main(){
	runLaxWendroff();
	
}
