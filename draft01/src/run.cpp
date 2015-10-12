/* MAIN PROGRAM */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include "../include/shoot.h"
#include "../include/read.h"
#include "../include/write.h"

void init(int, int, int, int, double&, double&, double*, double*);
void getAreaVlme(int, int, double&, double&);

int main(){
	int    i, j, k, l;
	int    n    = 9;		// size of solution vector
	int    m    = 1001;	// number of shooting points
	int    nrk;					// number of Runge-Kutta steps
	int    v;						// reduced volume
	int    rmax;				// maximum radius
	int    flag;				// error flag
	bool   info;				// boolean for file check

	// output directory
	string opath = "../output";

	// abscissa and solution vectors
	vector<double> t(m), si(n*m), sm(n*m), sf(n*m);
	
	// parameters
	double Ca, area, vlme;
	double dC0, dC1, slope;
	vector<double> vecCa;
	vector<int   > vecVr;
	vector<int   > vecRm;

	// range of capillary numbers
	vecCa.push_back(0.0);

	// range of reduced volumes
	vecVr.push_back(95);
	vecVr.push_back(90);
	vecVr.push_back(85);
	vecVr.push_back(80);
	vecVr.push_back(75);

	// range of maximum (projected) radii
	vecRm.push_back(90);
	//vecRm.push_back(91);
	vecRm.push_back(92);
	//vecRm.push_back(93);
	vecRm.push_back(94);
	//vecRm.push_back(95);
	vecRm.push_back(96);
	//vecRm.push_back(97);
	//vecRm.push_back(98);
	//vecRm.push_back(99);
	
	// loop over reduced volume
	for (l = 0; l < vecVr.size(); l++){
		v = vecVr[l];

		// loop over maximum radius
		for (k = 0; k < vecRm.size(); k++){
			rmax = vecRm[k];
	
			// initialize
			nrk = 30;
			cout << "Initializing... " << endl;
			init(n, m, v, rmax, area, vlme, t.data(), si.data());
			for (j = 0; j < m*n; j++){
				sm[j] = si[j];
			}
			cout << "Initialization complete." << endl;
			
			// loop over capillary numbers
			for (i = 0; i < vecCa.size(); i++){
				if (Ca > 50)
					nrk = 60;

				// choose capillary number
				Ca = vecCa[i];

		//		// check if file exists
		//		fileCheck(v, rmax, Ca, info);
		//		if (info){ // file exists
		//			// update solution
		//			readOutput(n, m, v, rmax, Ca, t.data(), si.data());
		//			for (j = 0; j < m*n; j++){
		//				sm[j] = si[j];
		//			}
		//		//	for (j = 0; j < m; j++){
		//		//		for (int jj = 0; jj < n; jj++){
		//		//			cout << sm[j*n + jj] << " ";
		//		//		}
		//		//		cout << endl;
		//		//	}
		//			
		//			// skip to next step
		//			cout << "Output file already exists. Skipping..." << endl;
		//			continue;
		//		}

				// multiple shooting method
				cout << "Shooting for v = 0." << v << ", rmax = 0." 
				     << rmax << ", Ca = " << Ca << "." << endl;
				mshoot(n, m, nrk, Ca, area, vlme, t.data(), si.data(), sf.data(), flag);
			
				// write to file
				if (flag == 0)
					writeSoln(n, m, v, rmax, Ca, t.data(), sf.data(), opath);
			
				/* update next initial guess using
				 * first-order continuation */
				if (flag == 0){
					if (i != vecCa.size() - 1){
						if (i == 0) {
							dC0 = vecCa[i  ] - 0;
							dC1 = vecCa[i+1] - 0;
						}
						else {
							dC0 = vecCa[i  ] - vecCa[i-1];
							dC1 = vecCa[i+1] - vecCa[i-1];
						}
						slope = dC1/dC0;
						for (j = 0; j < m*n; j++){
							//si[j] = sm[j] + slope*(sf[j] - sm[j]);
							si[j] = sm[j];
							sm[j] = sf[j];
						}
					}
				}
			}
		}
	}

	return(0);
}

void init(int n, int m, int v, int rmax, 
          double &area, double &vlme, double *t, double *s){
	// declare variables
	int    i, j;
	vector<double> T(m), S(m*n);
	
	// read file
	readEquil(n, m, v, rmax, T.data(), S.data());

	// copy abscissa and solution vectors
	for (i = 0; i < m; i++){
		t[i] = T[i];
		for (j = 0; j < n; j++){
			s[i*n + j] = S[i*n + j];
		}
	}

	// hard-code area and volume, for now
	area = 14.907;
	vlme = 5.103;
	getAreaVlme(v, rmax, area, vlme);

//	// parameters
//	area = s[(m-1)*n + 6];
//	vlme = s[(m-1)*n + 7];
}

void getAreaVlme(int v, int rmax, double &area, double &vlme){
	// lookup table (from BEM simulation results)
	if (v == 95){
		if      (rmax == 90){
			area = 14.907;
			vlme =  5.103;
		}
		else if (rmax == 92){
			area = 15.581;
			vlme =  5.449; 
		}
		else if (rmax == 94){
			area = 16.285;
			vlme =  5.822;
		}
		else if (rmax == 96){
			area = 17.023;
			vlme =  6.215;
		}
	}
	else if (v == 90){
		if      (rmax == 90){
			area = 18.825 ;
			vlme =  6.8555;
		}
		else if (rmax == 92){
			area = 19.690 ;
			vlme =  7.3285; 
		}
		else if (rmax == 94){
			area = 20.607;
			vlme =  7.841;
		}
		else if (rmax == 96){
			area = 21.693;
			vlme =  8.442;
		}
	}
	else if (v == 85){
		if      (rmax == 90){
			area = 22.615;
			vlme =  8.523;
		}
		else if (rmax == 92){
			area = 23.694;
			vlme =  9.128;
		}
		else if (rmax == 94){
			area = 24.9;
			vlme =  9.8;
		}
		else if (rmax == 96){
			area = 26.08 ;
			vlme = 10.494;
		}
	}
	else if (v == 80){
		if      (rmax == 90){
			area = 25.285;
			vlme =  9.487;
		}
		else if (rmax == 92){
			area = 26.46;
			vlme = 10.15;
		}
		else if (rmax == 94){
			area = 27.72;
			vlme = 10.87;
		}
		else if (rmax == 96){
			area = 29.28 ;
			vlme = 11.737;
		}
	}
	else if (v == 75){
		if      (rmax == 90){
			area = 27.55;
			vlme = 10.11;
		}
		else if (rmax == 92){
			area = 28.63;
			vlme = 10.72; 
		}
		else if (rmax == 94){
			area = 30.3 ;
			vlme = 11.61;
		}
		else if (rmax == 96){
			area = 31.9 ;
			vlme = 12.49;
		}
	}
	else{
		area = 0;
		vlme = 0;
	}
}
