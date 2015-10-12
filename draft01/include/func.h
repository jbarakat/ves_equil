/* FUNCTION EVALUATION
 *  Evaluate right-hand side vector f in the system of ODEs dy/dt = f(y,t).
 *
 * REFERENCES
 *  Secomb et al, J Fluid Mech (1986)
 *  
 * PARAMETERS
 *  ny					[input]		number of ODEs
 *  t						[input]		abscissas
 *  y						[input]		solution
 *  f						[input]		function
 *  r   = y[0]  [input]		radius
 *  psi = y[1]  [input]		tilt angle
 *  cs  = y[2]  [input]		meridional curvature
 *  qs  = y[3]  [input]		transverse shear tension
 *  p   = y[4]  [input]		pressure
 *  sig = y[5]  [input]		mean tension
 *  A   = y[6]  [input]		total surface area
 *  V   = y[7]  [input]		total volume
 *  Q   = y[8]  [input]		leakback flow rate
 *  S   = y[9]  [input]		total meridional arc length
 */

#ifndef FUNC_H
#define FUNC_H


/* HEADER FILES */
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

using namespace std;

/* PROTOTYPES */

/* IMPLEMENTATIONS */
void func(int ny, double Ca,
          double t, double *y, double *f){
  int i;
  
  if (ny != 9)
    cout << "Error: ny should equal 9." << endl;
  
  // define variables
  double r   = y[0];
  double psi = y[1];
  double cs  = y[2];
  double qs  = y[3];
  double A   = y[4];
  double V   = y[5];
	double p   = y[6];
	double sig = y[7];
  double S   = y[8];
	
	// check bounds on radius
	if (r > 0.99999999)
		r = 0.99999999;
	
	if (r <= 0)
		r = 1e-12;
	
	if (r > 0.02){ // check if far from end caps
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);
  	double cosr = cos/r;
  	double sinr = sin/r;

  	// calculate function
  	f[0] = -sin;
  	f[1] = -cs;
  	f[2] = sinr*cosr
  	     + cs*sinr 
  	     - qs;
  	f[3] = p
  	     - sig*(cs - cosr)
  	     + 0.5*(cs + cosr)*(cs*cs - cosr*cosr)
  	     + qs*sinr;
  	f[4] = 2*M_PI*r;
  	f[5] = M_PI*r*r*cos;
  	f[6] = 0;
  	f[7] = 0;
  	f[8] = 0;
	}
	else { // approximate cs = cphi as constant (spherical cap)
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);
		
		f[0] = -sin;
		f[1] = -cs;
		f[2] = 0;
		f[3] = 0.5*(p - 2.0*cs*sig);
		f[4] = 2.0*M_PI*r;
		f[5] = M_PI*r*r*cos;
		f[6] = 0;
  	f[7] = 0;
  	f[8] = 0;
	}

	for (i = 0; i < ny; i++){
		f[i] *= S/2.0;
	}
}



#endif
