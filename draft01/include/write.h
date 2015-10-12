/* WRITE FILE
 *  Write .dat file containing the abscissa and solution vectors.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 *  n			[input]		number of ODEs
 *  m			[input]		number of shooting points
 *  v			[input]		reduced volume (x 100)
 *  rmax	[input]		maximum radius (x 100)
 *  t			[output]	abscissas (meridional arc length scaled to the [0,1] axis)
 *  r			[output]	cylindrical radius
 *  psi		[output]	tilt angle
 *  cs		[output]	meridional curvature (points outward from closed contour)
 *  qs		[output]	meridional component of transverse shear tension
 *  p			[output]	pressure difference between exterior and interior
 *  sig		[output]	mean tension
 *  A			[output]	surface area
 *  V			[output]	volume
 *	Q			[output]	leakback flow rate
 *  S			[output]	total meridional arc length (half-space)
 */

#ifndef WRITE_H
#define WRITE_H

/* HEADER FILES */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

/* PROTOTYPES */

/* IMPLEMENTATIONS */

/* Write .dat file. */
void writeSoln(int n, int mh, int v, int rmax, double Ca, 
               double *th, double *sh, string path){
	// error flag
	if (n != 9){
		cout << "Error: support only for n = 9." << endl;
		return;
	}

	// declare variables
	int i, j;
	int width = 12;
	int m = 2*mh - 1;

	ostringstream ssRedVol, ssMaxRad, ssCapNum;
	ostringstream header;
	ofstream file;
	string dir, fn, line;

	vector<double> t(m,0.0), r(m,0.0), psi(m,0.0), cs(m,0.0), qs(m,0.0), p(m,0.0), sig(m,0.0), A(m,0.0), V(m,0.0), Q(m,0.0), S(m,0.0);
	
	double area = 2.0*sh[(mh-1)*n + 4];
	double vlme = 2.0*sh[(mh-1)*n + 5];

	// convert integers to strings
	ssRedVol << v   ;
	ssMaxRad << rmax;

	// convert double to string
	ssCapNum << fixed << setprecision(0) << Ca;
	
	// set directory
//	dir = path + "/v"    + ssRedVol.str() 
//	           + "/rmax" + ssMaxRad.str();
	dir = path;
	
	// set filename
	fn = "/pro_v" + ssRedVol.str()
	   + "_rmax"  + ssMaxRad.str() 
		 + "_Ca0"   ;
	fn.erase(remove(fn.begin(), fn.end(), '.'), fn.end());
	fn = dir + fn + ".dat";

	// reflect solution about symmetry axis
	for (i = 0; i < mh; i++){
		t  [i] = 0.5*th[i]      ; t  [m-1-i] =  1.0 - 0.5*th[i];
		r  [i] = sh[i*n + 0]; r  [m-1-i] =  sh[i*n + 0]; // even function
		psi[i] = sh[i*n + 1]; psi[m-1-i] = -sh[i*n + 1]; // odd function
		cs [i] = sh[i*n + 2]; cs [m-1-i] =  sh[i*n + 2]; // even function
		qs [i] = sh[i*n + 3]; qs [m-1-i] = -sh[i*n + 3]; // odd function
		p  [i] = sh[i*n + 6]; p  [m-1-i] =  sh[i*n + 6]; // constant function
		sig[i] = sh[i*n + 7]; sig[m-1-i] =  sh[i*n + 7]; // constant function
		A  [i] = sh[i*n + 4]; A  [m-1-i] = -sh[i*n + 4] + area; // odd function
		V  [i] = sh[i*n + 5]; V  [m-1-i] = -sh[i*n + 5] + vlme; // odd function
		Q  [i] = 0.0        ; Q  [m-1-i] =  0.0        ; // constant function
		S  [i] = sh[i*n + 8]; S  [m-1-i] =  sh[i*n + 8]; // constant function
	}
//	for (i = 0; i < m; i++){
//		cout << t[i] << "    " << qs[i] << endl;
//	}

	// write to file
	file.open(fn.c_str());
	header << setw(width) << "t" 
	       << setw(width) << "r" 
				 << setw(width) << "psi"
				 << setw(width) << "cs"
				 << setw(width) << "qs"
	       << setw(width+4) << "p" 
				 << setw(width+4) << "sig"
				 << setw(width) << "A"
				 << setw(width) << "V"
				 << setw(width) << "Q"
				 << setw(width) << "S";
	file << header.str() << endl;
	for (i = 0; i < m; i++){
		ostringstream tt, rr, ppsi, ccs, qqs, pp, 
		                  ssig, AA, VV , QQ , SS;
		tt   << fixed << setw(width)   << t  [i];
		rr   << fixed << setw(width)   << r  [i];
		ppsi << fixed << setw(width)   << psi[i];
		ccs  << fixed << setw(width)   << cs [i];
		qqs  << fixed << setw(width)   << qs [i];
		pp   << fixed << setw(width+4) << p  [i];
		ssig << fixed << setw(width+4) << sig[i];
		AA   << fixed << setw(width)   << A  [i];
		VV   << fixed << setw(width)   << V  [i];
		QQ   << fixed << setw(width)   << Q  [i];
		SS   << fixed << setw(width)   << S  [i];
		
		line = tt.str() + rr.str() + ppsi.str() + ccs.str() + qqs.str()
		     + pp.str() + ssig.str() + AA.str() + VV.str() + QQ.str() + SS.str();

		file << line << endl;
	}
	file.close();

	cout << "Solution written to " << fn << "." << endl;
}


#endif
