/* READ FILE
 *  Read .dat file containing the abscissa and solution vectors.
 *
 * REFERENCES
 *  Seifert et al, Phys Rev A (1991)
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

#ifndef READ_H
#define READ_H

/* HEADER FILES */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

/* PROTOTYPES */
void fileCheck(int, int, double, bool &);
void readEquil(int, int, int, int, double *, double *);
void readOutput(int, int, int, int, double, double *, double *);

/* IMPLEMENTATIONS */
///* Check if output file exists. */
//void fileCheck(int v, int rmax, double Ca, bool &info){
//	string filename, dir;
//	
//	// convert Ca double to string
//	ostringstream capnum;
//  if (Ca >= 0.0001 && Ca < 0.001)
//    capnum << fixed << setprecision(4) << Ca; 
//  if (Ca >= 0.001 && Ca < 0.01)
//    capnum << fixed << setprecision(3) << Ca; 
//  if (Ca >= 0.01 && Ca < 0.1)
//    capnum << fixed << setprecision(2) << Ca; 
//  if (Ca >= 0.1 && Ca < 1.0)
//    capnum << fixed << setprecision(1) << Ca; 
//  if (Ca >= 1.0)
//    capnum << fixed << setprecision(0) << Ca;
// 
//	// get filename
//	ostringstream redvol, maxrad;
//	redvol << fixed << v;
//	maxrad << fixed << rmax;
//	dir = "../output";
//	filename = "/v" + redvol.str() + "/rmax" + maxrad.str()
//	           + "/sln_v" + redvol.str() + "_rmax" + maxrad.str()
//	           + "_Ca" + capnum.str();
//	filename.erase(remove(filename.begin(), filename.end(), '.'), filename.end());
//	filename = dir + filename + ".dat";
//	
//	// check if file exists
//	ifstream file(filename.c_str());
//	if (file)
//		info = true; // file exists
//	else
//		info = false; // file does not exist
//}

/* Read .dat file containing the equilibrium solution. */
void readEquil(int n, int mh, int v, int rmax, double *th, double *sh){
	// error flags
	if (n != 9){
		cout << "Error: support only for n = 9." << endl;
		return;
	}

//	if (v != 65 && v != 70 && v != 75 && v!=80 && v != 85 && v != 90 && v != 95){
//		cout << "Error: support only for v = 85, 90, 95." << endl;
//		return;
//	}
//
//	if (rmax < 90 || rmax > 99){
//		cout << "Error: support only for 90 <= rmax <= 99." << endl;
//		return;
//	}

	typedef vector<double> vdouble;
	int m = 2*mh - 1; // guaranteed to be odd

	// declare variables
	int    i, j, k, i1, j1, k1;
	string filename, line;
	ifstream file;
	vdouble T;
	vdouble R   , PSI   , CS   , QS   , P   , SIG   , AA   , VV   , QQ   , SS   ;
	vdouble t(m);
	vdouble r(m), psi(m), cs(m), qs(m), p(m), sig(m), A (m), V (m), Q (m), S (m);
 
	// get filename
	ostringstream redvol, maxrad, capnum;
	redvol << fixed << v;
	maxrad << fixed << rmax;
	capnum << 0;
	filename = "../equil/pro_v" + redvol.str() + "_rmax" + maxrad.str()
	           + "_Ca" + capnum.str() + ".dat";
	cout << "Reading " << filename << "." << endl;

	// open file
	file.open(filename.c_str());
	while (file){
		getline(file, line);
		istringstream iss(line);
		vdouble row;
		
		for (i = 0; i < 11; i++){
			double data;
			iss >> data;
			row.push_back(data);
		}
		
		if (file.eof())
			break;

		T  .push_back(row[0 ]);
		R  .push_back(row[1 ]);
		PSI.push_back(row[2 ]);
		CS .push_back(row[3 ]);
		QS .push_back(row[4 ]);
		P  .push_back(row[5 ]);
		SIG.push_back(row[6 ]);
		AA .push_back(row[7 ]);
		VV .push_back(row[8 ]);
		QQ .push_back(row[9 ]);
		SS .push_back(row[10]);
	}

	file.close();

	// interpolate
	int    nt  = T.size() - 1;
	double t0 = 0.0;
	double t1 = 1.0;
	double dt = (t1 - t0)/(m-1);
	double ti;
	double DT, Dt;

	// assemble the abscissas
	for (i = 0; i < m; i++){
		t[i] = i*dt;
	}

	for (i = 0; i < m; i++){
		ti = t[i];

		// find ti in T
		for (j = 0; j < nt; j++){
			j1 = j + 1;
			if (ti >= T[j] && ti < T[j1] ){
				k = j;
				break;
			}
			else if (ti == T[n]){
				k = nt-1;
				break;
			}
		}

		k1 = k + 1;
		DT = T[k1] - T[k];
		Dt = ti    - T[k];

		// linear interpolation
		r  [i] = R  [k] + (R  [k1] - R  [k])*(Dt/DT);
		psi[i] = PSI[k] + (PSI[k1] - PSI[k])*(Dt/DT);
		cs [i] = CS [k] + (CS [k1] - CS [k])*(Dt/DT);
		qs [i] = QS [k] + (QS [k1] - QS [k])*(Dt/DT);
		p  [i] = P  [k] + (P  [k1] - P  [k])*(Dt/DT);
		sig[i] = SIG[k] + (SIG[k1] - SIG[k])*(Dt/DT);
		A  [i] = AA [k] + (AA [k1] - AA [k])*(Dt/DT);
		V  [i] = VV [k] + (VV [k1] - VV [k])*(Dt/DT);
		Q  [i] = QQ [k] + (QQ [k1] - QQ [k])*(Dt/DT);
		S  [i] = SS [k] + (SS [k1] - SS [k])*(Dt/DT);
	}

	/* assemble the solution vector
	 * NOTE: special routine for equilibrium case - use only
	 *       one quadrant (t = 0 to 0.5) for the solution */
	double tmid = t[mh-1];
	for (i = 0; i < m; i++){
		t[i] = t[i]/tmid;
	}

	for (i = 0; i < mh; i++){
		th[i]       = t[i];
    sh[i*n + 0] = r  [i];
    sh[i*n + 1] = psi[i];
    sh[i*n + 2] = cs [i];
    sh[i*n + 3] = qs [i];
    sh[i*n + 4] = A  [i];
    sh[i*n + 5] = V  [i];
    sh[i*n + 6] = p  [i];
    sh[i*n + 7] = sig[i];
    sh[i*n + 8] = S  [i];		
	}
}

///* Read output .dat file containing the solution s. */
//void readOutput(int n, int m, int v, int rmax, double Ca, double *t, double *s){
//	// error flags
//	if (n != 10){
//		cout << "Error: support only for n = 10." << endl;
//		return;
//	}
//
//	if (v != 65 && v != 70 && v != 75 && v!=80 && v != 85 && v != 90 && v != 95){
//		cout << "Error: support only for v = 85, 90, 95." << endl;
//		return;
//	}
//
//	if (rmax < 90 || rmax > 99){
//		cout << "Error: support only for 90 <= rmax <= 99." << endl;
//		return;
//	}
//
//	typedef vector<double> vdouble;
//
//	// declare variables
//	int    i, j, k, i1, j1, k1;
//	string filename, line, dir;
//	ifstream file;
//	vdouble T;
//	vdouble R   , PSI   , CS   , QS   , P   , SIG   , AA   , VV   , QQ   , SS   ;
//	vdouble r(m), psi(m), cs(m), qs(m), p(m), sig(m), A (m), V (m), Q (m), S (m);
//
//	// convert Ca double to string
//	ostringstream capnum;
//  if (Ca >= 0.0001 && Ca < 0.001)
//    capnum << fixed << setprecision(4) << Ca; 
//  if (Ca >= 0.001 && Ca < 0.01)
//    capnum << fixed << setprecision(3) << Ca; 
//  if (Ca >= 0.01 && Ca < 0.1)
//    capnum << fixed << setprecision(2) << Ca; 
//  if (Ca >= 0.1 && Ca < 1.0)
//    capnum << fixed << setprecision(1) << Ca; 
//  if (Ca >= 1.0)
//    capnum << fixed << setprecision(0) << Ca;
// 
//	// get filename
//	ostringstream redvol, maxrad;
//	redvol << fixed << v;
//	maxrad << fixed << rmax;
//	dir = "../output";
//	filename = "/v" + redvol.str() + "/rmax" + maxrad.str()
//	           + "/sln_v" + redvol.str() + "_rmax" + maxrad.str()
//	           + "_Ca" + capnum.str();
//	filename.erase(remove(filename.begin(), filename.end(), '.'), filename.end());
//	filename = dir + filename + ".dat";
//	cout << "Reading " << filename << "." << endl;
//
//	// open file
//	file.open(filename.c_str());
//	while (file){
//		getline(file, line);
//		istringstream iss(line);
//		vdouble row;
//		
//		for (i = 0; i < 11; i++){
//			double data;
//			iss >> data;
//			row.push_back(data);
//		}
//		
//		if (file.eof())
//			break;
//
//		T  .push_back(row[0 ]);
//		R  .push_back(row[1 ]);
//		PSI.push_back(row[2 ]);
//		CS .push_back(row[3 ]);
//		QS .push_back(row[4 ]);
//		P  .push_back(row[5 ]);
//		SIG.push_back(row[6 ]);
//		AA .push_back(row[7 ]);
//		VV .push_back(row[8 ]);
//		QQ .push_back(row[9 ]);
//		SS .push_back(row[10]);
//	}
//
//	file.close();
//
//	// interpolate
//	int    nt  = T.size() - 1;
//	double t0 = 0.0;
//	double t1 = 1.0;
//	double dt = (t1 - t0)/(m-1);
//	double ti;
//	double DT, Dt;
//
//	// assemble the abscissas
//	for (i = 0; i < m; i++){
//		t[i] = i*dt;
//	}
//
//	for (i = 0; i < m; i++){
//		ti = t[i];
//
//		// find ti in T
//		for (j = 0; j < nt; j++){
//			j1 = j + 1;
//			if (ti >= T[j] && ti < T[j1] ){
//				k = j;
//				break;
//			}
//			else if (ti == T[n]){
//				k = nt-1;
//				break;
//			}
//		}
//
//		k1 = k + 1;
//		DT = T[k1] - T[k];
//		Dt = ti    - T[k];
//
//		// linear interpolation
//		r  [i] = R  [k] + (R  [k1] - R  [k])*(Dt/DT);
//		psi[i] = PSI[k] + (PSI[k1] - PSI[k])*(Dt/DT);
//		cs [i] = CS [k] + (CS [k1] - CS [k])*(Dt/DT);
//		qs [i] = QS [k] + (QS [k1] - QS [k])*(Dt/DT);
//		p  [i] = P  [k] + (P  [k1] - P  [k])*(Dt/DT);
//		sig[i] = SIG[k] + (SIG[k1] - SIG[k])*(Dt/DT);
//		A  [i] = AA [k] + (AA [k1] - AA [k])*(Dt/DT);
//		V  [i] = VV [k] + (VV [k1] - VV [k])*(Dt/DT);
//		Q  [i] = QQ [k] + (QQ [k1] - QQ [k])*(Dt/DT);
//		S  [i] = SS [k] + (SS [k1] - SS [k])*(Dt/DT);
//	}
//
//	// assemble the solution vector
//	for (i = 0; i < m; i++){
//    s[i*n + 0] = r  [i];
//    s[i*n + 1] = psi[i];
//    s[i*n + 2] = cs [i];
//    s[i*n + 3] = qs [i];
//    s[i*n + 4] = p  [i];
//    s[i*n + 5] = sig[i];
//    s[i*n + 6] = A  [i];
//    s[i*n + 7] = V  [i];
//    s[i*n + 8] = Q  [i];
//    s[i*n + 9] = S  [i];		
//	}
//}

#endif
