//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef MYGAMMA_H
#define MYGAMMA_H

#include <cmath>

class MyGamma {
	public:
	MyGamma ();
	~MyGamma();
	
	double Gamma(double z);
	double  Gamma(double a,double x);
	private:
	double GamCf(double a,double x);
	double GamSer(double a,double x);
	
	// Abs
	static short  Abs(short d)  { return (d > 0) ? d : -d; }
	static int    Abs(int d)   { return (d > 0) ? d : -d; }
	static long   Abs(long d)  { return (d > 0) ? d : -d; }
	static float  Abs(float d){ return (d > 0) ? d : -d; }
	static double Abs(double d)   { return (d > 0) ? d : -d; }
	static double LnGamma(double z); 
	static double Log(double x){ return std::log(x); }
	static double Exp(double x){ return std::exp(x); }
	
	
};

#endif
