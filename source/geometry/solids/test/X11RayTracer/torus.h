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
/*
  Torus header file
 */

typedef enum Localisation {
		kOutside,
		kSurface,
		kInside
} Localisation;


typedef struct Intersect {
		double x,y,z,dx,dy,dz;
		double R0,R1;
		double phi,deltaphi;
} Intersect ;

void PrintIteration (void);
/**/inline double TorusEquation (double x,double y,double z,double R0,double R1);
inline Localisation Inside (double x,double y,double z,double R0,double R1);
inline double TorusDerivativeX (double x,double y,double z,double R0,double R1);
inline double TorusDerivativeY (double x,double y,double z,double R0,double R1);
inline double TorusDerivativeZ (double x,double y,double z,double R0,double R1);
inline double TorusGuess (double x,double y,double z,double R0,double R1);
/*
double DistanceToTorus (double x,double y,double z,double dx,double dy,double dz
			,double R0,double R1);
*/
double DistanceToTorus (Intersect * Inter);
double Newton (double guess, double Lmin, double Lmax,
			   double x, double y, double z,
			   double dx, double dy, double dz,
			   double Rmax, double Rmin);

extern int BVM_ONLY ;


