//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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


