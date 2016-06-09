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
//
//
// Author:  Luciano Pandola (Luciano.Pandola@cern.ch)
//
// History:
// -----------
// 17 Feb 2003   LP        Created
//
// -------------------------------------------------------------------
// Class description:
// Cubic spline interpolation and momentum calculation of given order
// obtained by linear interpolation 
// -------------------------------------------------------------------

#ifndef G4PENELOPEINTERPOLATOR_HH
#define G4PENELOPEINTERPOLATOR_HH 1

#include "globals.hh"

class G4DataVector;
class G4PenelopeInterpolator
{
public:

  G4PenelopeInterpolator(G4double* pX, G4double* pY, G4int nbOfData, G4double S1=0, G4double SN=0);

  ~G4PenelopeInterpolator(){;}
 
  G4double CubicSplineInterpolation(G4double xx);
  G4double FirstDerivative(G4double xx);
  G4double CalculateMomentum(G4double UpperLimit,G4int MomentumOrder);
  
private:
  G4int FindBin(G4double xx);
  //coefficients
  G4DataVector* a;
  G4DataVector* b;
  G4DataVector* c;
  G4DataVector* d;
  //stored data
  G4DataVector* x;
  G4DataVector* y;
};
 
#endif
 
