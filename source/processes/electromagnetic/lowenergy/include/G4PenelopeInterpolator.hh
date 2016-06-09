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
//
//
// Author:  Luciano Pandola (Luciano.Pandola@cern.ch)
//
// History:
// -----------
// 17 Feb 2003   LP        Created
// 17 Dec 2003   LP        Bug fixed (removed memory leak)
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

  ~G4PenelopeInterpolator();
 
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
 
