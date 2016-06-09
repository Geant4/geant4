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
// 30 Nov 2002   LP        Created
//
// -------------------------------------------------------------------
// Class description:
// Numerical quadrature of a function
// G4double function (G4double x)
// using Gauss adaptive-bipartition formula. 
// -------------------------------------------------------------------

#ifndef G4PENELOPEINTEGRATOR_HH
#define G4PENELOPEINTEGRATOR_HH 1

#include "globals.hh"
#include <cmath>

template <class T,class F> 
class G4PenelopeIntegrator 
{
 public:

  G4PenelopeIntegrator(){;}

  ~G4PenelopeIntegrator(){;}
 
  G4double Calculate(T& typeT,F f,G4double LowPoint,G4double HighPoint,
		     G4double MaxError);
  G4double Calculate(T* ptrT,F f,G4double LowPoint,G4double HighPoint,
		     G4double MaxError);

};

#include "G4PenelopeIntegrator.icc"
 
#endif
 
