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
#include <math.h>

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
 
