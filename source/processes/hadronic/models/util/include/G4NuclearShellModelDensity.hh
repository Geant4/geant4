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
#ifndef G4NuclearShellModelDensity_h
#define G4NuclearShellModelDensity_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VNuclearDensity.hh"

#include <CLHEP/Units/PhysicalConstants.h>	// pi, fermi,..
#include <cmath>				// pow,..

class G4NuclearShellModelDensity : public G4VNuclearDensity
{

  public:
    G4NuclearShellModelDensity(G4double anA, G4double aZ);
    ~G4NuclearShellModelDensity();
    
    G4double GetRelativeDensity(const G4ThreeVector & aPosition) const;
    G4double GetRadius(const G4double maxRelativeDenisty) const;
    G4double GetDeriv(const G4ThreeVector & aPosition) const;    

  private:
    G4int    theA;
    G4int    theZ;
    G4double theRsquare;
  
};

#endif

