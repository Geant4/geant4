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
#ifndef G4NuclearShellModelDensity_h
#define G4NuclearShellModelDensity_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VNuclearDensity.hh"

#include <CLHEP/Units/PhysicalConstants.h>	// pi, fermi,..
//#include <cmath>				// pow,..

class G4NuclearShellModelDensity : public G4VNuclearDensity
{

  public:
    G4NuclearShellModelDensity(G4int anA, G4int aZ);
    ~G4NuclearShellModelDensity();
    
    G4double GetRelativeDensity(const G4ThreeVector & aPosition) const;
    G4double GetRadius(const G4double maxRelativeDenisty) const;
    G4double GetDeriv(const G4ThreeVector & aPosition) const;    

  private:
    G4int    theA;
    //G4int    theZ;
    G4double theRsquare;
  
};

#endif

