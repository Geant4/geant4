// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NuclearFermiDensity.hh,v 1.1 1999/01/07 16:12:21 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
#ifndef G4NuclearFermiDensity_h
#define G4NuclearFermiDensity_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VNuclearDensity.hh"

#include <CLHEP/Units/PhysicalConstants.h>	// pi, fermi,..
#include <math.h>				// pow

class G4NuclearFermiDensity : public G4VNuclearDensity
{

  public:
    G4NuclearFermiDensity(G4double anA, G4double aZ);
    ~G4NuclearFermiDensity();
    
    G4double GetRelativeDensity(G4ThreeVector aPosition);
    G4double GetRadius(const G4double maxRelativeDenisty);
   
  private:
  
    G4int    theA;
    G4int    theZ;
    G4double theR;      // Nuclear Radius 
    const G4double a;	// Determines the nuclear surface thickness
  
};

#endif

