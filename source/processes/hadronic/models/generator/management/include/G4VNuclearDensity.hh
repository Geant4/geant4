// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VNuclearDensity.hh,v 1.3 2000/05/23 13:41:56 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
#ifndef G4VNuclearDensity_h
#define G4VNuclearDensity_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"


class G4VNuclearDensity 
{

  public:
    G4VNuclearDensity();
    virtual ~G4VNuclearDensity();
    
    inline G4double GetDensity(G4ThreeVector aPosition)
    {
	return rho0*GetRelativeDensity(aPosition);
    };
    
    virtual G4double GetRelativeDensity(G4ThreeVector aPosition) = 0;
    virtual G4double GetRadius(const G4double maxRelativeDenisty) = 0;
    virtual G4double GetDeriv(const G4ThreeVector & point) = 0;    

  protected:    
    inline void Setrho0(G4double arho0) { rho0=arho0; };
    inline G4double Getrho0() { return rho0; };
   
  private:
  
    G4double rho0;
};

#endif

