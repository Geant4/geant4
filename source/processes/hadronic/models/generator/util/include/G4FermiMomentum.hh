// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FermiMomentum.hh,v 1.1.10.1 1999/12/07 20:51:57 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
#ifndef G4FermiMomentum_h
#define G4FermiMomentum_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

class G4FermiMomentum 
{

  public:
    G4FermiMomentum();
    ~G4FermiMomentum();
    
    inline void Init(G4double anA, G4double aZ) {theA = anA; theZ = aZ;}
    
    inline G4double GetFermiMomentum(G4double density)
    {
	return constofpmax * cbrt(density * theA);
    }
    
    inline G4ThreeVector GetMomentum(G4double density)
    { 
	G4ThreeVector p;
	
	do {
	    p=G4ThreeVector(2.*G4UniformRand()-1.,
	    		    2.*G4UniformRand()-1.,
	    		    2.*G4UniformRand()-1.);
	    } while ( p.mag() > 1. );
	return p*GetFermiMomentum(density); 
     }
    
  private:
  
    G4double theA;
    G4double theZ;
      // pmax= hbar * c * ( 3* pi**2 * rho )**(1/3) =
      //       hbar * c * ( 3* pi**2 )**(1/3) * rho**(1/3)=
      //       constofpmax * rho**(1/3)
    G4double constofpmax;
  
};

#endif

