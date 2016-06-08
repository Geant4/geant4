// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NuclearFermiDensity.cc,v 1.1 1999/01/07 16:12:24 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#include "G4NuclearFermiDensity.hh"

G4NuclearFermiDensity::G4NuclearFermiDensity(G4double anA, G4double aZ) 
  :  a(0.545 * fermi) 
{
        const G4double r0=1.14*fermi;
	theA = anA;
	theZ = aZ;
	theR= r0 * pow(theA, 1./3. );
	Setrho0(3./4./pi/pow(r0,3.)/theA * ( 1. + sqr(a/theR)*pi2 ));
}

G4NuclearFermiDensity::~G4NuclearFermiDensity() {}
    
G4double G4NuclearFermiDensity::GetRelativeDensity(G4ThreeVector aPosition)
{
	return 1./(1.+exp((aPosition.mag()-theR)/a));
}
    
G4double G4NuclearFermiDensity::GetRadius(const G4double maxRelativeDensity)
{
     return (maxRelativeDensity>0 && maxRelativeDensity <= 1 ) ?
           (theR + a*log((1-maxRelativeDensity+exp(-1*theR/a))/maxRelativeDensity))  : DBL_MAX;
}
   
