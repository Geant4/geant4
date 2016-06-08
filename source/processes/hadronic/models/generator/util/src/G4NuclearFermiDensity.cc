// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NuclearFermiDensity.cc,v 1.6 2000/05/23 13:41:42 gunter Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

#include "G4NuclearFermiDensity.hh"

G4NuclearFermiDensity::G4NuclearFermiDensity(G4double anA, G4double aZ) 
  :  a(0.545 * fermi) 
{
//        const G4double r0=1.14*fermi;
	const G4double r0=1.16 * ( 1. - 1.16 * pow(anA, -2./3.)) * fermi;
	theA = G4int(anA);
	theZ = G4int(aZ);
	theR= r0 * pow(theA, 1./3. );
	Setrho0(3./ (4. * pi * pow(r0,3.) * theA * ( 1. + sqr(a/theR)*pi2 )));
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

G4double   G4NuclearFermiDensity::GetDeriv(const G4ThreeVector & aPosition)
{
     return -exp((aPosition.mag()-theR)/a) * sqr(GetDensity(aPosition)) / (a*Getrho0());
}
