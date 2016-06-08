// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NuclearShellModelDensity.cc,v 1.1 1999/01/07 16:12:24 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#include "G4NuclearShellModelDensity.hh"

G4NuclearShellModelDensity::G4NuclearShellModelDensity(G4double anA, G4double aZ) 
{
        const G4double r0sq=0.8133*fermi*fermi;
	theA = anA;
	theZ = aZ;
	theRsquare= r0sq * pow(theA, 2./3. );
	Setrho0(pow(1./(pi*theRsquare),3./2.));
}

G4NuclearShellModelDensity::~G4NuclearShellModelDensity() {}
    
G4double G4NuclearShellModelDensity::GetRelativeDensity(G4ThreeVector aPosition)
{
	return exp(-1*aPosition.mag2()/theRsquare);
}
    
G4double G4NuclearShellModelDensity::GetRadius(const G4double maxRelativeDensity)
{

     return (maxRelativeDensity>0 && maxRelativeDensity <= 1 ) ?
             sqrt(theRsquare * log(1/maxRelativeDensity) ) : DBL_MAX;
}
   
