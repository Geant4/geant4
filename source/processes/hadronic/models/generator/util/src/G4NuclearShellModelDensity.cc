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
// $Id: G4NuclearShellModelDensity.cc,v 1.8 2002/12/12 19:17:58 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//

#include "G4NuclearShellModelDensity.hh"

G4NuclearShellModelDensity::G4NuclearShellModelDensity(G4double anA, G4double aZ) 
{
        const G4double r0sq=0.8133*fermi*fermi;
	theA = G4int(anA);
	theZ = G4int(aZ);
	theRsquare= r0sq * pow(theA, 2./3. );
	Setrho0(pow(1./(pi*theRsquare),3./2.));
}

G4NuclearShellModelDensity::~G4NuclearShellModelDensity() {}
    
G4double G4NuclearShellModelDensity::GetRelativeDensity(const G4ThreeVector & aPosition) const
{
	return exp(-1*aPosition.mag2()/theRsquare);
}
    
G4double G4NuclearShellModelDensity::GetRadius(const G4double maxRelativeDensity) const
{

     return (maxRelativeDensity>0 && maxRelativeDensity <= 1 ) ?
             sqrt(theRsquare * log(1/maxRelativeDensity) ) : DBL_MAX;
}
   
G4double   G4NuclearShellModelDensity::GetDeriv(const G4ThreeVector & aPosition) const
{
     return -2* aPosition.mag() / theRsquare * GetDensity(aPosition);
}
