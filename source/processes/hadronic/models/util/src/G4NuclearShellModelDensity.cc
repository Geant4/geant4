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

#include "G4NuclearShellModelDensity.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

G4NuclearShellModelDensity::G4NuclearShellModelDensity(G4int anA, G4int /*aZ*/) 
: theA(anA)//, theZ(aZ)
{
        const G4double r0sq=0.8133*fermi*fermi;
	theRsquare= r0sq * G4Pow::GetInstance()->Z23(theA);
	G4double x = 1./(pi*theRsquare);
	Setrho0(x*std::sqrt(x));
}

G4NuclearShellModelDensity::~G4NuclearShellModelDensity() {}
    
G4double G4NuclearShellModelDensity::GetRelativeDensity(const G4ThreeVector & aPosition) const
{
	return G4Exp(-1*aPosition.mag2()/theRsquare);
}
    
G4double G4NuclearShellModelDensity::GetRadius(const G4double maxRelativeDensity) const
{

     return (maxRelativeDensity>0 && maxRelativeDensity <= 1 ) ?
             std::sqrt(theRsquare * G4Log(1/maxRelativeDensity) ) : DBL_MAX;
}
   
G4double   G4NuclearShellModelDensity::GetDeriv(const G4ThreeVector & aPosition) const
{
     return -2* aPosition.mag() / theRsquare * GetDensity(aPosition);
}
