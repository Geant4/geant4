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

#ifndef G4NuclearFermiDensity_h
#define G4NuclearFermiDensity_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VNuclearDensity.hh"

#include <CLHEP/Units/PhysicalConstants.h>	// pi, fermi,..
#include "G4Exp.hh"
#include "G4Log.hh"
//#include <cmath>				// pow

class G4NuclearFermiDensity : public G4VNuclearDensity
{
  public:
    G4NuclearFermiDensity(G4int anA, G4int aZ);
    ~G4NuclearFermiDensity();
    
    G4double GetRelativeDensity(const G4ThreeVector & aPosition) const
    {
      return 1./(1.+G4Exp((aPosition.mag()-theR)/a));
    }
    
    G4double GetRadius(const G4double maxRelativeDenisty) const
    {
      return (maxRelativeDenisty>0 && maxRelativeDenisty <= 1 ) ?
             (theR + a*G4Log((1-maxRelativeDenisty+G4Exp(-1*theR/a))/maxRelativeDenisty))  : DBL_MAX;
    }
    
    G4double GetDeriv(const G4ThreeVector & aPosition) const
    {
      G4double currentR=aPosition.mag();
      if (currentR > 40*theR  ) {return 0;}
      else return -G4Exp((currentR-theR)/a) * sqr(GetDensity(aPosition)) / (a*Getrho0());
    }   
   
  private:
    G4int theA;
    G4double theR;      // Nuclear Radius 
    const G4double a;	// Determines the nuclear surface thickness
};

#endif

