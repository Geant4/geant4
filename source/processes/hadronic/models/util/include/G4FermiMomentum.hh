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
#ifndef G4FermiMomentum_h
#define G4FermiMomentum_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4Pow.hh"

class G4FermiMomentum 
{

  public:
    G4FermiMomentum();
    ~G4FermiMomentum();
    
    inline void Init(G4int anA, G4int aZ) {theA = anA; theZ = aZ;}
    
    inline G4double GetFermiMomentum(G4double density)
    {
	return constofpmax * cbrt(density * theA);
    }
    
    inline G4ThreeVector GetMomentum(G4double density, 
    				     G4double maxMomentum=-1.)
    { 
	if (maxMomentum < 0 ) maxMomentum=GetFermiMomentum(density);
	G4ThreeVector p;
	
	do {
	    p=G4ThreeVector(2.*G4UniformRand()-1.,
	    		    2.*G4UniformRand()-1.,
	    		    2.*G4UniformRand()-1.);
	    } while ( p.mag() > 1. );  /* Loop checking, 30-Oct-2015, G.Folger */
	return p*maxMomentum; 
     }

  private:

  G4double cbrt(G4double x) { return G4Pow::GetInstance()->A13(x); }

  private:
  
    G4int theA;
    G4int theZ;
      // pmax= hbar * c * ( 3* pi**2 * rho )**(1/3) =
      //       hbar * c * ( 3* pi**2 )**(1/3) * rho**(1/3)=
      //       constofpmax * rho**(1/3)
    G4double constofpmax;
  
};

#endif

