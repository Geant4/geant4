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
// $Id: G4FermiMomentum.hh,v 1.6 2002/12/12 19:17:57 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
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

    G4double cbrt(G4double x) { return pow(x,1./3.); }

  private:
  
    G4double theA;
    G4double theZ;
      // pmax= hbar * c * ( 3* pi**2 * rho )**(1/3) =
      //       hbar * c * ( 3* pi**2 )**(1/3) * rho**(1/3)=
      //       constofpmax * rho**(1/3)
    G4double constofpmax;
  
};

#endif

