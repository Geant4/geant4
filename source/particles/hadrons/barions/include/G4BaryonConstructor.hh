// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BaryonConstructor.hh,v 1.1.6.1 1999/12/07 20:49:25 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
#ifndef G4BaryonConstructor_h
#define G4BaryonConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4BaryonConstructor
{
  //This class is a utility class for construction 
  //short lived particles

  public:
    G4BaryonConstructor();
    ~G4BaryonConstructor();
  
  public:
    void ConstructParticle();

  protected:
    void ConstructNucleons();
    void ConstructStrangeBaryons();
    void ConstructCharmBaryons();
    void ConstructBottomBaryons();
};

#endif
