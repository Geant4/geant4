// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BosonConstructor.hh,v 1.2 1999-12-15 14:50:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementatBoson file 
//
#ifndef G4BosonConstructor_h
#define G4BosonConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4BosonConstructor
{
  //This class is a utility class for constructBoson 
  //short lived particles

  public:
    G4BosonConstructor();
    ~G4BosonConstructor();
  
  public:
    void ConstructParticle();
};

#endif
