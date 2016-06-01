// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LeptonConstructor.hh,v 2.1 1998/07/13 17:17:25 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementatLepton file 
//
#ifndef G4LeptonConstructor_h
#define G4LeptonConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4LeptonConstructor
{
  //This class is a utility class for constructLepton 
  //short lived particles

  public:
    G4LeptonConstructor();
    ~G4LeptonConstructor();
  
  public:
    void ConstructParticle();

  protected:
    void ConstructELeptons();
    void ConstructMuLeptons();
    void ConstructTauLeptons();
};

#endif
