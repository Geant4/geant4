// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LeptonConstructor.hh,v 1.2 1999-12-15 14:51:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
