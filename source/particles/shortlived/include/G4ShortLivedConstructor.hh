// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ShortLivedConstructor.hh,v 1.2 1999-10-03 09:11:47 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
#ifndef G4ShortLivedConstructor_h
#define G4ShortLivedConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4ShortLivedConstructor
{
  //This class is a utility class for construction 
  //short lived particles

  public:
    G4ShortLivedConstructor();
    ~G4ShortLivedConstructor();
  
  public:
    void ConstructParticle();
 
  protected:
    void ConstructResonances();
    void ConstructBaryons();
    void ConstructMesons();
    void ConstructQuarks();

  private:
    static G4bool isConstructed;
    // flag for checking whether resonces exist or not
};

#endif











