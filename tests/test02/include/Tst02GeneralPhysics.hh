// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst02GeneralPhysics.hh,v 1.1 2001-01-06 07:04:24 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
// Class Description:
//      This class is an derived class of G4VPhysicsConstructor
//
// ------------------------------------------- 
//	History
//        first version                   12 Nov. 2000 by H.Kurashige 
// ------------------------------------------------------------
#ifndef Tst02GeneralPhysics_h
#define Tst02GeneralPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"


#include "G4Decay.hh"

class Tst02GeneralPhysics : public G4VPhysicsConstructor
{
  public: 
    Tst02GeneralPhysics(const G4String& name = "general");
    virtual ~Tst02GeneralPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  protected:
    G4Decay fDecayProcess;
};


#endif








