// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst02MuonPhysics.hh,v 1.2 2001-03-06 10:38:35 kurasige Exp $
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
#ifndef Tst02MuonPhysics_h
#define Tst02MuonPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MultipleScattering.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuIonisation.hh"
#include "G4hIonisation.hh"
#include "G4MuonMinusCaptureAtRest.hh"

class Tst02MuonPhysics : public G4VPhysicsConstructor
{
  public: 
    Tst02MuonPhysics(const G4String& name="muon");
    virtual ~Tst02MuonPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  protected:
   // Muon physics
   G4MultipleScattering   fMuPlusMultipleScattering;
   G4MuBremsstrahlung     fMuPlusBremsstrahlung ;
   G4MuPairProduction     fMuPlusPairProduction;
   G4MuIonisation         fMuPlusIonisation;

   G4MultipleScattering   fMuMinusMultipleScattering;
   G4MuBremsstrahlung     fMuMinusBremsstrahlung ;
   G4MuPairProduction     fMuMinusPairProduction;
   G4MuIonisation         fMuMinusIonisation;
   G4MuonMinusCaptureAtRest fMuonMinusCaptureAtRest;
   // Tau physics
   G4MultipleScattering   fTauPlusMultipleScattering;
   G4hIonisation          fTauPlusIonisation;

   G4MultipleScattering   fTauMinusMultipleScattering;
   G4hIonisation          fTauMinusIonisation;

};


#endif

