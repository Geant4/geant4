// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN04EMPhysics.hh,v 1.1 2001-01-06 06:56:16 kurasige Exp $
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
#ifndef ExN04EMPhysics_h
#define ExN04EMPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

class ExN04EMPhysics : public G4VPhysicsConstructor
{
  public: 
    ExN04EMPhysics(const G4String& name ="EM");
    virtual ~ExN04EMPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  protected:
   // Gamma physics
    G4PhotoElectricEffect thePhotoEffect;
    G4ComptonScattering theComptonEffect;
    G4GammaConversion thePairProduction;
  
    // Electron physics
    G4MultipleScattering theElectronMultipleScattering;
    G4eIonisation theElectronIonisation;
    G4eBremsstrahlung theElectronBremsStrahlung;
  
    //Positron physics
    G4MultipleScattering thePositronMultipleScattering;
    G4eIonisation thePositronIonisation; 
    G4eBremsstrahlung thePositronBremsStrahlung;  
    G4eplusAnnihilation theAnnihilation;
};


#endif





