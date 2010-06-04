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
// $Id: PhysicsListLHad.cc,v 1.4 2010-06-04 05:43:47 kmura Exp $
// ====================================================================
//   PhysicsListLHad.cc
//
//   Light package of hadron physics
// ====================================================================
#include "PhysicsListLHad.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"

// Hadronic Processes
#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"

// Low energy models
#include "G4LElastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"

// High-energy Models
#include "G4HEProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"

// Stopping processes
#include "G4AntiProtonAnnihilationAtRest.hh"

// Binary Cascade
#include "G4BinaryCascade.hh"
#include "G4ProtonInelasticCrossSection.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////
PhysicsListLHad::PhysicsListLHad()
  :  G4VPhysicsConstructor("LHad")
//////////////////////////////////
{
}


///////////////////////////////////
PhysicsListLHad::~PhysicsListLHad()
///////////////////////////////////
{
}


/////////////////////////////////////////
void PhysicsListLHad::ConstructParticle()
/////////////////////////////////////////
{
}


////////////////////////////////////////
void PhysicsListLHad::ConstructProcess()
////////////////////////////////////////
{
  G4ProcessManager* pManager;

  // ---------------------------------------------------------------
  // proton
  // ---------------------------------------------------------------
  pManager= G4Proton::Proton()-> GetProcessManager();

  // elastic
  G4HadronElasticProcess* thepElasticProcess = new G4HadronElasticProcess();
  G4LElastic* thepElasticModel = new G4LElastic();
  thepElasticProcess->RegisterMe(thepElasticModel);
  pManager-> AddDiscreteProcess(thepElasticProcess);

  // inelastic
  G4ProtonInelasticProcess* theProtonInelasticProcess 
    = new G4ProtonInelasticProcess(); 

  G4HEProtonInelastic* theProtonHEPModel = new G4HEProtonInelastic();  

  G4LEProtonInelastic* theProtonLEPModel = new G4LEProtonInelastic();
  theProtonLEPModel->SetMinEnergy(2.8*GeV);
   
  G4BinaryCascade* theProtonBICModel = new G4BinaryCascade();
  theProtonBICModel->SetMaxEnergy(3.2*GeV);

  theProtonInelasticProcess-> RegisterMe(theProtonHEPModel);
  theProtonInelasticProcess-> RegisterMe(theProtonLEPModel);
  theProtonInelasticProcess-> RegisterMe(theProtonBICModel);

  // add Xsection data of BIC 
  G4ProtonInelasticCrossSection* theProtonInelasticData 
    = new G4ProtonInelasticCrossSection();
  theProtonInelasticProcess-> AddDataSet( theProtonInelasticData );

  pManager-> AddDiscreteProcess(theProtonInelasticProcess);

  // QED
  G4VProcess* thepMultipleScattering = new G4hMultipleScattering();
  G4VProcess* thepIonisation         = new G4hIonisation();

  pManager-> AddProcess(thepIonisation);
  pManager-> AddProcess(thepMultipleScattering);

  pManager-> SetProcessOrdering(thepMultipleScattering, idxAlongStep, 1);
  pManager-> SetProcessOrdering(thepIonisation,         idxAlongStep, 2);

  pManager-> SetProcessOrdering(thepMultipleScattering, idxPostStep, 1);
  pManager-> SetProcessOrdering(thepIonisation,         idxPostStep, 2);

  // ---------------------------------------------------------------
  // anti-proton
  // ---------------------------------------------------------------
  pManager= G4AntiProton::AntiProton()-> GetProcessManager();

  // elastic
  G4HadronElasticProcess* theapElasticProcess = new G4HadronElasticProcess();
  G4LElastic* theapElasticModel = new G4LElastic();
  theapElasticProcess-> RegisterMe(theapElasticModel);
  pManager-> AddDiscreteProcess(theapElasticProcess);

  // inelastic
  G4AntiProtonInelasticProcess* theAntiProtonInelasticProcess 
    = new G4AntiProtonInelasticProcess(); 

  G4LEAntiProtonInelastic* theAntiProtonLEPModel 
    = new G4LEAntiProtonInelastic();
  G4HEAntiProtonInelastic* theAntiProtonHEPModel 
    = new G4HEAntiProtonInelastic();

  theAntiProtonInelasticProcess-> RegisterMe(theAntiProtonLEPModel);
  theAntiProtonInelasticProcess-> RegisterMe(theAntiProtonHEPModel);
  pManager-> AddDiscreteProcess(theAntiProtonInelasticProcess);

  G4AntiProtonAnnihilationAtRest* theAntiProtonAnnihilation
    =  new G4AntiProtonAnnihilationAtRest();
  pManager-> AddRestProcess(theAntiProtonAnnihilation);

  // QED
  G4VProcess* theapMultipleScattering = new G4hMultipleScattering();
  G4VProcess* theapIonisation         = new G4hIonisation();

  pManager-> AddProcess(theapIonisation);
  pManager-> AddProcess(theapMultipleScattering);

  pManager->SetProcessOrdering(theapMultipleScattering, idxAlongStep, 1);
  pManager->SetProcessOrdering(theapIonisation,         idxAlongStep, 2);

  pManager->SetProcessOrdering(theapMultipleScattering, idxPostStep, 1);
  pManager->SetProcessOrdering(theapIonisation,         idxPostStep, 2);


  // ---------------------------------------------------------------
  // mesons...
  // ---------------------------------------------------------------

  // ---------------------------------------------------------------
  // pi+
  // ---------------------------------------------------------------
  pManager= G4PionPlus::PionPlus()-> GetProcessManager();

  // elastic
  G4HadronElasticProcess* theppElasticProcess = new G4HadronElasticProcess();
  G4LElastic* theppElasticModel = new G4LElastic();
  theppElasticProcess-> RegisterMe(theppElasticModel);
  pManager-> AddDiscreteProcess(theppElasticProcess);

  // inelastic
  G4PionPlusInelasticProcess* thePionPlusInelasticProcess
    = new G4PionPlusInelasticProcess(); 

  G4LEPionPlusInelastic* thePionPlusLEPModel = new G4LEPionPlusInelastic();
  G4HEPionPlusInelastic* thePionPlusHEPModel = new G4HEPionPlusInelastic();
  thePionPlusInelasticProcess->RegisterMe(thePionPlusLEPModel);
  thePionPlusInelasticProcess->RegisterMe(thePionPlusHEPModel);
  pManager-> AddDiscreteProcess(thePionPlusInelasticProcess);

  // QED
  G4VProcess* theppMultipleScattering = new G4hMultipleScattering();
  G4VProcess* theppIonisation         = new G4hIonisation();

  pManager-> AddProcess(theppIonisation);
  pManager-> AddProcess(theppMultipleScattering);

  pManager-> SetProcessOrdering(theppMultipleScattering, idxAlongStep, 1);
  pManager-> SetProcessOrdering(theppIonisation,         idxAlongStep, 2);

  pManager-> SetProcessOrdering(theppMultipleScattering, idxPostStep, 1);
  pManager-> SetProcessOrdering(theppIonisation,         idxPostStep, 2);
  
  // ---------------------------------------------------------------
  // pi-
  // ---------------------------------------------------------------
  pManager= G4PionMinus::PionMinus()-> GetProcessManager();

  // elastic
  G4HadronElasticProcess* thepmElasticProcess = new G4HadronElasticProcess();
  G4LElastic* thepmElasticModel = new G4LElastic();
  thepmElasticProcess->RegisterMe(thepmElasticModel);
  pManager-> AddDiscreteProcess(thepmElasticProcess);

  // inelastic
  G4PionMinusInelasticProcess* thePionMinusInelasticProcess 
    = new G4PionMinusInelasticProcess(); 

  G4LEPionMinusInelastic* thePionMinusLEPModel = new G4LEPionMinusInelastic();
  G4HEPionMinusInelastic* thePionMinusHEPModel = new G4HEPionMinusInelastic();
  thePionMinusInelasticProcess-> RegisterMe(thePionMinusLEPModel);
  thePionMinusInelasticProcess-> RegisterMe(thePionMinusHEPModel);
  pManager-> AddDiscreteProcess(thePionMinusInelasticProcess);

  // QED
  G4VProcess* thepmMultipleScattering = new G4hMultipleScattering();
  G4VProcess* thepmIonisation         = new G4hIonisation();

  pManager-> AddProcess(thepmIonisation);
  pManager-> AddProcess(thepmMultipleScattering);

  pManager-> SetProcessOrdering(thepmMultipleScattering, idxAlongStep, 1);
  pManager-> SetProcessOrdering(thepmIonisation,         idxAlongStep, 2);

  pManager-> SetProcessOrdering(thepmMultipleScattering, idxPostStep, 1);
  pManager-> SetProcessOrdering(thepmIonisation,         idxPostStep, 2);

  // ---------------------------------------------------------------
  // K+
  // ---------------------------------------------------------------
  pManager= G4KaonPlus::KaonPlus()-> GetProcessManager();

  // elastic
  G4HadronElasticProcess* thekpElasticProcess = new G4HadronElasticProcess();
  G4LElastic* thekpElasticModel = new G4LElastic();
  thekpElasticProcess->RegisterMe(thekpElasticModel);
  pManager->AddDiscreteProcess(thekpElasticProcess);

  // inelastic
  G4KaonPlusInelasticProcess* theKaonPlusInelasticProcess 
    = new G4KaonPlusInelasticProcess(); 

  G4LEKaonPlusInelastic* theKaonPlusLEPModel = new G4LEKaonPlusInelastic();
  G4HEKaonPlusInelastic* theKaonPlusHEPModel = new G4HEKaonPlusInelastic();
  theKaonPlusInelasticProcess-> RegisterMe(theKaonPlusLEPModel);
  theKaonPlusInelasticProcess-> RegisterMe(theKaonPlusHEPModel);
  pManager-> AddDiscreteProcess(theKaonPlusInelasticProcess);

  // QED
  G4VProcess* thekpMultipleScattering = new G4hMultipleScattering();
  G4VProcess* thekpIonisation         = new G4hIonisation();

  pManager-> AddProcess(thekpIonisation);
  pManager-> AddProcess(thekpMultipleScattering);

  pManager-> SetProcessOrdering(thekpMultipleScattering, idxAlongStep, 1);
  pManager-> SetProcessOrdering(thekpIonisation,         idxAlongStep, 2);

  pManager-> SetProcessOrdering(thekpMultipleScattering, idxPostStep, 1);
  pManager-> SetProcessOrdering(thekpIonisation,         idxPostStep, 2);

  // ---------------------------------------------------------------
  // K-
  // ---------------------------------------------------------------
  pManager= G4KaonMinus::KaonMinus()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* thekmElasticProcess = new G4HadronElasticProcess();
  G4LElastic* thekmElasticModel = new G4LElastic();
  thekmElasticProcess->RegisterMe(thekmElasticModel);
  pManager->AddDiscreteProcess(thekmElasticProcess);

  // inelastic
  G4KaonMinusInelasticProcess* theKaonMinusInelasticProcess 
    = new G4KaonMinusInelasticProcess(); 

  G4LEKaonMinusInelastic* theKaonMinusLEPModel = new G4LEKaonMinusInelastic();
  G4HEKaonMinusInelastic* theKaonMinusHEPModel = new G4HEKaonMinusInelastic();
  theKaonMinusInelasticProcess->RegisterMe(theKaonMinusLEPModel);
  theKaonMinusInelasticProcess->RegisterMe(theKaonMinusHEPModel);
  pManager->AddDiscreteProcess(theKaonMinusInelasticProcess);

  // QED
  G4VProcess* thekmMultipleScattering = new G4hMultipleScattering();
  G4VProcess* thekmIonisation         = new G4hIonisation();

  pManager-> AddProcess(thekmIonisation);
  pManager-> AddProcess(thekmMultipleScattering);

  pManager->SetProcessOrdering(thekmMultipleScattering, idxAlongStep, 1);
  pManager->SetProcessOrdering(thekmIonisation,         idxAlongStep, 2);

  pManager->SetProcessOrdering(thekmMultipleScattering, idxPostStep, 1);
  pManager->SetProcessOrdering(thekmIonisation,         idxPostStep, 2);

  // ---------------------------------------------------------------
  // K0L
  // ---------------------------------------------------------------
  pManager= G4KaonZeroLong::KaonZeroLong()-> GetProcessManager();

  // elastic
  G4HadronElasticProcess* thek0lElasticProcess = new G4HadronElasticProcess();
  G4LElastic* thek0lElasticModel = new G4LElastic();
  thek0lElasticProcess-> RegisterMe(thek0lElasticModel);
  pManager-> AddDiscreteProcess(thek0lElasticProcess);

  // inelastic
  G4KaonZeroLInelasticProcess* theKaonZeroLInelasticProcess 
    = new G4KaonZeroLInelasticProcess(); 

  G4LEKaonZeroLInelastic* theKaonZeroLLEPModel = new G4LEKaonZeroLInelastic();
  G4HEKaonZeroInelastic* theKaonZerolHEPModel  = new G4HEKaonZeroInelastic();

  theKaonZeroLInelasticProcess-> RegisterMe(theKaonZeroLLEPModel);
  theKaonZeroLInelasticProcess-> RegisterMe(theKaonZerolHEPModel);

  pManager-> AddDiscreteProcess(theKaonZeroLInelasticProcess);

  // ---------------------------------------------------------------
  // K0S
  // ---------------------------------------------------------------
  pManager= G4KaonZeroShort::KaonZeroShort()-> GetProcessManager();

  // elastic
  G4HadronElasticProcess* thek0sElasticProcess = new G4HadronElasticProcess();
  G4LElastic* thek0sElasticModel = new G4LElastic();
  thek0sElasticProcess-> RegisterMe(thek0sElasticModel);
  pManager-> AddDiscreteProcess(thek0sElasticProcess);

  // inelastic
  G4KaonZeroSInelasticProcess* theKaonZeroSInelasticProcess 
    = new G4KaonZeroSInelasticProcess(); 

  G4LEKaonZeroSInelastic* theKaonZeroSLEPModel = new G4LEKaonZeroSInelastic();
  G4HEKaonZeroInelastic* theKaonZerosHEPModel  = new G4HEKaonZeroInelastic();

  theKaonZeroSInelasticProcess-> RegisterMe(theKaonZeroSLEPModel);
  theKaonZeroSInelasticProcess-> RegisterMe(theKaonZerosHEPModel);

  pManager-> AddDiscreteProcess(theKaonZeroSInelasticProcess);

}

