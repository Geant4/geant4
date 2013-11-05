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
/// \file runAndEvent/RE02/src/RE02IonPhysics.cc
/// \brief Implementation of the RE02IonPhysics class
//
// $Id$
// --------------------------------------------------------------
//
// 13-Oct-2003 Add Comment for Ionisation of Generic Ion by T. Koi 
// 05-Jan-2004 Change G. Ion Ionisation from G4hIonisation
//                                      to G4ionIonisation T. Koi
// 18-Nov-2005 Add Inelastic process with G4BinaryLightIonReaction.
//                                            T. Aso

#include "RE02IonPhysics.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE02IonPhysics::RE02IonPhysics(const G4String& name)
                 :  G4VPhysicsConstructor(name)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE02IonPhysics::~RE02IonPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4ProcessManager.hh"
#include "G4HadronElastic.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE02IonPhysics::ConstructProcess()
{
   G4ProcessManager * pManager = 0;

   // Generic Ion
   pManager = G4GenericIon::GenericIon()->GetProcessManager();

   // add process
   G4VProcess* thegionMultipleScattering = new G4hMultipleScattering();
   //
   // G4hIonization may be not able to use for Geanric Ion in future 
   // Please take care using this physics list after v5.2.p02 
   // G4VProcess* thegionIonisation        = new G4hIonisation();
   //
   // From V6.0 hIonisation does not work for GenericIon 
   G4VProcess* thegionIonisation        = new G4ionIonisation();

   // Inelastic process
   G4HadronInelasticProcess* thegionInelastic =
       new G4HadronInelasticProcess("IonInelastic",G4GenericIon::GenericIon());
   thegionInelastic->AddDataSet(new G4TripathiCrossSection);
   thegionInelastic->AddDataSet(new G4IonsShenCrossSection);
   G4BinaryLightIonReaction* thegionBCModel = new G4BinaryLightIonReaction;
   thegionInelastic->RegisterMe(thegionBCModel);

   // Inelastic for light ions
   G4BinaryLightIonReaction* binaryCascade = new G4BinaryLightIonReaction;
   binaryCascade->SetMinEnergy(0.0);
   binaryCascade->SetMaxEnergy(110*MeV);

   G4QMDReaction* qmd = new G4QMDReaction;
   qmd->SetMinEnergy(100*MeV);
   qmd->SetMaxEnergy(10*GeV);

   G4IonsShenCrossSection* shenXS = new G4IonsShenCrossSection;
   G4TripathiCrossSection* tripXS = new G4TripathiCrossSection;
   G4TripathiLightCrossSection* tripLightXS = new G4TripathiLightCrossSection;

   //
   pManager->AddProcess(thegionIonisation);
   pManager->AddProcess(thegionMultipleScattering);
   pManager->AddDiscreteProcess(thegionInelastic);
   //
   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thegionMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thegionIonisation,         idxAlongStep,2);
   //
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thegionMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thegionIonisation,         idxPostStep,2);

   // Deuteron
   pManager = G4Deuteron::Deuteron()->GetProcessManager();

   // add process
   G4HadronElasticProcess* thedueElasticProcess = new G4HadronElasticProcess();
   thedueElasticProcess->RegisterMe(new G4HadronElastic());
   pManager->AddDiscreteProcess(thedueElasticProcess);

   G4DeuteronInelasticProcess* theDeuteronInelasticProcess
                                            = new G4DeuteronInelasticProcess();
   theDeuteronInelasticProcess->RegisterMe(binaryCascade);
   theDeuteronInelasticProcess->RegisterMe(qmd);
   theDeuteronInelasticProcess->AddDataSet(shenXS);
   theDeuteronInelasticProcess->AddDataSet(tripXS);
   theDeuteronInelasticProcess->AddDataSet(tripLightXS);
   pManager->AddDiscreteProcess(theDeuteronInelasticProcess);

   G4VProcess* thedueMultipleScattering = new G4hMultipleScattering();
   G4VProcess* thedueIonisation         = new G4hIonisation();
   //
   pManager->AddProcess(thedueIonisation);
   pManager->AddProcess(thedueMultipleScattering);
   //
   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thedueMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thedueIonisation,         idxAlongStep,2);
   //
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thedueMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thedueIonisation,         idxPostStep,2);

   // Triton
   pManager = G4Triton::Triton()->GetProcessManager();

   // add process
   G4HadronElasticProcess* thetriElasticProcess = new G4HadronElasticProcess();
   thetriElasticProcess->RegisterMe(new G4HadronElastic());
   pManager->AddDiscreteProcess(thetriElasticProcess);

   G4TritonInelasticProcess* theTritonInelasticProcess =
                                                new G4TritonInelasticProcess();
   theTritonInelasticProcess->RegisterMe(binaryCascade);
   theTritonInelasticProcess->RegisterMe(qmd);
   theTritonInelasticProcess->AddDataSet(shenXS);
   theTritonInelasticProcess->AddDataSet(tripXS);
   theTritonInelasticProcess->AddDataSet(tripLightXS);
   pManager->AddDiscreteProcess(theTritonInelasticProcess);

   G4VProcess* thetriMultipleScattering = new G4hMultipleScattering();
   G4VProcess* thetriIonisation         = new G4hIonisation();
   //
   pManager->AddProcess(thetriIonisation);
   pManager->AddProcess(thetriMultipleScattering);
   //
   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thetriMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thetriIonisation,         idxAlongStep,2);
   //
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thetriMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thetriIonisation,         idxPostStep,2);

   // Alpha
   pManager = G4Alpha::Alpha()->GetProcessManager();

   // add processes
   G4HadronElasticProcess* thealElasticProcess = new G4HadronElasticProcess();
   thealElasticProcess->RegisterMe(new G4HadronElastic());
   pManager->AddDiscreteProcess(thealElasticProcess);

   G4AlphaInelasticProcess* theAlphaInelasticProcess =
                                                new G4AlphaInelasticProcess();
   theAlphaInelasticProcess->RegisterMe(binaryCascade);
   theAlphaInelasticProcess->RegisterMe(qmd);
   theAlphaInelasticProcess->AddDataSet(shenXS);
   theAlphaInelasticProcess->AddDataSet(tripXS);
   theAlphaInelasticProcess->AddDataSet(tripLightXS);
   pManager->AddDiscreteProcess(theAlphaInelasticProcess);

   G4VProcess* thealpMultipleScattering = new G4hMultipleScattering();
   G4VProcess* thealpIonisation         = new G4hIonisation();
   //
   pManager->AddProcess(thealpIonisation);
   pManager->AddProcess(thealpMultipleScattering);
   //
   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thealpMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thealpIonisation,         idxAlongStep,2);
   //
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thealpMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thealpIonisation,         idxPostStep,2);

   // He3
   pManager = G4He3::He3()->GetProcessManager();

   // add processes
   G4VProcess* thehe3MultipleScattering = new G4hMultipleScattering();
   G4VProcess* thehe3Ionisation         = new G4hIonisation();
   //
   pManager->AddProcess(thehe3Ionisation);
   pManager->AddProcess(thehe3MultipleScattering);
   //
   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thehe3MultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thehe3Ionisation,         idxAlongStep,2);
   //
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thehe3MultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thehe3Ionisation,         idxPostStep,2);

}



