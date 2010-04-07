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
//
// $Id: RE01IonPhysics.cc,v 1.3 2010-04-07 01:27:53 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RE01IonPhysics.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4IonConstructor.hh"

// processes
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4HadronElasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// models
#include "G4LElastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"

RE01IonPhysics::RE01IonPhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{;}


RE01IonPhysics::~RE01IonPhysics()
{;}


void RE01IonPhysics::ConstructParticle()
{ 
  // Construct light ions (d, t, 3He, alpha, and generic ion)
  G4IonConstructor ionConstruct;
  ionConstruct.ConstructParticle();
}


void RE01IonPhysics::ConstructProcess()
{
  // Hadronic Elastic Process and Model (for all ions except generic ion)

  G4HadronElasticProcess* elasticProcess = new G4HadronElasticProcess();
  G4LElastic* elasticModel = new G4LElastic();
  elasticProcess->RegisterMe(elasticModel);

  // Hadronic inelastic models

  G4ProcessManager * pManager = 0;
  
  ///////////////////
  //               //
  //   Deuteron    //
  //               //
  ///////////////////

  pManager = G4Deuteron::Deuteron()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4DeuteronInelasticProcess* dinelProc = new G4DeuteronInelasticProcess();
  G4LEDeuteronInelastic* LEPdModel = new G4LEDeuteronInelastic();
  dinelProc->RegisterMe(LEPdModel);
  pManager->AddDiscreteProcess(dinelProc);

  ///////////////////
  //               //
  //    Triton     //
  //               //
  ///////////////////

  pManager = G4Triton::Triton()->GetProcessManager(); 

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4TritonInelasticProcess* tinelProc = new G4TritonInelasticProcess();
  G4LETritonInelastic* LEPtModel = new G4LETritonInelastic();
  tinelProc->RegisterMe(LEPtModel);
  pManager->AddDiscreteProcess(tinelProc);

  ///////////////////
  //               //
  //      3He      //
  //               //
  ///////////////////

  pManager = G4He3::He3()->GetProcessManager(); 

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic 
  pManager->AddDiscreteProcess(elasticProcess);

  // NO INELASTIC PROCESS AVAILABLE FOR 3HE

  ///////////////////
  //               //
  //     Alpha     //
  //               //
  ///////////////////

  pManager = G4Alpha::Alpha()->GetProcessManager(); 

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);
  
  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AlphaInelasticProcess* ainelProc = new G4AlphaInelasticProcess();
  G4LEAlphaInelastic* LEPaModel = new G4LEAlphaInelastic();
  ainelProc->RegisterMe(LEPaModel);
  pManager->AddDiscreteProcess(ainelProc);

  ///////////////////
  //               //
  //  generic ion  //
  //               //
  ///////////////////

  pManager = G4GenericIon::GenericIon()->GetProcessManager();

  // Only EM processes for generic ion
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4ionIonisation(),      -1, 2, 2);
 
}
