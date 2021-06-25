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
//
// 

#include <iomanip>   

#include "GammaRayTelIonPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicParameters.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronElastic.hh"

#include "G4HadronInelasticProcess.hh"

#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"

#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4FTFModel.hh"
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4CascadeInterface.hh"


GammaRayTelIonPhysics::GammaRayTelIonPhysics(const G4String& name)
                 :  G4VPhysicsConstructor(name)
{;}

GammaRayTelIonPhysics::~GammaRayTelIonPhysics()
{;}

void GammaRayTelIonPhysics::ConstructParticle()
{;}


#include "G4ProcessManager.hh"


void GammaRayTelIonPhysics::ConstructProcess()
{

  
  /*

   // Deuteron physics
   G4HadronInelasticProcess  fDeuteronProcess;

   // Triton physics
   G4HadronInelasticProcess    fTritonProcess;
  
   // Alpha physics
   G4HadronInelasticProcess     fAlphaProcess;


   */


  
  G4ProcessManager * pManager = 0;
  
  const G4double theBERTMin =   0.0*GeV;
  const G4double theBERTMax =   5.0*GeV;
  const G4double theFTFMin =    4.0*GeV;
  const G4double theFTFMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  
  G4FTFModel* theStringModel = new G4FTFModel;
  G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay();
  theStringModel->SetFragmentationModel( theStringDecay );
  G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface();

  G4TheoFSGenerator* theModel = new G4TheoFSGenerator( "FTFP" );
  theModel->SetHighEnergyGenerator( theStringModel );
  theModel->SetTransport( theCascade );
  theModel->SetMinEnergy( theFTFMin );
  theModel->SetMaxEnergy( theFTFMax ); 

  G4CascadeInterface * theBERTModel = new G4CascadeInterface;
  theBERTModel->SetMinEnergy( theBERTMin );
  theBERTModel->SetMaxEnergy( theBERTMax );

  // Elastic Process
  G4HadronElastic* theElasticModel = new G4HadronElastic;
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  theElasticProcess->RegisterMe(theElasticModel);

  // Generic Ion
  pManager = G4GenericIon::GenericIon()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);

  pManager->AddProcess(new G4hIonisation, ordInActive, 2, 2);

  G4hMultipleScattering* fIonMultipleScattering = new G4hMultipleScattering;
  pManager->AddProcess(fIonMultipleScattering);
  pManager->SetProcessOrdering(fIonMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fIonMultipleScattering, idxPostStep,  1);

  G4ComponentGGNuclNuclXsc * ggNuclNuclXsec = new G4ComponentGGNuclNuclXsc();
  G4VCrossSectionDataSet * theGGNuclNuclData = new G4CrossSectionInelastic(ggNuclNuclXsec);
  
  // Deuteron 
  pManager = G4Deuteron::Deuteron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* fDeuteronProcess =
    new G4HadronInelasticProcess( "dInelastic", G4Deuteron::Definition() );
  fDeuteronProcess->AddDataSet(theGGNuclNuclData);
  fDeuteronProcess->RegisterMe(theBERTModel);
  fDeuteronProcess->RegisterMe(theModel);
  pManager->AddDiscreteProcess(fDeuteronProcess);

  pManager->AddProcess(new G4hIonisation, ordInActive, 2, 2);

  G4hMultipleScattering* fDeuteronMultipleScattering = new G4hMultipleScattering;
  pManager->AddProcess(fDeuteronMultipleScattering);
  pManager->SetProcessOrdering(fDeuteronMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fDeuteronMultipleScattering, idxPostStep,  1);
 
  // Triton 
  pManager = G4Triton::Triton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* fTritonProcess =
    new G4HadronInelasticProcess( "tInelastic", G4Triton::Definition() );  
  fTritonProcess->AddDataSet(theGGNuclNuclData);
  fTritonProcess->RegisterMe(theBERTModel);
  fTritonProcess->RegisterMe(theModel);
  pManager->AddDiscreteProcess(fTritonProcess);

  pManager->AddProcess(new G4hIonisation, ordInActive, 2, 2);

  G4hMultipleScattering* fTritonMultipleScattering = new G4hMultipleScattering;
  pManager->AddProcess(fTritonMultipleScattering);
  pManager->SetProcessOrdering(fTritonMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fTritonMultipleScattering, idxPostStep,  1);
 
  // Alpha 
  pManager = G4Alpha::Alpha()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* fAlphaProcess =
    new G4HadronInelasticProcess( "alphaInelastic", G4Alpha::Definition() );  
  fAlphaProcess->AddDataSet(theGGNuclNuclData);
  fAlphaProcess->RegisterMe(theBERTModel);
  fAlphaProcess->RegisterMe(theModel);
  pManager->AddDiscreteProcess(fAlphaProcess);

  pManager->AddProcess(new G4hIonisation, ordInActive, 2, 2);

  G4hMultipleScattering* fAlphaMultipleScattering = new G4hMultipleScattering;
  pManager->AddProcess(fAlphaMultipleScattering);
  pManager->SetProcessOrdering(fAlphaMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fAlphaMultipleScattering, idxPostStep,  1);
 
  // He3
  pManager = G4He3::He3()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);

  pManager->AddProcess(new G4hIonisation, ordInActive, 2, 2);

  G4hMultipleScattering* fHe3MultipleScattering = new G4hMultipleScattering;
  pManager->AddProcess(fHe3MultipleScattering);
  pManager->SetProcessOrdering(fHe3MultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(fHe3MultipleScattering, idxPostStep,  1);
   
}



