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

#include <iomanip>   

#include "GammaRayTelHadronPhysics.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4HadronicParameters.hh"
#include "G4ProcessManager.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"

#include "G4HadronElasticProcess.hh"
#include "G4NeutronFissionProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4HadronElastic.hh"

#include "G4HadronInelasticProcess.hh"

// Low-energy Models
#include "G4LFission.hh"

// Cross section handlers and high energy models
#include "G4VCrossSectionDataSet.hh"
#include "G4CascadeInterface.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4BGGNucleonInelasticXS.hh"

// Stopping processes
#include "G4HadronStoppingProcess.hh"
#include "G4HadronicAbsorptionBertini.hh"
#include "G4HadronicAbsorptionFritiof.hh"

// quark gluon string model with chips afterburner.
#include "G4FTFModel.hh"
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"


GammaRayTelHadronPhysics::GammaRayTelHadronPhysics(const G4String& name)
                    :  G4VPhysicsConstructor(name)
{;}

GammaRayTelHadronPhysics::~GammaRayTelHadronPhysics()
{}


void GammaRayTelHadronPhysics::ConstructProcess()
{  
  G4ProcessManager * pManager = 0;
  /*
  G4cout << "" << G4endl;
  G4cout << "You are using the GammaRayTelHadronPhysics" << G4endl;
  G4cout << " - Note that this hadronic physics list is not optimized for any particular usage" << G4endl;
  G4cout << " - If you wish to have a starting point tailored for a particular area of work," << G4endl;
  G4cout << "   please use one of the available physics lists by use-case." << G4endl;
  G4cout << "" << G4endl;
  */

  // Elastic Process
  G4HadronElastic* theElasticModel = new G4HadronElastic;
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  theElasticProcess->RegisterMe(theElasticModel);

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

  G4TheoFSGenerator* theModelDownToZero =  new G4TheoFSGenerator( "FTFP" );
  theModelDownToZero->SetHighEnergyGenerator( theStringModel );
  theModelDownToZero->SetTransport( theCascade );
  theModelDownToZero->SetMinEnergy(0*eV);
  theModelDownToZero->SetMaxEnergy(theFTFMax );

  G4CascadeInterface * theBERTModel = new G4CascadeInterface;
  theBERTModel->SetMinEnergy( theBERTMin );
  theBERTModel->SetMaxEnergy( theBERTMax );

  // pi+ and pi-      

  G4VCrossSectionDataSet * thePiPlusData = new G4BGGPionInelasticXS( G4PionPlus::Definition() ); 
  G4VCrossSectionDataSet * thePiMinusData = new G4BGGPionInelasticXS( G4PionMinus::Definition() ); 
  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );

  // PionPlus
  G4ParticleDefinition* pion = G4PionPlus::PionPlusDefinition();
  pManager = pion ->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* thePionPlusInelastic =
    new G4HadronInelasticProcess( "piInelastic", G4PionPlus::Definition() );
  thePionPlusInelastic->AddDataSet(thePiPlusData);
  thePionPlusInelastic->RegisterMe(theBERTModel);
  thePionPlusInelastic->RegisterMe(theModel);
  pManager->AddDiscreteProcess(thePionPlusInelastic);

  pManager->AddProcess(new G4hIonisation, ordInActive,2, 2);

  G4hMultipleScattering* thePionPlusMult = new G4hMultipleScattering;
  pManager->AddProcess(thePionPlusMult);
  pManager->SetProcessOrdering(thePionPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(thePionPlusMult, idxPostStep, 1);

  // PionMinus
 G4ParticleDefinition* pionMinus = G4PionMinus::PionMinusDefinition();
  pManager = pionMinus -> GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* thePionMinusInelastic =
    new G4HadronInelasticProcess( "piInelastic", G4PionMinus::Definition() );  
  thePionMinusInelastic->AddDataSet(thePiMinusData);
  thePionMinusInelastic->RegisterMe(theBERTModel);
  thePionMinusInelastic->RegisterMe(theModel);
  pManager->AddDiscreteProcess(thePionMinusInelastic);

  pManager->AddProcess(new G4hIonisation, ordInActive,2, 2);

  G4hMultipleScattering* thePionMinusMult = new G4hMultipleScattering;
  pManager->AddProcess(thePionMinusMult);
  pManager->SetProcessOrdering(thePionMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(thePionMinusMult, idxPostStep, 1);

  pManager->AddRestProcess(new G4HadronicAbsorptionBertini(G4PionMinus::Definition()), ordDefault);

  // KaonPlus
  G4ParticleDefinition* kaonPlus = G4KaonPlus::KaonPlusDefinition();
  pManager = kaonPlus->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* theKaonPlusInelastic =
    new G4HadronInelasticProcess( "kaon+Inelastic", G4KaonPlus::Definition() );
  theKaonPlusInelastic->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
			 	    GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
  theKaonPlusInelastic->RegisterMe(theBERTModel);
  theKaonPlusInelastic->RegisterMe(theModel);
  pManager->AddDiscreteProcess(theKaonPlusInelastic);

  pManager->AddProcess(new G4hIonisation, ordInActive,2, 2);

  G4hMultipleScattering* theKaonPlusMult = new G4hMultipleScattering;
  pManager->AddProcess(theKaonPlusMult);
  pManager->SetProcessOrdering(theKaonPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theKaonPlusMult, idxPostStep, 1);

  // KaonMinus
  G4ParticleDefinition* kaonMinus = G4KaonMinus::KaonMinusDefinition();
  pManager = kaonMinus->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* theKaonMinusInelastic =
    new G4HadronInelasticProcess( "kaon-Inelastic", G4KaonMinus::Definition() );
  theKaonMinusInelastic->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
			  	     GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
  theKaonMinusInelastic->RegisterMe(theBERTModel);
  theKaonMinusInelastic->RegisterMe(theModel);
  pManager->AddDiscreteProcess(theKaonMinusInelastic);

  pManager->AddProcess(new G4hIonisation, ordInActive,2, 2);

  G4hMultipleScattering* theKaonMinusMult = new G4hMultipleScattering;
  pManager->AddProcess(theKaonMinusMult);
  pManager->SetProcessOrdering(theKaonMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theKaonMinusMult, idxPostStep, 1);

  pManager->AddRestProcess(new G4HadronicAbsorptionBertini(G4KaonMinus::Definition()), ordDefault);

  // KaonZeroL
  pManager = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* theKaonZeroLInelastic =
    new G4HadronInelasticProcess( "kaon0LInelastic", G4KaonZeroLong::Definition() );  
  theKaonZeroLInelastic->AddDataSet( G4CrossSectionDataSetRegistry::Instance()
				     ->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
  theKaonZeroLInelastic->RegisterMe(theBERTModel);
  theKaonZeroLInelastic->RegisterMe(theModel);
  pManager->AddDiscreteProcess(theKaonZeroLInelastic);
 
  // KaonZeroS
  pManager = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* theKaonZeroSInelastic =
    new G4HadronInelasticProcess( "kaon0SInelastic", G4KaonZeroShort::Definition() );  
  theKaonZeroSInelastic->AddDataSet( G4CrossSectionDataSetRegistry::Instance()
				     ->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
  theKaonZeroSInelastic->RegisterMe(theBERTModel);
  theKaonZeroSInelastic->RegisterMe(theModel); 
  pManager->AddDiscreteProcess(theKaonZeroSInelastic);

  // Proton
  pManager = G4Proton::Proton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* theProtonInelastic =
    new G4HadronInelasticProcess( "protonInelastic", G4Proton::Definition() );  
  theProtonInelastic->AddDataSet(new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
  theProtonInelastic->RegisterMe(theBERTModel);
  theProtonInelastic->RegisterMe(theModel);
  pManager->AddDiscreteProcess(theProtonInelastic);

  pManager->AddProcess(new G4hIonisation, ordInActive,2, 2);

  G4hMultipleScattering* theProtonMult = new G4hMultipleScattering;
  pManager->AddProcess(theProtonMult);
  pManager->SetProcessOrdering(theProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theProtonMult, idxPostStep, 1);

  // anti-Proton
  pManager = G4AntiProton::AntiProton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* theAntiProtonInelastic =
    new G4HadronInelasticProcess( "anti_protonInelastic", G4AntiProton::Definition() );  
  theAntiProtonInelastic->AddDataSet( theAntiNucleonData );
  theAntiProtonInelastic->RegisterMe(theModelDownToZero);
  pManager->AddDiscreteProcess(theAntiProtonInelastic);

  pManager->AddProcess(new G4hIonisation, ordInActive,2, 2);

  G4hMultipleScattering* theAntiProtonMult = new G4hMultipleScattering;
  pManager->AddProcess(theAntiProtonMult);
  pManager->SetProcessOrdering(theAntiProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(theAntiProtonMult, idxPostStep, 1);

  pManager->AddRestProcess(new G4HadronicAbsorptionFritiof(G4AntiProton::Definition()));

  // Neutron
  pManager = G4Neutron::Neutron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* theNeutronInelastic =
    new G4HadronInelasticProcess( "neutronInelastic", G4Neutron::Definition() );  
  theNeutronInelastic->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ));
  theNeutronInelastic->RegisterMe(theBERTModel);
  theNeutronInelastic->RegisterMe(theModel);
  pManager->AddDiscreteProcess(theNeutronInelastic);

  G4NeutronFissionProcess* theNeutronFission = new G4NeutronFissionProcess;
  G4LFission* theNeutronFissionModel = new G4LFission;
  theNeutronFission->RegisterMe(theNeutronFissionModel);
  pManager->AddDiscreteProcess(theNeutronFission);

  G4NeutronCaptureProcess* theNeutronCapture = new G4NeutronCaptureProcess;
  theNeutronCapture->AddDataSet(new G4NeutronCaptureXS()); 
  pManager->AddDiscreteProcess(theNeutronCapture);

  // AntiNeutron
  pManager = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theElasticProcess);
  G4HadronInelasticProcess* theAntiNeutronInelastic =
    new G4HadronInelasticProcess( "anti_neutronInelastic", G4AntiNeutron::Definition() );    
  theAntiNeutronInelastic->AddDataSet( theAntiNucleonData );
  theAntiNeutronInelastic->RegisterMe(theModelDownToZero);
  pManager->AddDiscreteProcess(theAntiNeutronInelastic);

}
