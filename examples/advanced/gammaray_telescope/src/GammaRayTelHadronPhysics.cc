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


GammaRayTelHadronPhysics::GammaRayTelHadronPhysics(const G4String& name)
                    :  G4VPhysicsConstructor(name)
{;}

GammaRayTelHadronPhysics::~GammaRayTelHadronPhysics()
{
  delete theStringDecay;
}


#include "G4ProcessManager.hh"


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
  theElasticModel = new G4HadronElastic();
  theElasticProcess.RegisterMe(theElasticModel);

  const G4double theBERTMin =   0.0*GeV;
  const G4double theBERTMax =   5.0*GeV;
  const G4double theFTFMin =    4.0*GeV;
  const G4double theFTFMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  
  theStringModel = new G4FTFModel;
  theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  theModel = new G4TheoFSGenerator( "FTFP" );
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

  G4VCrossSectionDataSet * thePiData = new G4CrossSectionPairGG( new G4PiNuclearCrossSection, 91*GeV );
  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );


  // PionPlus
  G4ParticleDefinition* pion = G4PionPlus::PionPlusDefinition();
  pManager = pion ->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  thePionPlusInelastic.AddDataSet(thePiData);
  thePionPlusInelastic.RegisterMe(theBERTModel);
  thePionPlusInelastic.RegisterMe(theModel);
  pManager->AddDiscreteProcess(&thePionPlusInelastic);

  pManager->AddProcess(&thePionPlusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&thePionPlusMult);
  pManager->SetProcessOrdering(&thePionPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&thePionPlusMult, idxPostStep, 1);

  // PionMinus
 G4ParticleDefinition* pionMinus = G4PionMinus::PionMinusDefinition();
  pManager = pionMinus -> GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  thePionMinusInelastic.AddDataSet(thePiData);
  thePionMinusInelastic.RegisterMe(theBERTModel);
  thePionMinusInelastic.RegisterMe(theModel);
  pManager->AddDiscreteProcess(&thePionMinusInelastic);

  pManager->AddProcess(&thePionMinusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&thePionMinusMult);
  pManager->SetProcessOrdering(&thePionMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&thePionMinusMult, idxPostStep, 1);

  pManager->AddRestProcess(&thePionMinusAbsorption, ordDefault);

  // KaonPlus
  G4ParticleDefinition* kaonPlus = G4KaonPlus::KaonPlusDefinition();
  pManager = kaonPlus->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  theKaonPlusInelastic.AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
				   GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
  theKaonPlusInelastic.RegisterMe(theBERTModel);
  theKaonPlusInelastic.RegisterMe(theModel);
  pManager->AddDiscreteProcess(&theKaonPlusInelastic);

  pManager->AddProcess(&theKaonPlusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theKaonPlusMult);
  pManager->SetProcessOrdering(&theKaonPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theKaonPlusMult, idxPostStep, 1);

  // KaonMinus
  G4ParticleDefinition* kaonMinus = G4KaonMinus::KaonMinusDefinition();
  pManager = kaonMinus->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  theKaonMinusInelastic.AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
				   GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
  theKaonMinusInelastic.RegisterMe(theBERTModel);
  theKaonMinusInelastic.RegisterMe(theModel);
  pManager->AddDiscreteProcess(&theKaonMinusInelastic);

  pManager->AddProcess(&theKaonMinusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theKaonMinusMult);
  pManager->SetProcessOrdering(&theKaonMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theKaonMinusMult, idxPostStep, 1);

  pManager->AddRestProcess(&theKaonMinusAbsorption, ordDefault);

  // KaonZeroL
  pManager = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  theKaonZeroLInelastic.AddDataSet( G4CrossSectionDataSetRegistry::Instance()
				    ->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
  theKaonZeroLInelastic.RegisterMe(theBERTModel);
  theKaonZeroLInelastic.RegisterMe(theModel);
  pManager->AddDiscreteProcess(&theKaonZeroLInelastic);
 
  // KaonZeroS
  pManager = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  theKaonZeroSInelastic.AddDataSet( G4CrossSectionDataSetRegistry::Instance()
				    ->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
  theKaonZeroSInelastic.RegisterMe(theBERTModel);
  theKaonZeroSInelastic.RegisterMe(theModel); 
  pManager->AddDiscreteProcess(&theKaonZeroSInelastic);

  // Proton
  pManager = G4Proton::Proton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  theProtonInelastic.AddDataSet(new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
  theProtonInelastic.RegisterMe(theBERTModel);
  theProtonInelastic.RegisterMe(theModel);
  pManager->AddDiscreteProcess(&theProtonInelastic);

  pManager->AddProcess(&theProtonIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theProtonMult);
  pManager->SetProcessOrdering(&theProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theProtonMult, idxPostStep, 1);

  // anti-Proton
  pManager = G4AntiProton::AntiProton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  theAntiProtonInelastic.AddDataSet( theAntiNucleonData );
  theAntiProtonInelastic.RegisterMe(theModelDownToZero);
  pManager->AddDiscreteProcess(&theAntiProtonInelastic);

  pManager->AddProcess(&theAntiProtonIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theAntiProtonMult);
  pManager->SetProcessOrdering(&theAntiProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theAntiProtonMult, idxPostStep, 1);

  pManager->AddRestProcess(&theAntiProtonAnnihilation);

  // Neutron
  pManager = G4Neutron::Neutron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  theNeutronInelastic.AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ));
  theNeutronInelastic.RegisterMe(theBERTModel);
  theNeutronInelastic.RegisterMe(theModel);
  pManager->AddDiscreteProcess(&theNeutronInelastic);
  
  theNeutronFissionModel = new G4LFission();
  theNeutronFission.RegisterMe(theNeutronFissionModel);
  pManager->AddDiscreteProcess(&theNeutronFission);

  theNeutronCapture = new G4HadronCaptureProcess();
  theNeutronCapture->AddDataSet(new G4NeutronCaptureXS()); 
  pManager->AddDiscreteProcess(theNeutronCapture);

  // AntiNeutron
  pManager = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);  
  theAntiNeutronInelastic.AddDataSet( theAntiNucleonData );
  theAntiNeutronInelastic.RegisterMe(theModelDownToZero);
  pManager->AddDiscreteProcess(&theAntiNeutronInelastic);

}
