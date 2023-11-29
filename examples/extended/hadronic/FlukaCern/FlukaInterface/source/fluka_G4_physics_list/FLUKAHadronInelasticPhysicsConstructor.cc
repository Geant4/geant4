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
#ifdef G4_USE_FLUKA


#include "FLUKAHadronInelasticPhysicsConstructor.hh"

// G4
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
// G4
#include "G4ParticleDefinition.hh"
#include "G4PhysicsListHelper.hh"
// G4
#include "G4HadParticles.hh"
#include "G4HadronicParameters.hh"
// G4
#include "G4NeutronCaptureProcess.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPCapture.hh"
#include "G4NeutronRadCapture.hh"
// G4
#include "G4HadronInelasticProcess.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPInelastic.hh"
// G4
#include "G4NeutronFissionProcess.hh"
#include "G4ParticleHPFissionData.hh"
#include "G4ParticleHPFission.hh"
#include "G4LFission.hh"

#include "build_G4_process_helpers.hh"

#include "FLUKAInelasticScatteringXS.hh"
#include "FLUKANuclearInelasticModel.hh"


// DEBUG: CHERRY PICK FTFP_BERT XS / MODELS INSTEAD
//#include "G4BGGNucleonInelasticXS.hh"
//#include "G4TheoFSGenerator.hh"
//#include "G4GeneratorPrecompoundInterface.hh"
//#include "G4FTFModel.hh"
//#include "G4ExcitedStringDecay.hh"
//#include "G4QuasiElasticChannel.hh"


// ***************************************************************************
// Construct hadron inelastic physics processes with FLUKA.CERN XS and models.
// ***************************************************************************
FLUKAHadronInelasticPhysicsConstructor::FLUKAHadronInelasticPhysicsConstructor(G4int verbose)
:  G4VPhysicsConstructor("hInelastic FLUKA"),
// HP only:
  fNeutronHPMaxE(20*CLHEP::MeV),
//fNeutronFLUKAMinE(19.9*CLHEP::MeV),
  fNeutronFLUKAMinE(20*CLHEP::MeV)
{
  if (verbose > 1) { 
    G4cout << "### FLUKA Hadron Inelastic Physics" << G4endl;
  }

  const auto param = G4HadronicParameters::Instance();
  param->SetEnableBCParticles(true);
}


// ***************************************************************************
// Construct particles. 
// ***************************************************************************
void FLUKAHadronInelasticPhysicsConstructor::ConstructParticle() {
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}


// ***************************************************************************
// For each particle of interest, 
// processes are created, assigned XS and models, and registered to the process manager.
//
// IMPORTANT NB: The XS (G4VCrossSectionDataSet), models (G4HadronicInteraction), 
// and even processes (G4HadronInelasticProcess) 
// are constructed in a similar way as for any other G4 physics list in G4 source code.
// They do not seem to be OWNED by the G4CrossSectionDataStore, G4EnergyRangeManager 
// and G4ProcessManager respectively, HENCE THEY ARE NEVER DELETED
// (true for ANY physics list).
// Should not matter too much though, because the destructions 
// should have happened at the very end of the run anyway.
// ***************************************************************************
void FLUKAHadronInelasticPhysicsConstructor::ConstructProcess() {

  const auto helper = G4PhysicsListHelper::GetPhysicsListHelper();

  // FLUKA hadron - nucleus inelastic XS
  const auto flukaInelasticScatteringXS = new FLUKAInelasticScatteringXS();

  // FLUKA hadron - nucleus model
  const auto flukaModel = new FLUKANuclearInelasticModel();


  // PROTON
  build_G4_process_helpers::buildInelasticProcess(G4Proton::Proton(),
                                                  helper,
                                                  flukaInelasticScatteringXS,
                                                  flukaModel);
	
  // DEBUG: CHERRY PICK G4 XS / MODELS INSTEAD
  /*auto protonInelasticProcess = new G4HadronInelasticProcess("protonInelastic", G4Proton::Proton());
    helper->RegisterProcess(protonInelasticProcess, G4Proton::Proton());	
    protonInelasticProcess->AddDataSet(flukaInelasticScatteringXS);
    //auto BGG = new G4BGGNucleonInelasticXS(G4Proton::Proton());
    //protonInelasticProcess->AddDataSet(BGG);

    protonInelasticProcess->RegisterMe(flukaModel);
    //auto theModel = new G4TheoFSGenerator("FTFP");
    //auto theStringModel = new G4FTFModel();
    //theStringModel->SetFragmentationModel(new G4ExcitedStringDecay());
    //theModel->SetHighEnergyGenerator(theStringModel);
    //theModel->SetQuasiElasticChannel(new G4QuasiElasticChannel());
    //auto theCascade = new G4GeneratorPrecompoundInterface();
    //theModel->SetTransport(theCascade);
    //theModel->SetMinEnergy(G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade());
    //theModel->SetMaxEnergy(G4HadronicParameters::Instance()->GetMaxEnergy());
    //protonInelasticProcess->RegisterMe(theModel);
    */
	       

  // NEUTRON
  const auto neutron = G4Neutron::Neutron();

  // NEUTRON INELASTIC
  const auto neutronInelasticProcess = new G4HadronInelasticProcess("neutronInelastic", neutron);
  // NB: No XS is set by default in the G4HadronInelasticProcess constructor.
  helper->RegisterProcess(neutronInelasticProcess, neutron);

  // Also non-HP: FLUKA neutron inelastic
  // IMPORTANT NB: Since flukaInelasticScatteringXS is SetForAllAtomsAndEnergies,
  // it needs to be set first (would erase any previously defined dataset, 
  // see G4CrossSectionDataStore::AddDataSet).
  neutronInelasticProcess->AddDataSet(flukaInelasticScatteringXS);
  const auto flukaNeutronModel = new FLUKANuclearInelasticModel();
  flukaNeutronModel->SetMinEnergy(fNeutronFLUKAMinE);
  neutronInelasticProcess->RegisterMe(flukaNeutronModel);

  // HP only: G4 neutron HP inelastic
  const auto neutronHPInelasticXS = new G4ParticleHPInelasticData(neutron);
  neutronInelasticProcess->AddDataSet(neutronHPInelasticXS);
  const auto neutronHPInelasticModel = new G4ParticleHPInelastic(neutron, "NeutronHPInelastic");
  neutronHPInelasticModel->SetMaxEnergy(fNeutronHPMaxE);
  neutronInelasticProcess->RegisterMe(neutronHPInelasticModel);

	
  // TO DO: Not elegant to have G4 neutron capture and fission included in FLUKA inelastic. 
  // Create a Physics constructor just for it? 
  // (NB: CANNOT use the neutron builders, 
  // because they would ALSO create a G4HadronInelasticProcess, while we use the FLUKA one).


  // NEUTRON CAPTURE
  const auto neutronCaptureProcess = new G4NeutronCaptureProcess();
  // NB: XS (G4NeutronCaptureXS) is already added, in G4NeutronCaptureProcess constructor.
  helper->RegisterProcess(neutronCaptureProcess, neutron);

  // Also non-HP: neutron rad capture
  const auto neutronRadCaptureModel = new G4NeutronRadCapture();
  neutronRadCaptureModel->SetMinEnergy(fNeutronFLUKAMinE);  // HP only
  neutronCaptureProcess->RegisterMe(neutronRadCaptureModel);

  // HP only: neutron HP capture
  const auto neutronHPCaptureXS = new G4ParticleHPCaptureData();
  neutronCaptureProcess->AddDataSet(neutronHPCaptureXS);
  const auto neutronHPCaptureModel = new G4ParticleHPCapture();
  neutronHPCaptureModel->SetMaxEnergy(fNeutronHPMaxE);
  neutronCaptureProcess->RegisterMe(neutronHPCaptureModel);
		
    
  // NEUTRON FISSION
  const auto neutronFissionProcess = new G4NeutronFissionProcess();
  // NB: XS (G4ZeroXS) is already added, in the process's constructor.
  helper->RegisterProcess(neutronFissionProcess, neutron);

  // HP only: neutron fission
  const auto neutronLEPFissionModel = new G4LFission();
  neutronLEPFissionModel->SetMinEnergy(fNeutronFLUKAMinE);
  neutronLEPFissionModel->SetMaxEnergy(G4HadronicParameters::Instance()->GetMaxEnergy());
  neutronFissionProcess->RegisterMe(neutronLEPFissionModel);

  // HP only: neutron HP fission	
  const auto neutronHPFissionXS = new G4ParticleHPFissionData();
  neutronFissionProcess->AddDataSet(neutronHPFissionXS);	
  const auto neutronHPFissionModel = new G4ParticleHPFission();
  neutronHPFissionModel->SetMaxEnergy(fNeutronHPMaxE);
  neutronFissionProcess->RegisterMe(neutronHPFissionModel);


  // PI+, PI-
  build_G4_process_helpers::buildInelasticProcess(G4PionPlus::PionPlus(),
                                                  helper,
                                                  flukaInelasticScatteringXS,
                                                  flukaModel);
  build_G4_process_helpers::buildInelasticProcess(G4PionMinus::PionMinus(),
                                                  helper,
                                                  flukaInelasticScatteringXS,
                                                  flukaModel);


  // KAONS
  build_G4_process_helpers::buildInelasticProcessForEachParticle(G4HadParticles::GetKaons(), 
                                                                 helper,
                                                                 flukaInelasticScatteringXS,
                                                                 flukaModel);


  const auto param = G4HadronicParameters::Instance();
  if (param->GetMaxEnergy() > param->EnergyThresholdForHeavyHadrons()) {

    // HYPERONS, ANTI-HYPERONS
    build_G4_process_helpers::buildInelasticProcessForEachParticle(G4HadParticles::GetHyperons(),
                                                                   helper,
                                                                   flukaInelasticScatteringXS,
                                                                   flukaModel);
    build_G4_process_helpers::buildInelasticProcessForEachParticle(G4HadParticles::GetAntiHyperons(),
                                                                   helper,
                                                                   flukaInelasticScatteringXS,
                                                                   flukaModel);

    // LIGHT ANTI-IONS: PBAR, NBAR, ANTI LIGHT IONS
    build_G4_process_helpers::buildInelasticProcessForEachParticle(G4HadParticles::GetLightAntiIons(),
                                                                   helper,
                                                                   flukaInelasticScatteringXS,
                                                                   flukaModel);

    // B-, C- BARYONS AND MESONS
    if (G4HadronicParameters::Instance()->EnableBCParticles() ) {
      build_G4_process_helpers::buildInelasticProcessForEachParticle(G4HadParticles::GetBCHadrons(),
                                                                     helper,
                                                                     flukaInelasticScatteringXS,
                                                                     flukaModel);
    }
  }
}


#endif // G4_USE_FLUKA
