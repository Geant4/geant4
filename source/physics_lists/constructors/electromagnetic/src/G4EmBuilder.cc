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
// Geant4 class G4EmBuilder
//
// Author V.Ivanchenko 22.05.2020
//

#include "G4EmBuilder.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysListUtil.hh"
#include "G4EmParameters.hh"
#include "G4VEnergyLossProcess.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4CoulombScattering.hh"
#include "G4WentzelVIModel.hh"

#include "G4ProcessManager.hh"
#include "G4TransportationWithMsc.hh"
#include "G4TransportationProcessType.hh"

#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4hBremsstrahlungModel.hh"
#include "G4hPairProductionModel.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4NuclearStopping.hh"

#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4CoulombScattering.hh"

#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Lambda.hh"
#include "G4AntiLambda.hh"

#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

#include "G4PhysicsListHelper.hh"
#include "G4HadParticles.hh"
#include "G4HadronicParameters.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

void G4EmBuilder::ConstructBasicEmPhysics(G4hMultipleScattering* hmsc,
				          const std::vector<G4int>& partList)
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();

  for( auto & pdg : partList ) {
    auto part = table->FindParticle( pdg );
    if ( part == nullptr || part->GetPDGCharge() == 0.0 ) { continue; }
    ph->RegisterProcess(hmsc, part);
    ph->RegisterProcess(new G4hIonisation(), part);
  }
}

void G4EmBuilder::ConstructIonEmPhysics(G4hMultipleScattering* hmsc, 
                                        G4NuclearStopping* nucStopping)
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  G4ParticleDefinition* part = G4Deuteron::Deuteron();
  ph->RegisterProcess(hmsc, part);
  ph->RegisterProcess(new G4hIonisation(), part);

  part = G4Triton::Triton();
  ph->RegisterProcess(hmsc, part);
  ph->RegisterProcess(new G4hIonisation(), part);

  part = G4Alpha::Alpha();
  ph->RegisterProcess(new G4hMultipleScattering(), part);
  ph->RegisterProcess(new G4ionIonisation(), part);
  if( nucStopping != nullptr ) {
    ph->RegisterProcess(nucStopping, part);
  }

  part = G4He3::He3();
  ph->RegisterProcess(new G4hMultipleScattering(), part);
  ph->RegisterProcess(new G4ionIonisation(), part);
  if( nucStopping != nullptr ) {
    ph->RegisterProcess(nucStopping, part);
  }
}

void G4EmBuilder::ConstructIonEmPhysicsSS()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  G4ParticleDefinition* part = G4Deuteron::Deuteron();
  ph->RegisterProcess(new G4hIonisation(), part);
  ph->RegisterProcess(new G4CoulombScattering(), part);

  part = G4Triton::Triton();
  ph->RegisterProcess(new G4hIonisation(), part);
  ph->RegisterProcess(new G4CoulombScattering(), part);

  part = G4Alpha::Alpha();
  ph->RegisterProcess(new G4ionIonisation(), part);
  ph->RegisterProcess(new G4CoulombScattering(), part);

  part = G4He3::He3();
  ph->RegisterProcess(new G4ionIonisation(), part);
  ph->RegisterProcess(new G4CoulombScattering(), part);
}

void G4EmBuilder::ConstructLightHadrons(G4ParticleDefinition* part1,
			   	        G4ParticleDefinition* part2,
				        G4bool isHEP, G4bool isProton, 
                                        G4bool isWVI)
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  G4hMultipleScattering* msc = new G4hMultipleScattering();
  if(isWVI) { msc->SetEmModel(new G4WentzelVIModel()); }
  G4CoulombScattering* ss = ( isWVI ) ? new G4CoulombScattering() : nullptr;

  ph->RegisterProcess(msc, part1);
  ph->RegisterProcess(new G4hIonisation(), part1);

  G4hBremsstrahlung* brem = ( isHEP ) ? new G4hBremsstrahlung() : nullptr;
  G4hPairProduction* pair = ( isHEP ) ? new G4hPairProduction() : nullptr;

  if( isHEP ) {
    ph->RegisterProcess(brem, part1);
    ph->RegisterProcess(pair, part1);
  }
  if( isWVI ) { ph->RegisterProcess(ss, part1); }

  if( isProton ) {
    msc = new G4hMultipleScattering();
    if(isWVI) { 
      msc->SetEmModel(new G4WentzelVIModel());
      ss = new G4CoulombScattering();
    }
  }
  ph->RegisterProcess(msc, part2);
  ph->RegisterProcess(new G4hIonisation(), part2);
  if( isHEP ) {
    ph->RegisterProcess(brem, part2);
    ph->RegisterProcess(pair, part2);
  }
  if( isWVI ) { ph->RegisterProcess(ss, part2); }
}

void G4EmBuilder::ConstructLightHadronsSS(G4ParticleDefinition* part1,
			   	          G4ParticleDefinition* part2,
					  G4bool isHEP)
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  ph->RegisterProcess(new G4hIonisation(), part1);

  G4hBremsstrahlung* brem = ( isHEP ) ? new G4hBremsstrahlung() : nullptr;
  G4hPairProduction* pair = ( isHEP ) ? new G4hPairProduction() : nullptr;

  if( isHEP ) {
    ph->RegisterProcess(brem, part1);
    ph->RegisterProcess(pair, part1);
  }
  ph->RegisterProcess(new G4CoulombScattering(), part1);

  ph->RegisterProcess(new G4hIonisation(), part2);
  if( isHEP ) {
    ph->RegisterProcess(brem, part2);
    ph->RegisterProcess(pair, part2);
  }
  ph->RegisterProcess(new G4CoulombScattering(), part2);
}

void G4EmBuilder::ConstructCharged(G4hMultipleScattering* hmsc,
                                   G4NuclearStopping* nucStopping,
                                   G4bool isWVI) 
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4EmParameters* param = G4EmParameters::Instance();
  G4HadronicParameters* hpar = G4HadronicParameters::Instance();
  G4bool isHEP = ( param->MaxKinEnergy() > hpar->EnergyThresholdForHeavyHadrons() );

  // muon multiple and single scattering
  G4MuMultipleScattering* mumsc = new G4MuMultipleScattering();
  if(isWVI) { mumsc->SetEmModel(new G4WentzelVIModel()); }
  G4CoulombScattering* muss = ( isWVI ) ? new G4CoulombScattering() : nullptr;

  // Add standard EM Processes
  // mu+-
  G4ParticleDefinition* part = G4MuonPlus::MuonPlus();
  ph->RegisterProcess(mumsc, part);
  ph->RegisterProcess(new G4MuIonisation(), part);

  // muon bremsstrahlung and pair production
  G4MuBremsstrahlung* mub = ( isHEP ) ? new G4MuBremsstrahlung() : nullptr;
  G4MuPairProduction* mup = ( isHEP ) ? new G4MuPairProduction() : nullptr;

  if( isHEP ) {
    ph->RegisterProcess(mub, part);
    ph->RegisterProcess(mup, part);
  }
  if( isWVI ) { ph->RegisterProcess(muss, part); }

  part = G4MuonMinus::MuonMinus();
  ph->RegisterProcess(mumsc, part);
  ph->RegisterProcess(new G4MuIonisation(), part);
  if( isHEP ) {
    ph->RegisterProcess(mub, part);
    ph->RegisterProcess(mup, part);
  }
  if( isWVI ) { ph->RegisterProcess(muss, part); }

  // pi+-
  ConstructLightHadrons(G4PionPlus::PionPlus(), G4PionMinus::PionMinus(), isHEP, false, isWVI);

  // K+-
  ConstructLightHadrons(G4KaonPlus::KaonPlus(), G4KaonMinus::KaonMinus(), isHEP, false, isWVI);
 
  // p, pbar
  ConstructLightHadrons(G4Proton::Proton(), G4AntiProton::AntiProton(), isHEP, true, isWVI);
  if( nucStopping != nullptr ) {
    ph->RegisterProcess(nucStopping, G4Proton::Proton());
  }

  // ions
  ConstructIonEmPhysics(hmsc, nucStopping);

  // hyperons and anti particles
  if( isHEP ) { 
    ConstructBasicEmPhysics(hmsc, G4HadParticles::GetHeavyChargedParticles());

    // b- and c- charged particles
    if( hpar->EnableBCParticles() ) {
      ConstructBasicEmPhysics(hmsc, G4HadParticles::GetBCChargedHadrons());
    }
    // light hyper-nuclei
    if( hpar->EnableHyperNuclei() ) {
      ConstructBasicEmPhysics(hmsc, G4HadParticles::GetChargedHyperNuclei());
    }
  }
}

void G4EmBuilder::ConstructChargedSS(G4hMultipleScattering* hmsc) 
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4EmParameters* param = G4EmParameters::Instance();
  G4HadronicParameters* hpar = G4HadronicParameters::Instance();
  G4bool isHEP = ( param->MaxKinEnergy() > hpar->EnergyThresholdForHeavyHadrons() );

  // muon multiple and single scattering
  G4CoulombScattering* muss = new G4CoulombScattering();

  // Add standard EM Processes
  // mu+-
  G4ParticleDefinition* part = G4MuonPlus::MuonPlus();
  ph->RegisterProcess(new G4MuIonisation(), part);

  // muon bremsstrahlung and pair production
  G4MuBremsstrahlung* mub = ( isHEP ) ? new G4MuBremsstrahlung() : nullptr;
  G4MuPairProduction* mup = ( isHEP ) ? new G4MuPairProduction() : nullptr;

  if( isHEP ) {
    ph->RegisterProcess(mub, part);
    ph->RegisterProcess(mup, part);
  }
  ph->RegisterProcess(muss, part);

  part = G4MuonMinus::MuonMinus();
  ph->RegisterProcess(new G4MuIonisation(), part);
  if( isHEP ) {
    ph->RegisterProcess(mub, part);
    ph->RegisterProcess(mup, part);
  }
  ph->RegisterProcess(muss, part);

  // pi+-
  ConstructLightHadronsSS(G4PionPlus::PionPlus(), G4PionMinus::PionMinus(), isHEP);

  // K+-
  ConstructLightHadronsSS(G4KaonPlus::KaonPlus(), G4KaonMinus::KaonMinus(), isHEP);
 
  // p, pbar
  ConstructLightHadronsSS(G4Proton::Proton(), G4AntiProton::AntiProton(), isHEP);
  // ions
  ConstructIonEmPhysicsSS();

  // hyperons and anti particles
  if( isHEP ) { 
    ConstructBasicEmPhysics(hmsc, G4HadParticles::GetHeavyChargedParticles());

    // b- and c- charged particles
    if( hpar->EnableBCParticles() ) {
      ConstructBasicEmPhysics(hmsc, G4HadParticles::GetBCChargedHadrons());
    }
    // light hyper-nuclei
    if( hpar->EnableHyperNuclei() ) {
      ConstructBasicEmPhysics(hmsc, G4HadParticles::GetChargedHyperNuclei());
    }
  }
}

void G4EmBuilder::ConstructMinimalEmSet()
{
  // instantiate singletones for physics
  G4PhysListUtil::InitialiseParameters();
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  G4NeutrinoMu::NeutrinoMu();
  G4AntiNeutrinoMu::AntiNeutrinoMu();
  G4NeutrinoE::NeutrinoE();
  G4AntiNeutrinoE::AntiNeutrinoE();
  // gamma
  G4Gamma::Gamma();
  // leptons
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();
  // mesons
  G4PionPlus::PionPlus();
  G4PionMinus::PionMinus();
  G4PionZero::PionZero();
  G4KaonPlus::KaonPlus();
  G4KaonMinus::KaonMinus();
  // barions
  G4Proton::Proton();
  G4AntiProton::AntiProton();
  G4Neutron::Neutron();
  G4AntiNeutron::AntiNeutron();
  G4Lambda::Lambda();
  G4AntiLambda::AntiLambda();
  // ions
  G4Deuteron::Deuteron();
  G4Triton::Triton();
  G4He3::He3();
  G4Alpha::Alpha();
  G4GenericIon::GenericIon();
}

void G4EmBuilder::PrepareEMPhysics()
{
  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* ad = man->AtomDeexcitation();
  if(nullptr == ad) {
    ad = new G4UAtomicDeexcitation();
    man->SetAtomDeexcitation(ad);
  }
}

void G4EmBuilder::ConstructElectronMscProcess(G4VMscModel* msc1, G4VMscModel* msc2,
                                              G4ParticleDefinition* particle)
{
  G4TransportationWithMscType type =
    G4EmParameters::Instance()->TransportationWithMsc();
  G4ProcessManager* procManager = particle->GetProcessManager();
  auto plist = procManager->GetProcessList();
  G4int ptype = (0 < plist->size()) ? (*plist)[0]->GetProcessSubType() : 0;
  if(type != G4TransportationWithMscType::fDisabled && ptype == TRANSPORTATION) {
    // Remove default G4Transportation and replace with G4TransportationWithMsc.
    procManager->RemoveProcess(0);
    G4TransportationWithMsc* transportWithMsc = new G4TransportationWithMsc(
      G4TransportationWithMsc::ScatteringType::MultipleScattering);
    if(type == G4TransportationWithMscType::fMultipleSteps) {
      transportWithMsc->SetMultipleSteps(true);
    }
    transportWithMsc->AddMscModel(msc1);
    if(msc2 != nullptr) {
      transportWithMsc->AddMscModel(msc2);
    }
    procManager->AddProcess(transportWithMsc, -1, 0, 0);
  } else {
    // Register as a separate process.
    G4eMultipleScattering* msc = new G4eMultipleScattering;
    msc->SetEmModel(msc1);
    if(msc2 != nullptr) {
      msc->SetEmModel(msc2);
    }
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(msc, particle);
  }
}
