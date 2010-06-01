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
//    ---------- fake Test19PhysicsList class -------
//    Created by Mikhail Kossov, 7 Dec 2004 
//
// **********************************************************************

#include "Test19PhysicsList.hh"

Test19PhysicsList::Test19PhysicsList():  G4VUserPhysicsList()
{
  // Default cut values
  defaultCutValue = 2.0*mm;
  cutForGamma     = 1.0*micrometer;
  cutForElectron  = 1.0*micrometer;
  cutForProton    = 1.0*micrometer;

  SetVerboseLevel(1);
}

Test19PhysicsList::~Test19PhysicsList()
{}

void Test19PhysicsList::ConstructParticle()
{
  // Here are constructed all particles
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructIons();
  ConstructAllShortLiveds();
}

// In this method, static member functions should be called for ALL particles to be used.

void Test19PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}
void Test19PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4TauPlus::TauPlusDefinition();
  G4TauMinus::TauMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();  
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();  
}
void Test19PhysicsList::ConstructMesons()
{
  // mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();

  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();

  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();

  G4DMesonPlus::DMesonPlusDefinition();
  G4DMesonMinus::DMesonMinusDefinition();
  G4DMesonZero::DMesonZeroDefinition();
  G4AntiDMesonZero::AntiDMesonZeroDefinition();
  G4DsMesonPlus::DsMesonPlusDefinition();
  G4DsMesonMinus::DsMesonMinusDefinition();
  G4JPsi::JPsiDefinition();

  G4BMesonPlus::BMesonPlusDefinition();
  G4BMesonMinus::BMesonMinusDefinition();
  G4BMesonZero::BMesonZeroDefinition();
  G4AntiBMesonZero::AntiBMesonZeroDefinition();
  G4BsMesonZero::BsMesonZeroDefinition();
  G4AntiBsMesonZero::AntiBsMesonZeroDefinition();
}
void Test19PhysicsList::ConstructBaryons()
{
  // barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

  G4Lambda::LambdaDefinition();
  G4SigmaPlus::SigmaPlusDefinition();
  G4SigmaZero::SigmaZeroDefinition();
  G4SigmaMinus::SigmaMinusDefinition();
  G4XiMinus::XiMinusDefinition();
  G4XiZero::XiZeroDefinition();
  G4OmegaMinus::OmegaMinusDefinition();

  G4AntiLambda::AntiLambdaDefinition();
  G4AntiSigmaPlus::AntiSigmaPlusDefinition();
  G4AntiSigmaZero::AntiSigmaZeroDefinition();
  G4AntiSigmaMinus::AntiSigmaMinusDefinition();
  G4AntiXiMinus::AntiXiMinusDefinition();
  G4AntiXiZero::AntiXiZeroDefinition();
  G4AntiOmegaMinus::AntiOmegaMinusDefinition();

  G4LambdacPlus::LambdacPlusDefinition();
  G4SigmacPlusPlus::SigmacPlusPlusDefinition();
  G4SigmacPlus::SigmacPlusDefinition();
  G4SigmacZero::SigmacZeroDefinition();
  G4XicPlus::XicPlusDefinition();
  G4XicZero::XicZeroDefinition();
  G4OmegacZero::OmegacZeroDefinition();

  G4AntiLambdacPlus::AntiLambdacPlusDefinition();
  G4AntiSigmacPlusPlus::AntiSigmacPlusPlusDefinition();
  G4AntiSigmacPlus::AntiSigmacPlusDefinition();
  G4AntiSigmacZero::AntiSigmacZeroDefinition();
  G4AntiXicPlus::AntiXicPlusDefinition();
  G4AntiXicZero::AntiXicZeroDefinition();
  G4AntiOmegacZero::AntiOmegacZeroDefinition();
}
void Test19PhysicsList::ConstructIons()
{
  // ions
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4He3::He3Definition();
  G4Alpha::AlphaDefinition();
  G4ParticleDefinition* genericIon = G4GenericIon::GenericIonDefinition();
  genericIon->SetProcessManager(new G4ProcessManager(genericIon));

  //  Construct light ions @@ Is that the same as above?
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}
void Test19PhysicsList::ConstructAllShortLiveds()
{
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  

}
void Test19PhysicsList::ConstructProcess()
{
  // Transportation, electromagnetic and general processes 
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

// Here are respective header files for chosen processes

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4hLowEnergyIonisation.hh"

void Test19PhysicsList::ConstructEM()
{
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "gamma")
	{
	  //gamma
	  pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
	  pmanager->AddDiscreteProcess(new G4ComptonScattering());
	  pmanager->AddDiscreteProcess(new G4GammaConversion());
	}
      else if (particleName == "e-")
	{
	  //electron
	  pmanager->AddProcess(new G4eIonisation(),       -1, 2,2);
	  pmanager->AddProcess(new G4eBremsstrahlung(),   -1,-1,3);
	}
      else if (particleName == "e+")
	{
	  //positron
	  pmanager->AddProcess(new G4eIonisation(),       -1, 2,2);
	  pmanager->AddProcess(new G4eBremsstrahlung(),   -1,-1,3);
	  pmanager->AddProcess(new G4eplusAnnihilation(),  0,-1,4);
	}
      else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) &&
	       (particle->GetParticleName() != "chargedgeantino"))
	{
	  //all others charged particles except geantino

	  G4double demax = 0.05;  // try to lose at most 5% of the energy in
	  //    a single step (in limit of large energies)
	  G4double stmin = 1.e-9 * m;  // length of the final step: 10 angstrom
	  // reproduced angular distribution of TRIM

	  G4hLowEnergyIonisation* lowEIonisation = new G4hLowEnergyIonisation();
	  pmanager->AddProcess( lowEIonisation, -1,2,2);
	  lowEIonisation->SetStepFunction( demax, stmin );
	}
    }
}

#include "G4Decay.hh"

void Test19PhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) {
      pmanager ->AddProcess(theDecayProcess);
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

void Test19PhysicsList::SetCuts()
{
  // defaultCutValue you have typed in is used

  if (verboseLevel >1){
    G4cout << "Test19PhysicsList::SetCuts:" << G4endl;
  }

  // set cut values for gamma at first and for e- second
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForElectron, "e+");

  // set cut values for proton
  SetCutValue(cutForProton, "proton"); 
  SetCutValue(cutForProton, "anti_proton");

  // SetCutValueForOthers(defaultCutValue); 
 
  if (verboseLevel >1) { 
    DumpCutValuesTable(); 
  } 
}

void Test19PhysicsList::SetCutForGamma(G4double cut)
{
  ResetCuts();
  cutForGamma = cut;
}

void Test19PhysicsList::SetCutForElectron(G4double cut)
{
  ResetCuts();
  cutForElectron = cut;
}

void Test19PhysicsList::SetCutForProton(G4double cut)
{
  ResetCuts();
  cutForProton = cut;
}

G4double Test19PhysicsList::GetCutForGamma() const
{
  return cutForGamma;
}

G4double Test19PhysicsList::GetCutForElectron() const
{
  return cutForElectron;
}

G4double Test19PhysicsList::GetCutForProton() const
{
  return cutForProton;
}





