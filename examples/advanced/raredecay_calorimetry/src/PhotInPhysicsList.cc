//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: PhotInPhysicsList.cc,v 1.2 2005-05-31 15:23:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#define debug

#include "PhotInPhysicsList.hh"

PhotInPhysicsList::PhotInPhysicsList():  G4VUserPhysicsList() { SetVerboseLevel(1); }

PhotInPhysicsList::~PhotInPhysicsList() {}

void PhotInPhysicsList::ConstructParticle()
{
  // In this method, static member functions for particles should be called
  // for all particles which user is going to use in the simulation. If not
  // defined particle appear in the simulation? it can cause a WORNING which
  // means that in the simulation appeard unexpected particles. Then add them.

  // @@ Word "Definition" can be skipped. - Old fashion (M.K.)

  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gammas
  G4Gamma::GammaDefinition();

  // leptons (without tau and it's neutrino)
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  //  mesons
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

  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

  // hyperons
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
}


void PhotInPhysicsList::ConstructProcess()
{
#ifdef debug
  G4cout<<"PhotInPhysicsList::ConstructProcess: is called "<<G4endl;
#endif
  AddTransportation();     // Transportation is a "process" and defined in the basic class

  // Add Electromagnetic interaction Processes and Decays
  G4Decay* theDecayProcess = new G4Decay(); // @@ When this class is decayed? (M.K.)
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
#ifdef debug
				G4cout<<"PhotInPhysList::ConstructProcess: Part="<<particle->GetParticleName()<<G4endl;
#endif
    G4ProcessManager* pmanager = particle->GetProcessManager();
    // Decays
    if (theDecayProcess->IsApplicable(*particle))
    { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }

    G4String particleName = particle->GetParticleName();
    // EM Interactions     
    if (particleName == "gamma")    // gamma
    {
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } 
    else if (particleName == "e-")    //electron
    {
      G4VProcess* theeminusMultipleScattering = new G4MultipleScattering();
      G4VProcess* theeminusIonisation         = new G4eIonisation();
      G4VProcess* theeminusBremsstrahlung     = new G4eBremsstrahlung();
      //
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);
      //      
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         idxAlongStep,2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung,     idxAlongStep,3);      
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         idxPostStep,2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung,     idxPostStep,3);

    }
    else if (particleName == "e+")    //positron
    {
      G4VProcess* theeplusMultipleScattering = new G4MultipleScattering();
      G4VProcess* theeplusIonisation         = new G4eIonisation();
      G4VProcess* theeplusBremsstrahlung     = new G4eBremsstrahlung();
      G4VProcess* theeplusAnnihilation       = new G4eplusAnnihilation();
      //
      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);
      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);
      //
      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation,         idxAlongStep,2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung,     idxAlongStep,3);      
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation,         idxPostStep,2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung,     idxPostStep,3);
      pmanager->SetProcessOrdering(theeplusAnnihilation,       idxPostStep,4);
    }
    else if( particleName == "mu+" || particleName == "mu-"    )    //muon  of both signs
    {
      G4VProcess* aMultipleScattering = new G4MultipleScattering();
      G4VProcess* aBremsstrahlung     = new G4MuBremsstrahlung();
      G4VProcess* aPairProduction     = new G4MuPairProduction();
      G4VProcess* anIonisation        = new G4MuIonisation();
      //
      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);
      pmanager->AddProcess(aBremsstrahlung);
      pmanager->AddProcess(aPairProduction);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxAlongStep,2);
      pmanager->SetProcessOrdering(aBremsstrahlung,     idxAlongStep,3);
      pmanager->SetProcessOrdering(aPairProduction,     idxAlongStep,4);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxPostStep,2);
      pmanager->SetProcessOrdering(aBremsstrahlung,     idxPostStep,3);
      pmanager->SetProcessOrdering(aPairProduction,     idxPostStep,4);
    }
    else if(!particle->IsShortLived() && particle->GetPDGCharge() && 
	            particle->GetParticleName()!="chargedgeantino")// all others charged particles
    {
      G4VProcess* aMultipleScattering = new G4MultipleScattering();
      G4VProcess* anIonisation        = new G4hIonisation();
      //
      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxAlongStep,2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxPostStep,2);
    }
  }
}

void PhotInPhysicsList::SetCuts()
{
  if(verboseLevel>0) G4cout<<"PhotInPhysicsList::SetCuts: default cut length : "
                           <<G4BestUnit(defaultCutValue,"Length")<<G4endl;
  // These values are used as the default production thresholds for the world volume.
  SetCutsWithDefault();

  // Production thresholds for detector regions

  G4double fact = 1.; // Multiplicative factor for default cuts
  for(G4int i=0; i< PhotInNumSections; i++)
  { 
    G4Region* reg = G4RegionStore::GetInstance()-> GetRegion(PhotInRegName[i]);
    G4ProductionCuts* cuts = new G4ProductionCuts;
    cuts->SetProductionCut(defaultCutValue*fact);
    reg->SetProductionCuts(cuts);
    fact *= 10.; // @@ Increment the multiplicative factor by order of magnitude
  }
}


