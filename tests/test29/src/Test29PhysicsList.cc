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
//    ---------- fake Test29PhysicsList class -------
//    Created by Mikhail Kossov, 7 Dec 2004 
//
// **********************************************************************

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"

#include "globals.hh"

#include "Test29PhysicsList.hh"

Test29PhysicsList::Test29PhysicsList():  G4VUserPhysicsList()
{
  // Default cut values
  defaultCutValue = 2.0*mm;
  cutForGamma     = 1.0*micrometer;
  cutForElectron  = 1.0*micrometer;
  cutForProton    = 1.0*micrometer;

  SetVerboseLevel(1);
}

Test29PhysicsList::~Test29PhysicsList()
{}

void Test29PhysicsList::ConstructParticle()
{
  // Here are constructed all particles
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructAllShortLiveds();
}

// In this method, static member functions should be called for ALL particles to be used.

void Test29PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}
void Test29PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}
void Test29PhysicsList::ConstructMesons()
{
}
void Test29PhysicsList::ConstructBaryons()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}
void Test29PhysicsList::ConstructAllShortLiveds()
{
}
void Test29PhysicsList::ConstructProcess()
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
#include "G4MultipleScattering.hh"
#include "G4hLowEnergyIonisation.hh"

void Test29PhysicsList::ConstructEM()
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
	  pmanager->AddProcess(new G4MultipleScattering(),-1, 1,1);
	  pmanager->AddProcess(new G4eIonisation(),       -1, 2,2);
	  pmanager->AddProcess(new G4eBremsstrahlung(),   -1,-1,3);
	}
      else if (particleName == "e+")
	{
	  //positron
	  pmanager->AddProcess(new G4MultipleScattering(),-1, 1,1);
	  pmanager->AddProcess(new G4eIonisation(),       -1, 2,2);
	  pmanager->AddProcess(new G4eBremsstrahlung(),   -1,-1,3);
	  pmanager->AddProcess(new G4eplusAnnihilation(),  0,-1,4);
	}
      else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) &&
	       (particle->GetParticleName() != "chargedgeantino"))
	{
	  //all others charged particles except geantino
	  pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);

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

void Test29PhysicsList::ConstructGeneral()
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

void Test29PhysicsList::SetCuts()
{
  // defaultCutValue you have typed in is used

  if (verboseLevel >1){
    G4cout << "Test29PhysicsList::SetCuts:" << G4endl;
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

void Test29PhysicsList::SetCutForGamma(G4double cut)
{
  ResetCuts();
  cutForGamma = cut;
}

void Test29PhysicsList::SetCutForElectron(G4double cut)
{
  ResetCuts();
  cutForElectron = cut;
}

void Test29PhysicsList::SetCutForProton(G4double cut)
{
  ResetCuts();
  cutForProton = cut;
}

G4double Test29PhysicsList::GetCutForGamma() const
{
  return cutForGamma;
}

G4double Test29PhysicsList::GetCutForElectron() const
{
  return cutForElectron;
}

G4double Test29PhysicsList::GetCutForProton() const
{
  return cutForProton;
}





