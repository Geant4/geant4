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
/// \file field/field02/src/F02PhysicsList.cc
/// \brief Implementation of the F02PhysicsList class
//
//
// $Id: F02PhysicsList.cc 77794 2013-11-28 09:54:17Z gcosmo $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Timer.hh"
 
#include "F02PhysicsList.hh"
#include "F02DetectorConstruction.hh"
#include "F02PhysicsListMessenger.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4EnergyLossTables.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F02PhysicsList::F02PhysicsList(F02DetectorConstruction* p)
:  G4VUserPhysicsList(), fMaxChargedStep(DBL_MAX)
{
  fDet = p;

  defaultCutValue = 1.000*mm;

  fCutForGamma = defaultCutValue;
  fCutForElectron = defaultCutValue;

  SetVerboseLevel(1);
  fPhysicsListMessenger = new F02PhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F02PhysicsList::~F02PhysicsList()
{
  delete fPhysicsListMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBarions();

  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02PhysicsList::ConstructBosons()
{
  // gamma

  G4Gamma::GammaDefinition();

  // charged geantino

  G4ChargedGeantino::ChargedGeantinoDefinition();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02PhysicsList::ConstructLeptons()
{
  // leptons

  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02PhysicsList::ConstructMesons()
{
  //  mesons

  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02PhysicsList::ConstructBarions()
{
  //  barions

  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02PhysicsList::ConstructProcess()
{
  AddTransportation();

  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "F02StepCut.hh"

void F02PhysicsList::ConstructEM()
{
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma")
    {
      // Construct processes for gamma

      G4PhotoElectricEffect* fPhotoElectricEffect = new G4PhotoElectricEffect();
      G4ComptonScattering* fComptonScattering   = new G4ComptonScattering();
      G4GammaConversion* fGammaConversion     = new G4GammaConversion();

      pmanager->AddDiscreteProcess(fPhotoElectricEffect);
      pmanager->AddDiscreteProcess(fComptonScattering);
      pmanager->AddDiscreteProcess(fGammaConversion);
    }
    else if (particleName == "e-")
    {
      // Construct processes for electron

      G4eIonisation* feminusIonisation = new G4eIonisation();
      G4eBremsstrahlung* feminusBremsstrahlung = new G4eBremsstrahlung();
      F02StepCut* feminusStepCut = new F02StepCut();
      feminusStepCut->SetMaxStep(fMaxChargedStep);

      pmanager->AddProcess(feminusIonisation,-1,2,2);
      pmanager->AddProcess(feminusBremsstrahlung,-1,-1,3);
      pmanager->AddProcess(feminusStepCut,-1,-1,4);
    }
    else if (particleName == "e+")
    {
      // Construct processes for positron

      G4eIonisation* feplusIonisation = new G4eIonisation();
      G4eBremsstrahlung* feplusBremsstrahlung = new G4eBremsstrahlung();
      F02StepCut* feplusStepCut = new F02StepCut();
      feplusStepCut->SetMaxStep(fMaxChargedStep);

      pmanager->AddProcess(feplusIonisation,-1,2,2);
      pmanager->AddProcess(feplusBremsstrahlung,-1,-1,3);
      pmanager->AddProcess(feplusStepCut,-1,-1,5);
    }
    else if( particleName == "mu+" || particleName == "mu-" )
    {
      // Construct processes for muon+

      F02StepCut* muonStepCut = new F02StepCut();
      muonStepCut->SetMaxStep(fMaxChargedStep);
      G4MuIonisation* muIonisation = new G4MuIonisation();

      pmanager->AddProcess(new G4MuMultipleScattering(),-1,1,1);
      pmanager->AddProcess(muIonisation,-1,2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1,3);
      pmanager->AddProcess(new G4MuPairProduction(),-1,-1,4);
      pmanager->AddProcess(muonStepCut,-1,-1,3);
    }
    else if ( particleName == "proton"
              || particleName == "antiproton"
              || particleName == "pi+"
              || particleName == "pi-"
              || particleName == "kaon+"
              || particleName == "kaon-"
            )
    {
      F02StepCut* theHadronStepCut = new F02StepCut();
      theHadronStepCut->SetMaxStep(fMaxChargedStep);

      G4hIonisation* thehIonisation = new G4hIonisation();
      G4hMultipleScattering* thehMultipleScattering =
                                                  new G4hMultipleScattering();

      pmanager->AddProcess(thehMultipleScattering,-1,1,1);
      pmanager->AddProcess(thehIonisation,-1,2,2);
      pmanager->AddProcess(theHadronStepCut,-1,-1,3);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"

void F02PhysicsList::ConstructGeneral()
{
  // Add Decay Process

  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (theDecayProcess->IsApplicable(*particle))
    {
      pmanager ->AddProcess(theDecayProcess);

      // set ordering for PostStepDoIt and AtRestDoIt

      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02PhysicsList::SetCuts()
{
  G4Timer theTimer;
  theTimer.Start();
  if (verboseLevel >0)
  {
    G4cout << "F02PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
 
  SetCutValue(fCutForGamma,"gamma");

  SetCutValue(fCutForElectron,"e-");
  SetCutValue(fCutForElectron,"e+");

  if (verboseLevel>1)     DumpCutValuesTable();

  theTimer.Stop();
  G4cout.precision(6);
  G4cout << G4endl;
  G4cout << "total time(SetCuts)=" << theTimer.GetUserElapsed()
                                   << " s " <<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02PhysicsList::SetGammaCut(G4double val)
{
  fCutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02PhysicsList::SetElectronCut(G4double val)
{
  fCutForElectron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02PhysicsList::SetMaxStep(G4double step)
{
  fMaxChargedStep = step;
  G4cout << " MaxChargedStep=" << fMaxChargedStep << G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
