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
// $Id: Em8PhysicsList.cc,v 1.10 2003-08-28 09:29:37 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 

#include "G4Timer.hh"
   
#include "Em8PhysicsList.hh"
#include "Em8DetectorConstruction.hh"
#include "Em8PhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4EnergyLossTables.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include <iomanip>

#include "G4FastSimulationManagerProcess.hh"


/////////////////////////////////////////////////////////////
//
//

Em8PhysicsList::Em8PhysicsList(Em8DetectorConstruction* p)
:  G4VUserPhysicsList(),
 MaxChargedStep(DBL_MAX),
 thePhotoElectricEffect(NULL),theComptonScattering(NULL),
 theGammaConversion(NULL),
 theeminusMultipleScattering(NULL),theeminusIonisation(NULL),
 theeminusBremsstrahlung(NULL),
 theeplusMultipleScattering(NULL),theeplusIonisation(NULL),
 theeplusBremsstrahlung(NULL),
 theeplusAnnihilation(NULL),
 theeminusStepCut(NULL),theeplusStepCut(NULL)
{
  pDet = p;

  defaultCutValue = 1.000*mm ;

  cutForGamma = defaultCutValue ;
  cutForElectron = defaultCutValue ;
  cutForElectron = 10*mm;

  SetVerboseLevel(1);
  physicsListMessenger = new Em8PhysicsListMessenger(this);
}

/////////////////////////////////////////////////////////////////////////
//
//

Em8PhysicsList::~Em8PhysicsList()
{
  delete physicsListMessenger; 
}

///////////////////////////////////////////////////////////////////////////
//
//

void Em8PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBarions();
}

////////////////////////////////////////////////////////////////////////////
//
//

void Em8PhysicsList::ConstructBosons()
{
  // gamma
  G4Gamma::GammaDefinition();
}

void Em8PhysicsList::ConstructLeptons()
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

void Em8PhysicsList::ConstructMesons()
{
 //  mesons

  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
}


void Em8PhysicsList::ConstructBarions()
{
//  barions

  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
}


///////////////////////////////////////////////////////////////////////
//
//

void Em8PhysicsList::ConstructProcess()
{
  AddTransportation();
  AddParameterisation();

  ConstructEM();
  ConstructGeneral();
}

/////////////////////////////////////////////////////////////////////////////
//
//

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering52.hh"

#include "G4eIonisation52.hh"
#include "G4eBremsstrahlung52.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation52.hh"
#include "G4MuBremsstrahlung52.hh"
#include "G4MuPairProduction52.hh"

#include "G4hIonisation52.hh"
#include "G4PAIonisation.hh"
#include "G4ForwardXrayTR.hh"

#include "Em8StepCut.hh"

#include "G4IonisationByLogicalVolume.hh"

void Em8PhysicsList::ConstructEM()
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

      thePhotoElectricEffect = new G4PhotoElectricEffect();
      theComptonScattering   = new G4ComptonScattering();
      theGammaConversion     = new G4GammaConversion();

      pmanager->AddDiscreteProcess(thePhotoElectricEffect);
      pmanager->AddDiscreteProcess(theComptonScattering);

      pmanager->AddDiscreteProcess(theGammaConversion);

    }
    else if (particleName == "e-")
    {
      // Construct processes for electron

      // theeminusMultipleScattering = new G4MultipleScattering();
      //   theeminusIonisation = new G4eIonisation();
      theeminusBremsstrahlung = new G4eBremsstrahlung52();

     //   fPAIonisation = new G4PAIonisation("Xenon") ;
     // fForwardXrayTR = new G4ForwardXrayTR("Air","Polypropelene","XrayTR") ;

      theeminusStepCut = new Em8StepCut();

      //  pmanager->AddProcess(theeminusMultipleScattering,-1,1,1);

      //   pmanager->AddProcess(theeminusIonisation,-1,2,2);

       pmanager->AddProcess(new G4IonisationByLogicalVolume(particleName,
                                     pDet->GetLogicalAbsorber(),
                                    "IonisationByLogVol"),-1,2,-2);

       pmanager->AddProcess(theeminusBremsstrahlung,-1,-1,3);

       //   pmanager->AddProcess(fPAIonisation,-1,2,2);

       //  pmanager->AddProcess(fForwardXrayTR,-1,-1,2);

      pmanager->AddProcess(theeminusStepCut,-1,-1,4);
      theeminusStepCut->SetMaxStep(MaxChargedStep) ;

    }
    else if (particleName == "e+")
    {
      // Construct processes for positron

      theeplusIonisation = new G4eIonisation52();
      theeplusBremsstrahlung = new G4eBremsstrahlung52();
      // theeplusAnnihilation = new G4eplusAnnihilation();

      //  fPAIonisation = new G4PAIonisation("Xenon") ;
      // fForwardXrayTR = new G4ForwardXrayTR("Air","Polypropelene","XrayTR") ;

      theeplusStepCut = new Em8StepCut();

      //  pmanager->AddProcess(theeplusMultipleScattering,-1,1,1);
      pmanager->AddProcess(theeplusIonisation,-1,2,2);
      pmanager->AddProcess(theeplusBremsstrahlung,-1,-1,3);
      //  pmanager->AddProcess(theeplusAnnihilation,0,-1,4);

      //  pmanager->AddProcess(fPAIonisation,-1,2,2);


      // pmanager->AddProcess(fForwardXrayTR,-1,-1,2);

      pmanager->AddProcess(theeplusStepCut,-1,-1,5);
      theeplusStepCut->SetMaxStep(MaxChargedStep) ;

    }
    else if( particleName == "mu+" ||
               particleName == "mu-"    )
    {
     // Construct processes for muon+

      Em8StepCut* muonStepCut = new Em8StepCut();
       pmanager->AddProcess(new G4IonisationByLogicalVolume(particleName,
                                     pDet->GetLogicalAbsorber(),
                                    "IonisationByLogVol"),-1,2,-2);

      // G4MuIonisation* themuIonisation = new G4MuIonisation() ;
     pmanager->AddProcess(new G4MultipleScattering52(),-1,1,1);
     //  pmanager->AddProcess(themuIonisation,-1,2,2);
     pmanager->AddProcess(new G4MuBremsstrahlung52(),-1,-1,3);
     pmanager->AddProcess(new G4MuPairProduction52(),-1,-1,4);

      //  pmanager->AddProcess(new G4PAIonisation("Xenon"),-1,2,2) ;
      pmanager->AddProcess( muonStepCut,-1,-1,3);
     muonStepCut->SetMaxStep(MaxChargedStep) ;

    }
    else if (
                  particleName == "proton"
               || particleName == "antiproton"
               || particleName == "pi+"
               || particleName == "pi-"
               || particleName == "kaon+"
               || particleName == "kaon-"
            )
    {
        Em8StepCut* thehadronStepCut = new Em8StepCut();

      //  G4hIonisation* thehIonisation = new G4hIonisation() ;
      G4MultipleScattering52* thehMultipleScattering =
                     new G4MultipleScattering52() ;

        pmanager->AddProcess(new G4IonisationByLogicalVolume(particleName,
                                     pDet->GetLogicalAbsorber(),
                                    "IonisationByLogVolHadr"),-1,2,-2);

        pmanager->AddProcess(thehMultipleScattering,-1,1,1);
      //  pmanager->AddProcess(thehIonisation,-1,2,2);

      //  pmanager->AddProcess(new G4PAIonisation("Xenon"),-1,2,2) ;
      // pmanager->AddProcess(new G4PAIonisation("Argon"),-1,2,2) ;

        pmanager->AddProcess( thehadronStepCut,-1,-1,3);
        thehadronStepCut->SetMaxStep(MaxChargedStep) ;
      // thehadronStepCut->SetMaxStep(10*mm) ;

    }
  }
}


#include "G4Decay.hh"

void Em8PhysicsList::ConstructGeneral()
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

/////////////////////////////////////////////////////////////////////////////

void Em8PhysicsList::AddParameterisation()
{
  G4FastSimulationManagerProcess* theFastSimulationManagerProcess = 
                                  new G4FastSimulationManagerProcess() ;
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value() ;
    G4ProcessManager* pmanager = particle->GetProcessManager() ;

    // both postStep and alongStep action are required: because
    // of the use of ghost volumes. If no ghost, the postStep is sufficient.

    pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);
  }
}



/////////////////////////////////////////////////////////////////////////////

void Em8PhysicsList::SetCuts()
{
  G4Timer theTimer ;
  theTimer.Start() ;
  if (verboseLevel >0)
  {
    G4cout << "Em8PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma

   SetCutValue(cutForGamma,"gamma");

   SetCutValue(cutForElectron,"e-");
   SetCutValue(cutForElectron,"e+");

  if (verboseLevel>1)     DumpCutValuesTable();

  theTimer.Stop();
  G4cout.precision(6);
  G4cout << G4endl ;
  G4cout << "total time(SetCuts)=" << theTimer.GetUserElapsed() << " s " <<G4endl;

}

///////////////////////////////////////////////////////////////////////////

void Em8PhysicsList::SetGammaCut(G4double val)
{
  cutForGamma = val;
}

///////////////////////////////////////////////////////////////////////////

void Em8PhysicsList::SetElectronCut(G4double val)
{
  cutForElectron = val;
}

////////////////////////////////////////////////////////////////////////////

void Em8PhysicsList::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
  G4cout << " MaxChargedStep=" << MaxChargedStep << G4endl;
  G4cout << G4endl;
}

