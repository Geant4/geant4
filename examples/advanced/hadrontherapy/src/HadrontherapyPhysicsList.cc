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
// $Id: HadrontherapyPhysicsList.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#include "globals.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4VeLowEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "G4ios.hh"
#include <iomanip.h>                

// ---------------------------------------------------------------------
HadrontherapyPhysicsList::HadrontherapyPhysicsList(HadrontherapyDetectorConstruction* p)
  :G4VUserPhysicsList()
{
  pDet = p;

  currentDefaultCut = defaultCutValue = 10 *mm;
  cutForGamma       = defaultCutValue;
  cutForElectron    = defaultCutValue;
  cutForProton      = defaultCutValue;
  SetVerboseLevel(1);
}

// ----------------------------------------------------------------------
HadrontherapyPhysicsList::~HadrontherapyPhysicsList()
{
}

// -----------------------------------------------------------------------
void HadrontherapyPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBarions();
  ConstructIons();
}

// -----------------------------------------------------------------------
void HadrontherapyPhysicsList::ConstructBosons()
{
  // ******* 
  //  gamma
  // *******
  G4Gamma::GammaDefinition();
  
  // **************** 
  //  optical photons
  // ****************
  G4OpticalPhoton::OpticalPhotonDefinition();
}

// ------------------------------------------------------------------------
void HadrontherapyPhysicsList::ConstructLeptons()
{
  // *******  
  // leptons
  // *******
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

// ------------------------------------------------------------------------
void HadrontherapyPhysicsList::ConstructMesons()
{
  // ********
  //  mesons
  // ********
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
}

// ------------------------------------------------------------------------
void HadrontherapyPhysicsList::ConstructBarions()
{
  // **********
  //  barions
  // **********
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}

// -----------------------------------------------------------------------
void HadrontherapyPhysicsList::ConstructIons()
{
  // ******
  //  ions
  // ******
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4He3::He3Definition();
  G4Alpha::AlphaDefinition();
  G4GenericIon::GenericIonDefinition();
}

// -----------------------------------------------------------------------
void HadrontherapyPhysicsList::ConstructProcess()
{
  AddTransportation(); 
  ConstructEM();
  ConstructHad();
  ConstructGeneral();
}

// -----------------------------------------------------------------------
// Electromagnetic processes valid for all charged particles
// -----------------------------------------------------------------------
#include "G4MultipleScattering.hh"
//  >>>  gamma  <<
#include "G4LowEnergyRayleigh.hh" 
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"  
#include "G4LowEnergyGammaConversion.hh" 

// >>> e- <<<
#include "G4LowEnergyIonisation.hh" 
#include "G4LowEnergyBremsstrahlung.hh" 

// >>> e+ <<<
#include "G4eplusAnnihilation.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"

// >>>  Alpha and generic Ions, deuteron, triton, He3  <<<
#include "G4hLowEnergyIonisation.hh"  
// hLowEnergyIonisation uses Ziegler 1988 as the default
#include "G4EnergyLossTables.hh"

// >>> muon  <<<<
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh"

// >>> Standard process for gamma   <<<
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

// >>>  Standard processes for hadro ionisation   <<<
#include "G4hIonisation.hh"

// ------------------------------------------------------------------------------
void HadrontherapyPhysicsList::ConstructEM()
{
  G4int LowEnergy = 1;
  // to active the Low Energy processes make LowEnergy =1;
  
  // ***********
  //  processes
  // ***********
  G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
  G4LowEnergyPhotoElectric* lowePhot        = new G4LowEnergyPhotoElectric();
  G4LowEnergyIonisation* loweIon            = new G4LowEnergyIonisation();
  G4LowEnergyBremsstrahlung* loweBrem       = new G4LowEnergyBremsstrahlung();
  G4hLowEnergyIonisation* ahadronLowEIon    = new G4hLowEnergyIonisation();

  if (LowEnergy == 1) 
   {    
    ahadronLowEIon -> SetNuclearStoppingPowerModel("ICRU_R49") ; // ICRU49 models for nuclear SP
    ahadronLowEIon -> SetNuclearStoppingOn() ;
  
    // setting tables explicitly for electronic stopping power
    ahadronLowEIon -> SetElectronicStoppingPowerModel(G4GenericIon::GenericIonDefinition(), 
						    "ICRU_R49p") ;  // ICRU49 models for elettronic SP
    ahadronLowEIon -> SetElectronicStoppingPowerModel(G4Proton::ProtonDefinition(), 
						    "ICRU_R49p") ;
    // Switch off the Barkas and Bloch corrections
    ahadronLowEIon -> SetBarkasOff();
  }  

  theParticleIterator -> reset();
  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      G4String particleType = particle->GetParticleType();
      G4double charge = particle->GetPDGCharge();

      if (particleName == "gamma") 
	{
	  // >>>  gamma  <<<
	 
	  if (LowEnergy == 1) {
	    // Low Energy processes
	    pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh());  
	    pmanager->AddDiscreteProcess(lowePhot);
	    pmanager->AddDiscreteProcess(new G4LowEnergyCompton());
	    pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion());
	  }
	 
	  else {
	    // Standard processes
	    pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
	    pmanager->AddDiscreteProcess(new G4ComptonScattering);
	    pmanager->AddDiscreteProcess(new G4GammaConversion);
	  }
	} 
     
      else if (particleName == "e-") 
       
	{
	  // >>>  electron  <<<
	 
	  if (LowEnergy == 1) {
	    // Low energy process
	    pmanager->AddProcess(aMultipleScattering,     -1, 1, 1);
	    pmanager->AddProcess(loweIon,                 -1, 2, 2);
	    pmanager->AddProcess(loweBrem,                -1,-1, 3);
	  }
	 
	  else {
	    // Standard processes
	    pmanager->AddProcess(aMultipleScattering,     -1, 1, 1);
	    pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
	    pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
	  }
	} 
     
      else if (particleName == "e+") 
	{
	  // >>>  positron  <<<
	 
	  // Standard processes
	  pmanager->AddProcess(aMultipleScattering,     -1, 1,1);
	  pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
	  pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
	  pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
	 
	} 
      else if( particleName == "mu+" || 
	       particleName == "mu-"    ) 
	{
	  // >>>  muon  <<<  
	  pmanager->AddProcess(aMultipleScattering,     -1, 1,1);
	  pmanager->AddProcess(new G4MuIonisation,      -1, 2,2);
	  pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,3);
	  pmanager->AddProcess(new G4MuPairProduction,  -1,-1,4);
	 
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
	 
	  if (LowEnergy == 1) { 
	    pmanager->AddProcess(aMultipleScattering,     -1,1,1); 
	    pmanager->AddProcess(ahadronLowEIon,       -1,2,2);
	  }
	  else {
	    pmanager->AddProcess(aMultipleScattering,     -1,1,1);
	    pmanager->AddProcess(new G4hIonisation,       -1,2,2);
	  }
	}
     
      else if ((!particle->IsShortLived()) &&
	       (charge != 0.0) ) 
	{
	  //all others charged particles except geantino
	 
	  if (LowEnergy == 1) {
	    pmanager->AddProcess(aMultipleScattering,-1,1,1);
	    pmanager->AddProcess(ahadronLowEIon,       -1,2,2);
	   
	  }
	  else {
	    pmanager->AddProcess(aMultipleScattering,-1,1,1); 
	    pmanager->AddProcess(new G4hIonisation(),-1, 2,2);
	  }
	 
	}
    }
}

// ----------------------------------------------------------------------
#include "G4Decay.hh"
void HadrontherapyPhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
   
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    
    //add decay process
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

// ----------------------------------------------------------------------
#include "G4ProtonInelasticProcess.hh"
#include "G4HadronElasticProcess.hh"
#include "G4LElastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4HEProtonInelastic.hh"

void HadrontherapyPhysicsList:: ConstructHad()
{
  G4int Hadronic = 1;
  // to activate the hadronic processes make Hadronic =1

  if(Hadronic == 1)
    {
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      G4LElastic* theElasticModel = new G4LElastic;
      theElasticProcess->RegisterMe(theElasticModel);
  
      theParticleIterator->reset();
      while ((*theParticleIterator)()) 
	{
	  G4ParticleDefinition* particle = theParticleIterator->value();
	  G4ProcessManager* pmanager = particle->GetProcessManager();
	  G4String particleName = particle->GetParticleName();
      
	  
	  if (particleName == "proton") 
	    {
	      pmanager->AddDiscreteProcess(theElasticProcess);
	      G4ProtonInelasticProcess* theInelasticProcess = new G4ProtonInelasticProcess("inelastic");
	      G4LEProtonInelastic* theLEInelasticModel = new G4LEProtonInelastic;
	      theInelasticProcess->RegisterMe(theLEInelasticModel);
	      G4HEProtonInelastic* theHEInelasticModel = new G4HEProtonInelastic;
	      theInelasticProcess->RegisterMe(theHEInelasticModel);
	      pmanager->AddDiscreteProcess(theInelasticProcess);
	    }
	}
    }
}
#include "G4Region.hh"
#include "G4RegionStore.hh"
void HadrontherapyPhysicsList::SetCuts()
{
  // reactualise cutValues
  if (currentDefaultCut != defaultCutValue)
    {
      if(cutForGamma    == currentDefaultCut) cutForGamma    = defaultCutValue;
      if(cutForElectron == currentDefaultCut) cutForElectron = defaultCutValue;
      if(cutForProton   == currentDefaultCut) cutForProton   = defaultCutValue;
      currentDefaultCut = defaultCutValue;
    }
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma,"gamma");
  SetCutValue(cutForElectron,"e-");
  SetCutValue(cutForElectron,"e+");

  
  // Cut per region
  // in region dosemeter we need a very accurate precision
  G4Region* region;
  G4String regName;
  G4ProductionCuts* cuts;
    
  regName = "Dosemeter";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(0.02*mm,G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(0.02*mm,G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(0.02*mm,G4ProductionCuts::GetIndex("e+"));
  region->SetProductionCuts(cuts);
      
  //  SetCutValueForOthers(defaultCutValue);        
   
  if (verboseLevel >0){
    G4cout << "HadrontherapyPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
    
  if (verboseLevel>0) DumpCutValuesTable();
}

// ---------------------------------------------------------------------------
void HadrontherapyPhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}

// ---------------------------------------------------------------------------
void HadrontherapyPhysicsList::SetElectronCut(G4double val)
{
  ResetCuts();
  cutForElectron = val;
}

// ---------------------------------------------------------------------------
void HadrontherapyPhysicsList::SetProtonCut(G4double val)
{
  ResetCuts();
  cutForProton = val;
}

// ---------------------------------------------------------------------------
void HadrontherapyPhysicsList::GetRange(G4double val)
{
  G4ParticleTable* theParticleTable =  G4ParticleTable::GetParticleTable();
  G4Material* currMat = pDet->GetDosemeterMaterial();

  G4ParticleDefinition* part;
  G4double cut;
  part = theParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  G4cout << "material : " << currMat->GetName() << G4endl;
  G4cout << "particle : " << part->GetParticleName() << G4endl;
  G4cout << "energy   : " << G4BestUnit(val,"Energy") << G4endl;
  G4cout << "range    : " << G4BestUnit(cut,"Length") << G4endl;
}


