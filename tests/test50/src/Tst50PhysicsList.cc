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
// $Id: Tst50PhysicsList.cc,v 1.6 2003-01-17 17:14:14 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "Tst50PhysicsList.hh"
#include "Tst50PhysicsListMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include"G4VeLowEnergyLoss.hh"
#include "G4Material.hh"
#include "G4ios.hh"
#include "g4std/iomanip"                

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50PhysicsList::Tst50PhysicsList(G4bool LowEn,G4bool range,G4bool SP,G4bool RY):  G4VUserPhysicsList()
{ RangeOn=range;
 Stopping= SP;
  Low=LowEn;
  RadiationY=RY;
  rangeOn = "off";
  defaultCutValue = 1.*mm;
   SetVerboseLevel(1);
 physicsListMessenger = new Tst50PhysicsListMessenger(this);
 // cutForElectron = defaultCutValue ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50PhysicsList::~Tst50PhysicsList()
{
  delete  physicsListMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50PhysicsList::ConstructLeptons()
{
  // leptons
  //  e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  // mu+/-
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  // nu_e
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  // nu_mu
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50PhysicsList::ConstructMesons()
{
  //  mesons
  //    light mesons
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50PhysicsList::ConstructBaryons()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "G4UserSpecialCuts.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Low energy processes for gamma and e-
// gamma
#include "G4LowEnergyRayleigh.hh" 
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"  
#include "G4LowEnergyGammaConversion.hh" 
// e-
#include "G4LowEnergyIonisation.hh" 
#include "G4LowEnergyBremsstrahlung.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tst50PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
    
      if ( Low==true)
	{
 lowePhot = new  G4LowEnergyPhotoElectric("LowEnPhotoElec");
      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh);
      pmanager->AddDiscreteProcess(lowePhot);
      pmanager->AddDiscreteProcess(new G4LowEnergyCompton);
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion);
      G4cout<<"Low Energy processes for gamma ray"<<G4endl;}
      else{
         pmanager->AddDiscreteProcess(new G4GammaConversion());
         pmanager->AddDiscreteProcess(new G4ComptonScattering());      
         pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
           G4cout<<"Standard  processes for gamma ray"<<G4endl;
      }
     
    } else if (particleName == "e-") {
 

//------------- Standard for e- --------------------//
          
  if(Low==true)
	{    
	
     loweIon  = new G4LowEnergyIonisation("LowEnergyIoni");
     loweBrem = new G4LowEnergyBremsstrahlung("LowEnBrem");
   

 
     G4VProcess*  multipleScattering= new G4MultipleScattering();
						
       pmanager->AddProcess(multipleScattering, -1, 1,1);
     
      pmanager->AddProcess(loweIon,     -1, 2,2);
      pmanager->AddProcess(loweBrem,    -1,-1,3);      
 if(RangeOn==true || Stopping==true || RadiationY==true){
   pmanager->RemoveProcess(multipleScattering);
   G4VeLowEnergyLoss::SetEnlossFluc(false);
  
 }}
      else{
     
 	
	   G4VProcess* theeminusIonisation         = new G4eIonisation();
	  G4VProcess* theeminusBremsstrahlung     = new G4eBremsstrahlung();
	  G4VProcess* theeminusMultipleScattering=new G4MultipleScattering();
	   pmanager->AddProcess( theeminusIonisation, -1, 2,2);
	  pmanager->AddProcess(theeminusBremsstrahlung,-1,-1,3);
pmanager->AddProcess(theeminusMultipleScattering,-1,1,1);

 if(RangeOn==true || Stopping==true|| RadiationY==true){pmanager->RemoveProcess(theeminusMultipleScattering);}
      }
    } else if (particleName == "e+") {
    //positron
      G4VProcess* theeplusMultipleScattering = new G4MultipleScattering();
      G4VProcess* theeplusIonisation         = new G4eIonisation();
      G4VProcess* theeplusBremsstrahlung     = new G4eBremsstrahlung();
      G4VProcess* theeplusAnnihilation       = new G4eplusAnnihilation();
      //
      // add processes
      pmanager->AddProcess(theeplusMultipleScattering,-1,1,1);
      pmanager->AddProcess(theeplusIonisation, -1, 2,2);
      pmanager->AddProcess(theeplusBremsstrahlung, -1,-1,3);
      pmanager->AddProcess(theeplusAnnihilation,0,-1,4);
    }

  


  }}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"
void Tst50PhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50PhysicsList::SetCuts()
{ 
if (verboseLevel >0)
  {
    G4cout << "Tst50PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  } 
 G4VUserPhysicsList::SetCutsWithDefault();
// method sets 
   //  the default cut value for all particle types 
  
 /*
 G4VUserPhysicsList::SetCutValue(cutForElectron,"e-");
  G4cout<<"arrivo dopo gli e-"<<G4endl;

   G4VUserPhysicsList:: SetCutValueForOthers(defaultCutValue);
   G4cout<<"arrivo dopo tutte"<<G4endl;
 */
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void  Tst50PhysicsList::SetRangeConditions(G4String val)
{ 
 rangeOn = val;
 ConstructEM();
 /*
if (rangeOn=="off")
      { G4VeLowEnergyLoss::SetEnlossFluc(false); 
      G4cout<<"no energy loss fluct"<<G4endl;}

if (rangeOn=="on")
      { G4VeLowEnergyLoss::SetEnlossFluc(true); 
      G4cout<<"energy loss fluct on"<<G4endl;}
 */
}

void Tst50PhysicsList::SetElectronCut(G4double val)
{
  ResetCuts();
 defaultCutValue = val;

}









