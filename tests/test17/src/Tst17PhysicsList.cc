// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17PhysicsList.cc,v 1.1 1999-11-30 18:01:55 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
   
#include "G4Timer.hh"

#include "Tst17PhysicsList.hh"
#include "Tst17DetectorConstruction.hh"
#include "Tst17PhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4EnergyLossTables.hh"
#include "G4ios.hh"
#include <iomanip.h>                

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst17PhysicsList::Tst17PhysicsList(Tst17DetectorConstruction* p):  
 G4VUserPhysicsList(),
 theLEPhotoElectric(0), theLECompton(0), theLEGammaConversion(0), 
 theLERayleigh(0), theLEIonisation(0), theeminusMultipleScattering(0),theLEBremsstrahlung(0),
 theeplusMultipleScattering(0),theeplusIonisation(0),
 theeplusBremsstrahlung(0),
 theeplusAnnihilation(0),
 MaxChargedStep(DBL_MAX)
{
  pDet = p;

  defaultCutValue = 1.000*mm;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  cutForProton = defaultCutValue;

  SetVerboseLevel(1);
  physicsListMessenger = new Tst17PhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst17PhysicsList::~Tst17PhysicsList()
{
  delete physicsListMessenger; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17PhysicsList::ConstructBosons()
{
  // gamma
  G4Gamma::GammaDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"

#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17PhysicsList::ConstructEM()
{
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma

      theLEPhotoElectric   = new G4LowEnergyPhotoElectric();      
      theLECompton         = new G4LowEnergyCompton();
      theLEGammaConversion = new G4LowEnergyGammaConversion();
      theLERayleigh        = new G4LowEnergyRayleigh();

      pmanager->AddDiscreteProcess(theLEPhotoElectric);
      pmanager->AddDiscreteProcess(theLECompton);
      pmanager->AddDiscreteProcess(theLERayleigh);
      pmanager->AddDiscreteProcess(theLEGammaConversion);
      
    } 
    else if (particleName == "e-") {

      theLEIonisation = new G4LowEnergyIonisation();
      theLEBremsstrahlung = new G4LowEnergyBremsstrahlung();
      theeminusMultipleScattering = new G4MultipleScattering();

      pmanager->AddProcess(theeminusMultipleScattering,-1,1,1);
      pmanager->AddProcess(theLEIonisation,-1,2,2);
      pmanager->AddProcess(theLEBremsstrahlung,-1,-1,3);      

    }
    else if (particleName == "e+") {
    // Construct processes for positron <------------------------------
      theeplusMultipleScattering = new G4MultipleScattering();
      theeplusIonisation = new G4eIonisation();
      theeplusBremsstrahlung = new G4eBremsstrahlung();
      theeplusAnnihilation = new G4eplusAnnihilation();
      
      pmanager->AddProcess(theeplusMultipleScattering,-1,1,1);
      pmanager->AddProcess(theeplusIonisation,-1,2,2);
      pmanager->AddProcess(theeplusBremsstrahlung,-1,-1,3);
      pmanager->AddProcess(theeplusAnnihilation,0,-1,4);      
    } 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17PhysicsList::SetGELowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "Tst17PhysicsList::SetCuts:";
    G4cout << "Gamma and Electron cut in energy: " << lowcut*MeV << " (MeV)" << endl;
  }  

  G4Gamma::SetEnergyRange(lowcut*MeV,1e5*MeV);
  G4Electron::SetEnergyRange(lowcut*MeV,1e5*MeV);

}

void Tst17PhysicsList::SetGammaLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "Tst17PhysicsList::SetCuts:";
    G4cout << "Gamma cut in energy: " << lowcut*MeV << " (MeV)" << endl;
  }  

  G4Gamma::SetEnergyRange(lowcut*MeV,1e5*MeV);

}

void Tst17PhysicsList::SetElectronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){

    G4cout << "Tst17PhysicsList::SetCuts:";
    G4cout << "Electron cut in energy: " << lowcut*MeV << " (MeV)" << endl;

  }  

  G4Electron::SetEnergyRange(lowcut*MeV,1e5*MeV);

}

void Tst17PhysicsList::SetLowEnSecPhotCut(G4double cut){
  
  G4cout<<"Low energy secondary photons cut is now set to: "<<cut*MeV<<" (MeV)"<<endl;
  G4cout<<"for processes LowEnergyBremsstrahlung, LowEnergyPhotoElectric, LowEnergyIonisation"<<endl;
  theLEBremsstrahlung->SetCutForLowEnSecPhotons(cut);
  theLEPhotoElectric->SetCutForLowEnSecPhotons(cut);
  theLEIonisation->SetCutForLowEnSecPhotons(cut);
}

void Tst17PhysicsList::SetLowEnSecElecCut(G4double cut){
  
  G4cout<<"Low energy secondary electrons cut is now set to: "<<cut*MeV<<" (MeV)"<<endl;
  //  G4cout<<"for processes LowEnergyBremsstrahlung, LowEnergyPhotoElectric, LowEnergyIonisation"<<endl;
  G4cout<<"for processes LowEnergyIonisation"<<endl;
  //  G4LowEnergyPhotoElectric::SetCutForLowEnSecElectrons(cut);
  theLEIonisation->SetCutForLowEnSecElectrons(cut);
}

void Tst17PhysicsList::SetCuts()
{
  G4Timer theTimer ;
  theTimer.Start() ;
  if (verboseLevel >0){
    G4cout << "Tst17PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << endl;
  }  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma,"gamma");
  
  SetCutValue(cutForElectron,"e-");
  SetCutValue(cutForElectron,"e+");
  
  if (verboseLevel>1) {
    DumpCutValuesTable();
  }
  
  theTimer.Stop();
  G4cout.precision(6);
  G4cout << endl ;
  G4cout << "total time(SetCuts)=" << theTimer.GetUserElapsed() << " s " <<endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17PhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17PhysicsList::SetElectronCut(G4double val)
{
  ResetCuts();
  cutForElectron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17PhysicsList::SetProtonCut(G4double val)
{
  ResetCuts();
  cutForProton = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17PhysicsList::SetCutsByEnergy(G4double val)
{
  G4ParticleTable* theParticleTable =  G4ParticleTable::GetParticleTable();
  G4Material* currMat = pDet->GetAbsorberMaterial();

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  G4ParticleDefinition* part;
  G4double cut;

  part = theParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "e-");

  part = theParticleTable->FindParticle("e+");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "e+");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17PhysicsList::GetRange(G4double val)
{

  G4ParticleTable* theParticleTable =  G4ParticleTable::GetParticleTable();
  G4Material* currMat = pDet->GetAbsorberMaterial();

  G4ParticleDefinition* part;
  G4double cut;
  part = theParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  G4cout << "material : " << currMat->GetName() << endl;
  G4cout << "particle : " << part->GetParticleName() << endl;
  G4cout << "energy   : " << val / keV << " (keV)" << endl;
  G4cout << "range    : " << cut / mm << " (mm)" << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


