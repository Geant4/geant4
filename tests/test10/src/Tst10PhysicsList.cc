// $Id: Tst10PhysicsList.cc,v 1.1 1999-01-08 16:35:33 gunter Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This is a version for maximum particle set
//	History
//        first version              10  Jan. 1998 by H.Kurashige
//        add decay at rest          26  Feb. 1998 by H.Kurashige
// ------------------------------------------------------------

#include "globals.hh"
#include "Tst10PhysicsList.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip.h>                


Tst10PhysicsList::Tst10PhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

Tst10PhysicsList::~Tst10PhysicsList()
{
}

void Tst10PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();

}

void Tst10PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}



void Tst10PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}
#include "G4OpBoundaryProcess.hh"

void Tst10PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "opticalphoton") {
    //opticalphoton
      // Construct processes for opticalphoton
			pmanager->AddDiscreteProcess(new G4OpBoundaryProcess());
 
		} 
  }
}

#include "G4Decay.hh"
void Tst10PhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager->AddProcess(theDecayProcess, INT_MAX, -1, INT_MAX); 
   }
  }
}

void Tst10PhysicsList::SetCuts(G4double cut)
{
  if (verboseLevel >0){
    G4cout << "Tst10PhysicsList::SetCuts:";
    G4cout << "CutLength : " << cut/mm << " (mm)" << endl;
  }  
  // uppress error messages even in case e/gamma/proton do not exist            
  G4int temp = GetVerboseLevel();                                               
  SetVerboseLevel(0);                                                           
                                                                                
  // set cut values for gamma at first and for e- second and next for e+,       
  // because some processes for e+/e- need cut values for gamma                 
  SetCutValue(cut, "gamma");                                                    
  SetCutValue(cut, "e-");                                                       
  SetCutValue(cut, "e+");                                                       
                                                                                
  // set cut values for proton and anti_proton before all other hadrons         
  // because some processes for hadrons need cut values for proton/anti_proton  
  SetCutValue(cut, "proton");                                                   
  SetCutValue(cut, "anti_proton");                                              
                                                                                
  // set a cut value to all particles                                           
                                                                                
  SetCutValueForOthers(cut);                                                    
                                                                                
  SetVerboseLevel(temp);                                                        

  if (verboseLevel>1) {
    DumpCutValuesTable();
  }
}


