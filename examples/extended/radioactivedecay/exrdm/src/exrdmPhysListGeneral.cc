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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "exrdmPhysListGeneral.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmPhysListGeneral::exrdmPhysListGeneral(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmPhysListGeneral::~exrdmPhysListGeneral()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysListGeneral::ConstructProcess()
{
  // Add Decay Process

  G4Decay* fDecayProcess = new G4Decay();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (fDecayProcess->IsApplicable(*particle)) { 

      pmanager ->AddProcess(fDecayProcess);

      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(fDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(fDecayProcess, idxAtRest);

    }
  }
  // Declare radioactive decay to the GenericIon in the IonTable.
  //
  const G4IonTable *theIonTable = G4ParticleTable::GetParticleTable()->GetIonTable();

  G4RadioactiveDecay*  theRadioactiveDecay = new G4RadioactiveDecay();
  for (G4int i=0; i<theIonTable->Entries(); i++)
    {
      G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
      if (particleName == "GenericIon")
	{
	  G4ProcessManager* pmanager = theIonTable->GetParticle(i)->GetProcessManager();
	  pmanager->SetVerboseLevel(0);
	  pmanager ->AddProcess(theRadioactiveDecay);
	  pmanager ->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
	  pmanager ->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
	}
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

