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
// $Id: Tst23PhysicsList.cc,v 1.1 2001-12-14 14:53:42 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Tst23DetectorConstruction.hh"
#include "Tst23PhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   

#include "Tst23GeneralPhysics.hh"
#include "Tst23EMPhysics.hh"
#include "Tst23MuonPhysics.hh"
#include "Tst23HadronPhysics.hh"
#include "Tst23IonPhysics.hh"
#include "Tst23PhysicsListMessenger.hh"

Tst23PhysicsList::Tst23PhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  cutInRange[0] = defaultCutValue;
  cutInRange[1] = defaultCutValue;
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new Tst23GeneralPhysics("general") );

  // EM Physics
  RegisterPhysics( new Tst23EMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(  new Tst23MuonPhysics("muon"));

   // Hadron Physics
  RegisterPhysics(  new Tst23HadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics( new Tst23IonPhysics("ion"));
 
  myMessenger = new Tst23PhysicsListMessenger(this);
  
}

Tst23PhysicsList::~Tst23PhysicsList()
{
  if (myMessenger!=0) delete myMessenger;
}

void Tst23PhysicsList::SetCutInRangeForRegion( G4double cut, G4int region)
{
  if ((region>=0) && (region<N_Regions) && (cut>0)) {
    cutInRange[region] = cut; 
    ResetCuts();
  }
}

void Tst23PhysicsList::SetParticleCuts( G4double* cuts, G4ParticleDefinition* particle)
{
  if (fRetrievePhysicsTable && !fIsRestoredCutValues) {
#ifdef G4VERBOSE  
    if (verboseLevel>2){
      G4cout << "Tst23PhysicsList::SetParticleCuts  ";
      G4cout << " Retrieve Cut Values for ";
      G4cout << particle->GetParticleName() << G4endl;
    }
#endif
    RetrieveCutValues(directoryPhysicsTable, fStoredInAscii);
    fIsRestoredCutValues = true;
  }

  if (particle->GetEnergyCuts() == 0) {
    particle->SetCuts(defaultCutValue);
    for (size_t idx=0; idx<N_Regions ; idx++){
      G4String matName = Tst23DetectorConstruction::GetMaterialName(idx); 
      G4Material* mat = G4Material::GetMaterial(matName);
      particle->SetRangeCut(cuts[idx],mat);
    } 
  }    
}

void Tst23PhysicsList::SetCutValue(G4double* cuts, const G4String& name)
{
  G4ParticleDefinition* particle = theParticleTable->FindParticle(name);
#ifdef G4VERBOSE    
  if (particle != 0){
    if (verboseLevel >1) {
      G4cout << "Tst23PhysicsList::SetCutValue       :";
      G4cout << "Set cuts for " << name << G4endl;
    }
  } else {
    if (verboseLevel >0) {
      G4cout << "Tst23PhysicsList::SetCutValue       :";
      G4cout << name << " is not found in ParticleTable" << G4endl;
    }
  }
#endif

  if (particle != 0){
    if (!particle->IsShortLived()) {
      //set cut value
      SetParticleCuts( cuts ,particle );
      // build physics table
      BuildPhysicsTable(particle);
    }
  } 
}

void Tst23PhysicsList::SetCutValueForOthers(G4double* cuts)
{
 // check cut value is positive
  if ( (cuts[0] <= 0.0) ||  (cuts[1] <= 0.0)) {
#ifdef G4VERBOSE    
    if (verboseLevel >0){
      G4cout << "Tst23PhysicsList::SetCutValueForOthers: negative cut values";
      G4cout << "  :" << cuts[0]/mm ;
      G4cout << " , " << cuts[1]/mm << "[mm]" << G4endl;
    }
#endif
    return;
  }

#ifdef G4VERBOSE    
  if (verboseLevel >1) {
      G4cout << "Tst23PhysicsList::SetCutValueForOthers ";
      G4cout << "  :" << cuts[0]/mm ;
      G4cout << " , " << cuts[1]/mm << "[mm]" << G4endl;
  }
#endif

  // Sets a cut value to particle types which have not be called SetCuts() 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();

    if (!particle->IsShortLived()) {
      // check if the cut value has already been set
      if ((particle->GetLengthCuts()==0) ||(particle->GetEnergyCuts()==0)){
        // set cut value
        SetParticleCuts( cuts ,particle );
        // build physics table 
        BuildPhysicsTable(particle);

#ifdef G4VERBOSE    
        if (verboseLevel >1) 
          G4cout << "Set cuts for " << particle->GetParticleName() << G4endl;
#endif
      }
      // special treatment for particles for store cuts
      G4int idx = isParticleForStoreCuts(particle);
      if (idx >=0) {
        if (!isBuildPhysicsTable[idx])BuildPhysicsTable(particle);
      }
    }
  }
}



void Tst23PhysicsList::SetCuts()
{
#ifdef G4VERBOSE    
  if (verboseLevel >1){
    G4cout << "G4VUserPhysicsList::SetCutsWithDefault:";
    G4cout << "CutLength [0]: " << cutInRange[0]/mm << " (mm)" ;
    G4cout << "CutLength [1]: " << cutInRange[1]/mm << " (mm)" << G4endl;
  }  
#endif

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutInRange, "gamma");
  SetCutValue(cutInRange, "e-");
  SetCutValue(cutInRange, "e+");
 
  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton 
  SetCutValue(cutInRange, "proton");
  SetCutValue(cutInRange, "anti_proton");

  SetCutValueForOthers(cutInRange);

  if (verboseLevel>1) {
    DumpCutValuesTable();
  }

}







