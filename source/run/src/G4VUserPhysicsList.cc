// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VUserPhysicsList.cc,v 1.1 1999-01-07 16:14:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
//	History
//        first version                   09 Jan. 1998 by H.Kurashige 
//        modified                        24 Jan. 1998 by H.Kurashige 
//        modified                        06 June 1998  by H.Kurashige 
//        add G4ParticleWithCuts::SetEnergyRange
//                                        18 June 1998  by H.Kurashige 
//       modifeid for short lived particles 27  June 1998  by H.Kurashige
//       G4BestUnit on output             12 nov. 1998  mma  
// ------------------------------------------------------------

#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleWithCuts.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BarionConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Material.hh"
#include "G4UserPhysicsListMessenger.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include <iomanip.h>                


G4VUserPhysicsList::G4VUserPhysicsList()
                   :verboseLevel(1)
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;

  // set energy range for SetCut calcuration
  G4ParticleWithCuts::SetEnergyRange(0.99*keV, 100*TeV);

  // pointer to the particle table
  theParticleTable = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();
 
  // UI Messenger 
  theMessenger = new G4UserPhysicsListMessenger(this);
}

G4VUserPhysicsList::~G4VUserPhysicsList()
{
  if (theMessenger != NULL) {
    delete theMessenger;
    theMessenger = NULL;
  }
}

void G4VUserPhysicsList::AddProcessManager(G4ParticleDefinition* newParticle,
					   G4ProcessManager*     newManager)
{
  if (newParticle == NULL) return;
  if (newManager  == NULL){
    newManager = new G4ProcessManager(newParticle);
    if (newParticle->GetParticleType() == "nucleus") {
      G4ParticleDefinition* genericIon = 
	   (G4ParticleTable::GetParticleTable())->FindParticle("GenericIon");
      if (genericIon != NULL) {
	G4ProcessManager* ionMan = genericIon->GetProcessManager();
	if (ionMan != NULL) {
	  delete newManager;
	  newManager = new G4ProcessManager(*ionMan);
	}
      }
    }
  }
  
  if (verboseLevel >2){
    G4cerr << "G4VUserPhysicsList::AddProcessManager: ";
    G4cerr  << "adds ProcessManager to ";
    G4cerr  << newParticle->GetParticleName() << endl;
    newManager->DumpInfo();
  } 
  newManager->SetParticleType(newParticle);
  if (newParticle->GetProcessManager() == NULL) {
    newParticle->SetProcessManager(newManager);
  } else {
    if (verboseLevel >0){
      G4cerr << "G4VUserPhysicsList::AddProcessManager: ";
      G4cerr  << newParticle->GetParticleName();
      G4cerr << " already has ProcessManager " << endl;
    }
  }
}


void G4VUserPhysicsList::InitializeProcessManager()
{
  // loop over all particles in G4ParticleTable 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if  (pmanager==NULL) {
      pmanager = new G4ProcessManager(particle);
      particle->SetProcessManager(pmanager);
    }
  }
}

#include "G4Transportation.hh"
void G4VUserPhysicsList::AddTransportation()
{
  G4Transportation* theTransportationProcess= new G4Transportation();

  // loop over all particles in G4ParticleTable 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (!particle->IsShortLived()) {
      pmanager ->AddProcess(theTransportationProcess);
      pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
      pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
    }
  }
}


void G4VUserPhysicsList::SetDefaultCutValue(G4double value)
{
   if (value<=0.0) {
     if (verboseLevel >0){
       G4cerr << "G4VUserPhysicsList::SetDefaultCutValue: negative cut values";
       G4cerr << "  :" << value/mm << "[mm]" << endl;
     }
   } else { 
     if (verboseLevel >1){
       G4cout << "G4VUserPhysicsList::SetDefaultCutValue:";
       G4cout << "default cut value is changed to   :" ;
       G4cout << value/mm << "[mm]" << endl;
     }
     defaultCutValue = value;
     ResetCuts();
   }
}

void G4VUserPhysicsList::ResetCuts()
{
  if (verboseLevel >1) {
    G4cout << "G4VUserPhysicsList::ResetCuts()" << endl;
    G4cout << "  cut values in energy will be calculated later" << endl;
  }

  // Reset cut values for other particles
  // loop over all particles in G4ParticleTable 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    if (!particle->IsShortLived()) particle->ResetCuts();
  }
    
  // inform that cut values are modified to the run manager
  // i.e state will be changed and SetCuts will be invoked 
  // just before event loop
  G4UImanager::GetUIpointer()->ApplyCommand("/run/cutoffModified");
}

void G4VUserPhysicsList::SetCutValueForOtherThan(G4double cutValue,
                                 G4ParticleDefinition* first,
                                 G4ParticleDefinition* second,
				 G4ParticleDefinition* third,
                                 G4ParticleDefinition* fourth,
				 G4ParticleDefinition* fifth,
                                 G4ParticleDefinition* sixth,
				 G4ParticleDefinition* seventh,
                                 G4ParticleDefinition* eighth,
				 G4ParticleDefinition* nineth,
                                 G4ParticleDefinition* tenth  )
{
  // check cut value is positive
  if (cutValue <= 0.0) {
    if (verboseLevel >0){
      G4cerr << "G4VUserPhysicsList::SetCutValueForOtherThan: negative cut values";
      G4cerr << "  :" << cutValue/mm << "[mm]" << endl;
    }
    return;
  } else {
    if (verboseLevel >1) {
      G4cout << "G4VUserPhysicsList::SetCutValueForOtherThan ";
      G4cout << "  :" << cutValue/mm << "[mm]" << endl;
    }
  }

  // check specified particle types in arguments
  G4ParticleDefinition* specifiedParticles[10];
  G4int numberOfSpecifiedParticles = 0;
  if (first != 0) {
    specifiedParticles[numberOfSpecifiedParticles] = first;
    numberOfSpecifiedParticles +=1;
  }
  if (second != 0) {
    specifiedParticles[numberOfSpecifiedParticles] = second;
    numberOfSpecifiedParticles +=1;
  }
  if (third != 0) {
    specifiedParticles[numberOfSpecifiedParticles] = third;
    numberOfSpecifiedParticles +=1;
  }
  if (fourth != 0) {
    specifiedParticles[numberOfSpecifiedParticles] = fourth;
    numberOfSpecifiedParticles +=1;
  }
  if (fifth != 0) {
    specifiedParticles[numberOfSpecifiedParticles] = fifth;
    numberOfSpecifiedParticles +=1;
  }
  if (sixth != 0) {
    specifiedParticles[numberOfSpecifiedParticles] = sixth;
    numberOfSpecifiedParticles +=1;
  }
  if (seventh != 0) {
    specifiedParticles[numberOfSpecifiedParticles] = seventh;
    numberOfSpecifiedParticles +=1;
  }
  if (eighth != 0) {
    specifiedParticles[numberOfSpecifiedParticles] = eighth;
    numberOfSpecifiedParticles +=1;
  }
  if (nineth != 0) {
    specifiedParticles[numberOfSpecifiedParticles] = nineth;
    numberOfSpecifiedParticles +=1;
  }
  if (tenth != 0) {
    specifiedParticles[numberOfSpecifiedParticles] = tenth;
    numberOfSpecifiedParticles +=1;
  }

   // Set cut values for other particles
  G4bool isSpecified;
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    isSpecified = false;
    for (G4int index = 0; index <numberOfSpecifiedParticles; index++) {
      isSpecified = isSpecified || (particle == specifiedParticles[index]);
    }
    if (!isSpecified) {
       if (!particle->IsShortLived()) {
	 particle->SetCuts(cutValue);
	 if (verboseLevel >1) 
	   G4cout << "Set cuts for " << particle->GetParticleName() << endl;
	 BuildPhysicsTable(particle);
       }
    }
  }
}

void G4VUserPhysicsList::SetCutValue(G4double aCut, const G4String& name)
{
  G4ParticleDefinition* particle;
  if (particle = theParticleTable->FindParticle(name)){
     if (!particle->IsShortLived()) {
       particle->SetCuts( aCut );
       if (verboseLevel >1) G4cout << "Set cuts for " << name << endl;
       BuildPhysicsTable(particle);
     }
  } else {
    if (verboseLevel >0) 
      G4cout << name << " is not found in ParticleTable" << endl;
  }
}


void G4VUserPhysicsList::SetCutValueForOthers(G4double cutValue)
{
 // check cut value is positive
  if (cutValue <= 0.0) {
    if (verboseLevel >0){
      G4cerr << "G4VUserPhysicsList::SetCutValueForOthers: negative cut values";
      G4cerr << "  :" << cutValue/mm << "[mm]" << endl;
    }
    return;
  } else {
    if (verboseLevel >1) {
      G4cout << "G4VUserPhysicsList::SetCutValueForOthers ";
      G4cout << "  :" << cutValue/mm << "[mm]" << endl;
    }
  }

  // Sets a cut value to particle types which have not be called SetCuts() 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    if (!particle->IsShortLived()) {
      if ((particle->GetLengthCuts()<0.0) ||(particle->GetEnergyCuts()==NULL)) {
	particle->SetCuts(cutValue);
	if (verboseLevel >1) 
	  G4cout << "Set cuts for " << particle->GetParticleName() << endl;
	BuildPhysicsTable(particle);
      }
    }
  }
}

void G4VUserPhysicsList::ReCalcCutValue(const G4String& name)
{
  G4ParticleDefinition* particle;
  if (particle = theParticleTable->FindParticle(name)){
    if (!particle->IsShortLived()) {
      particle->ReCalcCuts();
      if (verboseLevel >1) G4cout << "Recalc cuts for " << name << endl;
      BuildPhysicsTable(particle);
    }
  } else {
    if (verboseLevel >0) 
      G4cout << name << " is not found in ParticleTable" << endl;
  }
}

void G4VUserPhysicsList::ReCalcCutValueForOthers()
{
  if (verboseLevel >1) {
    G4cout << "G4VUserPhysicsList::ReCalcCutValueForOthers ";
  }

  // Sets a cut value to particle types which have not be called SetCuts() 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    if (!particle->IsShortLived()) {
      if (particle->GetEnergyCuts()==NULL) {
	if (verboseLevel >1) 
	  G4cout << "ReCalc cuts for " << particle->GetParticleName() << endl;
	particle->ReCalcCuts();
	BuildPhysicsTable(particle);
      }
    }
  }
}

void G4VUserPhysicsList::BuildPhysicsTable(G4ParticleDefinition* particle)
{
  G4int j;
  
  // Rebuild the physics tables for every process for this particle type
  G4ProcessVector* pVector = (particle->GetProcessManager())->GetProcessList();
  for ( j=0; j < pVector->length(); ++j) {
    (*pVector)[j]->BuildPhysicsTable(*particle);
  }
  for ( j=0; j < pVector->length(); ++j) {

    //*********************************************************************
    // temporary addition to make the integral schema of electromagnetic
    // processes work.
    //
    if((*pVector)[j]->GetProcessName() == "Imsc")
      {
	(*pVector)[j]->BuildPhysicsTable(*particle);
      }
    if((*pVector)[j]->GetProcessName() == "IeIoni")
      {
	(*pVector)[j]->BuildPhysicsTable(*particle);
      }
    if((*pVector)[j]->GetProcessName() == "IeBrems")
      {
	(*pVector)[j]->BuildPhysicsTable(*particle);
      }
    if((*pVector)[j]->GetProcessName() == "Iannihil")
      {
        (*pVector)[j]->BuildPhysicsTable(*particle);
      }
    if((*pVector)[j]->GetProcessName() == "IhIoni")
      {
        (*pVector)[j]->BuildPhysicsTable(*particle);
      }
    if((*pVector)[j]->GetProcessName() == "IMuIoni")
      {
        (*pVector)[j]->BuildPhysicsTable(*particle);
      }
    if((*pVector)[j]->GetProcessName() == "IMuBrems")
      {
        (*pVector)[j]->BuildPhysicsTable(*particle);
      }
    if((*pVector)[j]->GetProcessName() == "IMuPairProd")
      {
        (*pVector)[j]->BuildPhysicsTable(*particle);
      }
    //*********************************************************************

  }
}

void G4VUserPhysicsList::ConstructAllParticles()
{
  ConstructAllBosons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllBarions();
  ConstructAllIons();
  ConstructAllShortLiveds();
}

void G4VUserPhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void G4VUserPhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void G4VUserPhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void G4VUserPhysicsList::ConstructAllBarions()
{
  //  Construct all barions
  G4BarionConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void G4VUserPhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4VUserPhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();  
}


void G4VUserPhysicsList::DumpList() const
{
  theParticleIterator->reset();
  G4int idx = 0;
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4cout << particle->GetParticleName();
    if ((idx++ % 4) == 3) {
      G4cout << endl;
    } else {
      G4cout << ", ";
    }      
  }
  G4cout << endl;
}

void G4VUserPhysicsList::DumpCutValues(const G4String &particle_name) const
{
  G4ParticleDefinition* particle;
  if ((particle_name == "ALL") || (particle_name == "all")) {
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
      particle = theParticleIterator->value();
      DumpCutValues(particle);
    }
  } else {
     particle = theParticleTable->FindParticle(particle_name);
     if (particle != NULL) DumpCutValues(particle);
  }
}

void G4VUserPhysicsList::DumpCutValues( G4ParticleDefinition* particle) const
{
  if (particle == NULL) return;
  
  G4int prec = G4cout.precision(3);

  if (particle->IsShortLived()) {
    G4cout << " --- " << particle->GetParticleName() << " is a short lived particle ------ " << endl;
  } else {
    G4cout << " --- " << particle->GetParticleName() << " ------ " << endl;
    G4cout << "   - Cut in range = " << G4BestUnit(particle->GetLengthCuts(),"Length") << endl;
    G4double*  theKineticEnergyCuts = particle->GetEnergyCuts();
    
    if (theKineticEnergyCuts != NULL) {
      const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
      G4cout << "   - Material ---------------- Energy Cut ---" << endl;
      for (G4int idx=0; idx<materialTable->length(); idx++){
	G4cout << "     " << setw(19) << (*materialTable)[idx]->GetName(); 
	G4cout << " : "   << setw(10) << G4BestUnit(theKineticEnergyCuts[idx],"Energy");
	G4cout << endl;
      }
    } else {
      G4cout << "   - Cuts in energy are not calculated yet --" << endl;
      G4cout << " Enter /run/initialize command to calculate cuts " << endl;
    }
  }
  G4cout.precision(prec);
}

void G4VUserPhysicsList::DumpCutValuesTable() const
{
  // This methods Print out a table of cut value information
  // for "e-", "gamma", "mu-", "proton" and "neutron"
  G4int prec = G4cout.precision(3);
  const G4int size = 5;
  G4String name[size] = {"gamma", "e-", "mu-", "proton", "neutron"};
  G4ParticleDefinition* particle[size];
  G4bool IsOK = true;
  G4int idx; 
  for (idx=0; idx <size; idx++) {
    particle[idx] = theParticleTable->FindParticle(name[idx]);
  }

  //line 1 //-- Commented out (M.Asai)
  //G4cout << "Default cut value in range :";
  //G4cout << defaultCutValue/mm << "[mm]" << endl;

  //line 2
  G4cout << "============= The cut Energy ==============================" <<endl;
  
  // line 3
  G4cout << "                     ";
  for (idx=0; idx <size; idx++) {
    G4cout << " " << setw(11) << name[idx] << "    ";
  }
  G4cout << endl;

  // line 4
  G4cout << "Cut in range       ";
  for (idx=0; idx <size; idx++) {
    if (particle[idx] == NULL) {
      G4cout << "            ";
    } else {
      G4cout << " " << setw(11) << G4BestUnit(particle[idx]->GetLengthCuts(),"Length");
    }
  }
  G4cout << endl;

  // line 5
  G4cout << "Cut in energy";
  G4cout << endl;

 // line 6 ..
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  for (G4int J=0; J<materialTable->length(); J++) {
    G4cout << " " << setw(18) << ((*materialTable)[J])->GetName();
    for (idx=0; idx <size; idx++) {
      if (particle[idx] == NULL) {
	G4cout << "            ";
      } else {
        if (particle[idx]->GetEnergyCuts() == NULL) {
	  G4cout << " ---------- ";
          IsOK = false;
	} else {
	  G4cout << " " << setw(11) << G4BestUnit((particle[idx]->GetEnergyCuts())[J],"Energy");
	}
      }
    }
    G4cout << endl;
  }
  
  if (!IsOK) {
    G4cout << " Cuts in energy have not calculated yet !!" << endl;
    G4cout << " Enter /run/initialize command to calculate cuts " << endl;
  }

  // last line 
  G4cout << "===================================================" << endl;
  G4cout.precision(prec);
}













