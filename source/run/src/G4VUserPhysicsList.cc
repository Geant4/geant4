// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VUserPhysicsList.cc,v 1.9 2000-10-19 13:30:23 kurasige Exp $
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
//       Added RemoveProcessManager        9 Feb. 1999 by H.Kurashige
//       Fixed RemoveProcessManager       15 Apr. 1999 by H.Kurashige
//       Removed ConstructAllParticles()  15 Apr. 1999 by H.Kurashige
// ------------------------------------------------------------

#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleWithCuts.hh"
#include "G4Material.hh"
#include "G4UserPhysicsListMessenger.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "g4std/iomanip"                


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
  if (theMessenger != 0) {
    delete theMessenger;
    theMessenger = 0;
  }
}

void G4VUserPhysicsList::AddProcessManager(G4ParticleDefinition* newParticle,
					   G4ProcessManager*     newManager)
{
  if (newParticle == 0) return;
  if (newParticle->GetProcessManager() != 0) {
#ifdef G4VERBOSE    
    if (verboseLevel >1){
      G4cout << "G4VUserPhysicsList::AddProcessManager: ";
      G4cout  << newParticle->GetParticleName();
      G4cout << " already has ProcessManager " << G4endl;
    }
#endif
    return;
  }

  // create new process manager if newManager  == 0
  if (newManager  == 0){
    // Add ProcessManager
    if (newParticle->GetParticleType() == "nucleus") {
      // Create a copy of the process manager of "GenericIon" in case of "nucleus"
      G4ParticleDefinition* genericIon = 
	   (G4ParticleTable::GetParticleTable())->FindParticle("GenericIon");

      if (genericIon != 0) {
	G4ProcessManager* ionMan = genericIon->GetProcessManager();
	if (ionMan != 0) {
	  newManager = new G4ProcessManager(*ionMan);
	} else {
	  // no process manager has been registered yet 
	  newManager = new G4ProcessManager(newParticle);
	}
      } else {
	// "GenericIon" does not exist
	newManager = new G4ProcessManager(newParticle);
      }

    } else {
      // create process manager for particles other than "nucleus"
      newManager = new G4ProcessManager(newParticle);
    }
  }

  // set particle type   
  newManager->SetParticleType(newParticle);

  // add the process manager
  newParticle->SetProcessManager(newManager);

#ifdef G4VERBOSE    
  if (verboseLevel >2){
    G4cout << "G4VUserPhysicsList::AddProcessManager: ";
    G4cout  << "adds ProcessManager to ";
    G4cout  << newParticle->GetParticleName() << G4endl;
    newManager->DumpInfo();
  } 
#endif

}


void G4VUserPhysicsList::InitializeProcessManager()
{
  // loop over all particles in G4ParticleTable 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if  (pmanager==0) {
      // create process manager if the particle has no its one
      pmanager = new G4ProcessManager(particle);
      particle->SetProcessManager(pmanager);
    }
  }
}

void G4VUserPhysicsList::RemoveProcessManager()
{
  // loop over all particles in G4ParticleTable 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if  (pmanager!=0) delete pmanager;
    particle->SetProcessManager(0);
#ifdef G4VERBOSE    
    if (verboseLevel >2){
      G4cout << "G4VUserPhysicsList::RemoveProcessManager: ";
      G4cout  << "remove ProcessManager from ";
      G4cout  << particle->GetParticleName() << G4endl;
    }
#endif
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
      // Add transportation process for all particles other than  "shortlived"
      if ( pmanager == 0) {
	// Error !! no process manager
	G4Exception("G4VUserPhysicsList::AddTransportation : no process manager!");
      } else {
	// add transportation with ordering = ( -1, "first", "first" )
	pmanager ->AddProcess(theTransportationProcess);
	pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
	pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
      }
    } else {
      // shortlived particle case
    }
  }
}


void G4VUserPhysicsList::SetDefaultCutValue(G4double value)
{
   if (value<=0.0) {
#ifdef G4VERBOSE    
     if (verboseLevel >0){
       G4cout << "G4VUserPhysicsList::SetDefaultCutValue: negative cut values";
       G4cout << "  :" << value/mm << "[mm]" << G4endl;
     }
#endif
   } else { 
#ifdef G4VERBOSE    
     if (verboseLevel >1){
       G4cout << "G4VUserPhysicsList::SetDefaultCutValue:";
       G4cout << "default cut value is changed to   :" ;
       G4cout << value/mm << "[mm]" << G4endl;
     }
#endif
     defaultCutValue = value;
     ResetCuts();
   }
}

void G4VUserPhysicsList::ResetCuts()
{
#ifdef G4VERBOSE    
  if (verboseLevel >1) {
    G4cout << "G4VUserPhysicsList::ResetCuts()" << G4endl;
    G4cout << "  cut values in energy will be calculated later" << G4endl;
  }
#endif

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
  
  //   set /control/verbose 0
  G4int tempVerboseLevel = G4UImanager::GetUIpointer()->GetVerboseLevel();
  G4UImanager::GetUIpointer()->SetVerboseLevel(0);
  //   issue /run/cutoffModified
  G4UImanager::GetUIpointer()->ApplyCommand("/run/cutoffModified");
  //   retreive  /control/verbose 
  G4UImanager::GetUIpointer()->SetVerboseLevel(tempVerboseLevel);
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
#ifdef G4VERBOSE    
    if (verboseLevel >0){
      G4cout << "G4VUserPhysicsList::SetCutValueForOtherThan: negative cut values";
      G4cout << "  :" << cutValue/mm << "[mm]" << G4endl;
    }
#endif
    return;
  } else {
#ifdef G4VERBOSE    
    if (verboseLevel >1) {
      G4cout << "G4VUserPhysicsList::SetCutValueForOtherThan ";
      G4cout << "  :" << cutValue/mm << "[mm]" << G4endl;
    }
#endif
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
    if ( (!isSpecified) && (!particle->IsShortLived()) ){
      // set cut value
      particle->SetCuts(cutValue);
      // build physics table
      BuildPhysicsTable(particle);

#ifdef G4VERBOSE    
      if (verboseLevel >1) G4cout << "Set cuts for " << particle->GetParticleName() << G4endl;
#endif
    }
  }
}

void G4VUserPhysicsList::SetCutValue(G4double aCut, const G4String& name)
{
  G4ParticleDefinition* particle = theParticleTable->FindParticle(name);
  if (particle != 0){
    if (!particle->IsShortLived()) {
      //set cut value
      particle->SetCuts( aCut );
      // build physics table
      BuildPhysicsTable(particle);
    }
  } 

#ifdef G4VERBOSE    
  if (particle != 0){
    if (verboseLevel >1) G4cout << "Set cuts for " << name << G4endl;
  } else {
    if (verboseLevel >0) 
      G4cout << name << " is not found in ParticleTable" << G4endl;
  }
#endif
}

void G4VUserPhysicsList::SetCutsWithDefault()
{
  // default cut value
   G4double cut = defaultCutValue;

#ifdef G4VERBOSE    
  if (verboseLevel >1){
    G4cout << "G4VUserPhysicsList::SetCutsWithDefault:";
    G4cout << "CutLength : " << cut/mm << " (mm)" << G4endl;
  }  
#endif

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cut, "gamma");
  SetCutValue(cut, "e-");
  SetCutValue(cut, "e+");
 
  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton 
  SetCutValue(cut, "proton");
  SetCutValue(cut, "anti_proton");
  
  SetCutValueForOthers(cut);

  if (verboseLevel>1) {
    DumpCutValuesTable();
  }
}  

void G4VUserPhysicsList::SetCutValueForOthers(G4double cutValue)
{
 // check cut value is positive
  if (cutValue <= 0.0) {
#ifdef G4VERBOSE    
    if (verboseLevel >0){
      G4cout << "G4VUserPhysicsList::SetCutValueForOthers: negative cut values";
      G4cout << "  :" << cutValue/mm << "[mm]" << G4endl;
    }
#endif
    return;
  }

#ifdef G4VERBOSE    
  if (verboseLevel >1) {
      G4cout << "G4VUserPhysicsList::SetCutValueForOthers ";
      G4cout << "  :" << cutValue/mm << "[mm]" << G4endl;
  }
#endif

  // Sets a cut value to particle types which have not be called SetCuts() 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();

    if (!particle->IsShortLived()) {
      // check if the cut value has already been set
      if ((particle->GetLengthCuts()<0.0) ||(particle->GetEnergyCuts()==0)) {
	// set cut value
	particle->SetCuts(cutValue);
	// build physics table 
	BuildPhysicsTable(particle);

#ifdef G4VERBOSE    
	if (verboseLevel >1) 
	  G4cout << "Set cuts for " << particle->GetParticleName() << G4endl;
#endif
      }
    }

  }
}

void G4VUserPhysicsList::ReCalcCutValue(const G4String& name)
{
  G4ParticleDefinition* particle = theParticleTable->FindParticle(name);
  if (particle != 0 ){
    if (!particle->IsShortLived()) {
      particle->ReCalcCuts();
      BuildPhysicsTable(particle);

#ifdef G4VERBOSE    
      if (verboseLevel >1) G4cout << "Recalc cuts for " << name << G4endl;
#endif
    }
  } 

#ifdef G4VERBOSE    
  if (( particle ==0 ) && (verboseLevel >0)) {
    G4cout << name << " is not found in ParticleTable" << G4endl;
  }
#endif

}

void G4VUserPhysicsList::ReCalcCutValueForOthers()
{
#ifdef G4VERBOSE    
  if (verboseLevel >1) {
    G4cout << "G4VUserPhysicsList::ReCalcCutValueForOthers ";
  }
#endif

  // Sets a cut value to particle types which have not be called SetCuts() 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();

    if (!particle->IsShortLived()) {
      if (particle->GetEnergyCuts()==0) {
	particle->ReCalcCuts();
	BuildPhysicsTable(particle);

#ifdef G4VERBOSE    
	if (verboseLevel >1) 
	  G4cout << "ReCalc cuts for " << particle->GetParticleName() << G4endl;
#endif
      }
    }

  }
}

void G4VUserPhysicsList::BuildPhysicsTable(G4ParticleDefinition* particle)
{
  G4int j;
  
  // Rebuild the physics tables for every process for this particle type
  G4ProcessVector* pVector = (particle->GetProcessManager())->GetProcessList();
  for ( j=0; j < pVector->entries(); ++j) {
    (*pVector)[j]->BuildPhysicsTable(*particle);
  }

  for ( j=0; j < pVector->entries(); ++j) {

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
void G4VUserPhysicsList::DumpList() const
{
  theParticleIterator->reset();
  G4int idx = 0;
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4cout << particle->GetParticleName();
    if ((idx++ % 4) == 3) {
      G4cout << G4endl;
    } else {
      G4cout << ", ";
    }      
  }
  G4cout << G4endl;
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
     if (particle != 0) DumpCutValues(particle);
  }
}

void G4VUserPhysicsList::DumpCutValues( G4ParticleDefinition* particle) const
{
  if (particle == 0) return;
  
  G4int prec = G4cout.precision(3);

  if (particle->IsShortLived()) {
    // name field
    G4cout << " --- " << particle->GetParticleName() << " is a short lived particle ------ " << G4endl;
  } else {
    // name field
    G4cout << " --- " << particle->GetParticleName() << " ------ " << G4endl;

    // cut value in range field
    G4cout << "   - Cut in range = " << G4BestUnit(particle->GetLengthCuts(),"Length") << G4endl;

    // material and energy cut value for the material 
    G4double*  theKineticEnergyCuts = particle->GetEnergyCuts();
    
    if (theKineticEnergyCuts != 0) {
      const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
      G4cout << "   - Material ---------------- Energy Cut ---" << G4endl;
      for (G4int idx=0; idx<materialTable->entries(); idx++){
	G4cout << "     " << G4std::setw(19) << (*materialTable)[idx]->GetName(); 
	G4cout << " : "   << G4std::setw(10) << G4BestUnit(theKineticEnergyCuts[idx],"Energy");
	G4cout << G4endl;
      }

    } else {
      G4cout << "   - Cuts in energy are not calculated yet --" << G4endl;
      G4cout << " Enter /run/initialize command to calculate cuts " << G4endl;
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
  G4int size_display=2;
  G4int idx; 
  for (idx=0; idx <size_display; idx++) {
    particle[idx] = theParticleTable->FindParticle(name[idx]);
  }

  //line 1 //-- Commented out (M.Asai)
  //G4cout << "Default cut value in range :";
  //G4cout << defaultCutValue/mm << "[mm]" << G4endl;

  //line 2
  G4cout << "============= The cut Energy ==============================" <<G4endl;
  
  // line 3
  G4cout << "                     ";
  for (idx=0; idx <size_display; idx++) {
    G4cout << " " << G4std::setw(11) << name[idx] << "    ";
  }
  G4cout << G4endl;

  // line 4
  G4cout << "Cut in range       ";
  for (idx=0; idx <size_display; idx++) {
    if (particle[idx] == 0) {
      G4cout << "            ";
    } else {
      G4cout << " " << G4std::setw(11) << G4BestUnit(particle[idx]->GetLengthCuts(),"Length");
    }
  }
  G4cout << G4endl;

  // line 5
  G4cout << "Cut in energy";
  G4cout << G4endl;

 // line 6 ..
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  for (G4int J=0; J<materialTable->entries(); J++) {
    G4cout << " " << G4std::setw(18) << ((*materialTable)[J])->GetName();
    for (idx=0; idx <size_display; idx++) {
      if (particle[idx] == 0) {
	G4cout << "            ";
      } else {
        if (particle[idx]->GetEnergyCuts() == 0) {
	  G4cout << " ---------- ";
          IsOK = false;
	} else {
	  G4cout << " " << G4std::setw(11) << G4BestUnit((particle[idx]->GetEnergyCuts())[J],"Energy");
	}
      }
    }
    G4cout << G4endl;
  }
  
  if (!IsOK) {
    G4cout << " Cuts in energy have not calculated yet !!" << G4endl;
    G4cout << " Enter /run/initialize command to calculate cuts " << G4endl;
  }

  // last line 
  G4cout << "===================================================" << G4endl;
  G4cout.precision(prec);
}













