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
// $Id: G4VUserPhysicsList.cc,v 1.28 2003-01-14 22:46:34 asaim Exp $
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
//       modified                         08, Mar 2001 by H.Kurashige
//       modified for CUTS per REGION     10, Oct 2002 by H.Kurashige
// ------------------------------------------------------------

#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"
#include "G4UserPhysicsListMessenger.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProductionCuts.hh"
#include "G4MaterialCutsCouple.hh"

#include "G4ios.hh"
#include "g4std/iomanip"                
#include "g4std/fstream"                

////////////////////////////////////////////////////////
G4VUserPhysicsList::G4VUserPhysicsList()
                   :verboseLevel(1), 
		    fRetrievePhysicsTable(false),
		    fStoredInAscii(true),
		    fIsCheckedForRetrievePhysicsTable(false),
		    fIsRestoredCutValues(false),
                    directoryPhysicsTable("."),
                    fDisplayThreshold(0)
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;


  // pointer to the particle table
  theParticleTable = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();
  
  // pointer to the cuts table
  fCutsTable =  G4ProductionCutsTable::GetProductionCutsTable();

  // set energy range for SetCut calcuration
  fCutsTable->SetEnergyRange(0.99*keV, 100*TeV);

  // UI Messenger 
  theMessenger = new G4UserPhysicsListMessenger(this);
}

////////////////////////////////////////////////////////
G4VUserPhysicsList::~G4VUserPhysicsList()
{
  if (theMessenger != 0) {
    delete theMessenger;
    theMessenger = 0;
  }
}

////////////////////////////////////////////////////////
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


////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////
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


////////////////////////////////////////////////////////
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


////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////
void G4VUserPhysicsList::ResetCuts()
{
  //   set /control/verbose 0
  G4int tempVerboseLevel = G4UImanager::GetUIpointer()->GetVerboseLevel();
  G4UImanager::GetUIpointer()->SetVerboseLevel(0);
  //   issue /run/cutoffModified
  G4UImanager::GetUIpointer()->ApplyCommand("/run/cutoffModified");
  //   retreive  /control/verbose 
  G4UImanager::GetUIpointer()->SetVerboseLevel(tempVerboseLevel);
}




////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetCutValue(G4double aCut, const G4String& name)
{
  G4ParticleDefinition* particle = theParticleTable->FindParticle(name);
  if (particle != 0){
    if (!particle->IsShortLived()) {
      //set cut value
      SetParticleCuts( aCut ,particle );
    }
  } 
}



////////////////////////////////////////////////////////
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

  // set cut values for gamma at first and for e- and e+
  SetCutValue(cut, "gamma");
  SetCutValue(cut, "e-");
  SetCutValue(cut, "e+");
 
  // set cut values for proton and anti_proton
  SetCutValue(cut, "proton");
  SetCutValue(cut, "anti_proton");
  SetCutValue(cut, "neutron");
  SetCutValue(cut, "anti_neutron");
  

  //if (verboseLevel>1) {
  //  DumpCutValuesTable();
  //}
}  


////////////////////////////////////////////////////////
void  G4VUserPhysicsList::SetApplyCuts(G4bool value, const G4String& name)
{
  // Sets a value to particle 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();

    // check if the cut value has already been set
    if (name =="all" || name =="ALL" ){
      particle->SetApplyCutsFlag(value);
    } else if (name == particle->GetParticleName()){
      particle->SetApplyCutsFlag(value);
#ifdef G4VERBOSE    
      if (verboseLevel >1) {
	G4cout << "Set ApplyCutsFLag ";
	if (value) {
	  G4cout << "TRUE ";
	} else {
	  G4cout << "FALSE ";
	}
	G4cout << "for " << particle->GetParticleName() << G4endl;
      }
#endif
      break;
    }
  }
}
////////////////////////////////////////////////////////
G4bool G4VUserPhysicsList::GetApplyCuts(const G4String& name) const
{
  // Sets a value to particle 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();

    // check if the cut value has already been set
    if (name == particle->GetParticleName()){
      return particle->GetApplyCutsFlag();
    }
  }
#ifdef G4VERBOSE    
  if (verboseLevel >0){
    G4cout << "G4VUserPhysicsList::GetApplyCuts: ";
    G4cout << " particle " <<  name  ;
    G4cout << " does not exist " << G4endl;
  }
#endif
  return false;
}



////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetParticleCuts( G4double cut, G4ParticleDefinition* particle)
{
  if (fRetrievePhysicsTable && !fIsRestoredCutValues) {
#ifdef G4VERBOSE  
    if (verboseLevel>2){
      G4cout << "G4VUserPhysicsList::SetParticleCuts  ";
      G4cout << " Retrieve Cut Values for ";
      G4cout << particle->GetParticleName() << G4endl;
    }
#endif
    fCutsTable->RetrieveCutsTable(directoryPhysicsTable, fStoredInAscii);
    fIsRestoredCutValues = true;
  } else {
    particle->SetCuts(cut);
  } 
}    

///////////////////////////////////////////////////////////////
void G4VUserPhysicsList::BuildPhysicsTable()
{
  // Sets a value to particle 
  // set cut values for gamma at first and for e- and e+
  G4String particleName;
  G4ParticleDefinition* GammaP = theParticleTable->FindParticle(particleName="gamma");
  if(GammaP) BuildPhysicsTable(GammaP);
  G4ParticleDefinition* EMinusP = theParticleTable->FindParticle(particleName="e-");
  if(EMinusP) BuildPhysicsTable(EMinusP);
  G4ParticleDefinition* EPlusP = theParticleTable->FindParticle(particleName="e+");
  if(EPlusP) BuildPhysicsTable(EPlusP);
  G4ParticleDefinition* ProtonP = theParticleTable->FindParticle(particleName="proton");
  if(ProtonP) BuildPhysicsTable(ProtonP);
  G4ParticleDefinition* AntiProtonP = theParticleTable->FindParticle(particleName="anti_proton");
  if(AntiProtonP) BuildPhysicsTable(AntiProtonP);
  G4ParticleDefinition* NeutronP = theParticleTable->FindParticle(particleName="newtron");
  if(NeutronP) BuildPhysicsTable(NeutronP);
  G4ParticleDefinition* AntiNeutronP = theParticleTable->FindParticle(particleName="anti_neutron");
  if(AntiNeutronP) BuildPhysicsTable(AntiNeutronP);
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    if(particle!=GammaP && particle!=EMinusP && particle!=EPlusP && particle!=ProtonP
    && particle!=AntiProtonP && particle!=NeutronP && particle!=AntiNeutronP)
    { BuildPhysicsTable(particle); }
  }
}
///////////////////////////////////////////////////////////////
void G4VUserPhysicsList::BuildPhysicsTable(G4ParticleDefinition* particle)
{
  if (fRetrievePhysicsTable) {
#ifdef G4VERBOSE  
    if (verboseLevel>2){
      G4cout << "G4VUserPhysicsList::BuildPhysicsTable  ";
      G4cout << " Retrieve Physics Table for ";
      G4cout << particle->GetParticleName() << G4endl;
    }
#endif
    if ( fIsRestoredCutValues){
      //  Retrieve PhysicsTable from files for proccesses
      RetrievePhysicsTable(particle, directoryPhysicsTable, fStoredInAscii);
      return;
    } else {
#ifdef G4VERBOSE  
      if (verboseLevel>2){
	G4cout << "CheckForRetrievePhysicsTable failed " << G4endl;
      }
#endif
    }
  }

#ifdef G4VERBOSE  
    if (verboseLevel>2){
      G4cout << "G4VUserPhysicsList::BuildPhysicsTable  ";
      G4cout << " for " << particle->GetParticleName() << G4endl;
    }
#endif
  G4int j;
  // Rebuild the physics tables for every process for this particle type
  G4ProcessVector* pVector = (particle->GetProcessManager())->GetProcessList();
  for ( j=0; j < pVector->size(); ++j) {
    (*pVector)[j]->BuildPhysicsTable(*particle);
  }
  for ( j=0; j < pVector->size(); ++j) {
    // temporary addition to make the integral schema
    BuildIntegralPhysicsTable((*pVector)[j], particle); 
  }
}

///////////////////////////////////////////////////////////////
void  G4VUserPhysicsList::BuildIntegralPhysicsTable(G4VProcess* process,
						    G4ParticleDefinition* particle) 
{
  //*********************************************************************
  // temporary addition to make the integral schema of electromagnetic
  // processes work.
  //

  if ( (process->GetProcessName() == "Imsc") ||
       (process->GetProcessName() == "IeIoni") ||
       (process->GetProcessName() == "IeBrems") ||
       (process->GetProcessName() == "Iannihil") ||
       (process->GetProcessName() == "IhIoni") ||
       (process->GetProcessName() == "IMuIoni") ||
       (process->GetProcessName() == "IMuBrems") ||
       (process->GetProcessName() == "IMuPairProd")  ) {
#ifdef G4VERBOSE  
    if (verboseLevel>2){
      G4cout << "G4VUserPhysicsList::BuildIntegralPhysicsTable  ";
      G4cout << " BuildPhysicsTable is invoked for ";
      G4cout << process->GetProcessName();
      G4cout << "(" << particle->GetParticleName() << ")" << G4endl;
    }
#endif
    process->BuildPhysicsTable(*particle);
  }
}

///////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////
//void G4VUserPhysicsList::DumpCutValues(const G4String &particle_name) const
//{
//  G4ParticleDefinition* particle;
//  if ((particle_name == "ALL") || (particle_name == "all")) {
//    theParticleIterator->reset();
//    while( (*theParticleIterator)() ){
//      particle = theParticleIterator->value();
//      DumpCutValues(particle);
//    }
//  } else {
//     particle = theParticleTable->FindParticle(particle_name);
//     if (particle != 0) DumpCutValues(particle);
//  }
//}

///////////////////////////////////////////////////////////////
//void G4VUserPhysicsList::DumpCutValues( G4ParticleDefinition* particle) const
//{
//  if (particle == 0) return;
//  
//  G4int prec = G4cout.precision(3);
//
//  if (particle->IsShortLived()) {
//    // name field
//    G4cout << " --- " << particle->GetParticleName() << " is a short lived particle ------ " << G4endl;
//  } else {
//    // name field
//    G4cout << " --- " << particle->GetParticleName() << " ------ " << G4endl;
//
//    // cut value in range
//    G4double* theRangeCuts = particle->GetLengthCuts();
//
//    // material and energy cut value for the material 
//    G4double*  theKineticEnergyCuts = particle->GetEnergyCuts();
//    
//    if (theRangeCuts != 0) {
//      const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
//      size_t numberOfMaterials = G4Material::GetNumberOfMaterials();
//      G4cout << "   - Material --- Cut in range ------------- Energy Cut ---" << G4endl;
//      for (size_t idx=0; idx<numberOfMaterials; idx++){
//	G4cout << "     " << G4std::setw(19) << (*materialTable)[idx]->GetName(); 
//	G4cout <<" : ";
//	G4cout << G4std::setw(12)<< G4BestUnit(theRangeCuts[idx],"Length");
//	if (theKineticEnergyCuts ==0) {
//	  G4cout << " : --------------------------------";
//	} else {
//	  G4cout << " : ";
//	  G4cout << G4std::setw(10)<< G4BestUnit(theKineticEnergyCuts[idx],"Energy");
//	}	
//	G4cout << G4endl;
//      }
//      if (theKineticEnergyCuts ==0) {
//	G4cout << "   - Cuts in energy are not calculated yet --" << G4endl;
//	G4cout << " Enter /run/initialize command to calculate cuts " << G4endl;
//      }
//    } else {
//      G4cout << "   - Cuts in range are not set yet --" << G4endl;
//    }
//
//  }
//  G4cout.precision(prec);
//}

///////////////////////////////////////////////////////////////
void G4VUserPhysicsList::DumpCutValuesTable(G4int nParticles) 
{ fDisplayThreshold = nParticles; }

void G4VUserPhysicsList::DumpCutValuesTableIfRequested()
{
  if(fDisplayThreshold==0) return;

  // This methods Print out a table of cut value information
  // for "gamma", "e-", "e+", "proton" and "neutron"
  G4int prec = G4cout.precision(3);
  const G4int size = 5;
  G4String name[size] = {"gamma", "e-", "e+", "proton", "neutron"};
  G4ParticleDefinition* particle[size];
  G4int size_display=fDisplayThreshold;
  if(size_display>size) size_display=size;
  G4int idx; 
  for (idx=0; idx <size_display; idx++) {
    particle[idx] = theParticleTable->FindParticle(name[idx]);
  }

  G4cout << G4endl;
  G4cout << "======== Table of production thresholds in range and energy ===========" << G4endl;

  typedef G4std::vector<G4Region*>::iterator regionIterator;
  G4RegionStore* fG4RegionStore = G4RegionStore::GetInstance();
  for(regionIterator rItr=fG4RegionStore->begin();rItr!=fG4RegionStore->end();rItr++)
  {
    G4cout << G4endl;
    G4cout << "-----------------------------------------------------------------------" << G4endl;
    G4cout << "Region : " << (*rItr)->GetName() << G4endl;
    G4cout << "-----------------------------------------------------------------------" << G4endl;

  
    G4cout << "                             ";
    for (idx=0; idx <size_display; idx++) {
      G4cout << " " << G4std::setw(11) << name[idx] << "    ";
    }
    G4cout << G4endl;

    const G4ProductionCuts* aCut = (*rItr)->GetProductionCuts();
    G4cout << "Threshold in range ";
    for (idx=0; idx <size_display; idx++) {
      if(aCut->GetProductionCut(name[idx])<=0.0)
      { G4cout << " ---------- "; }
      else
      { G4cout << " " << G4std::setw(11) << G4BestUnit(aCut->GetProductionCut(name[idx]),"Length"); }
    }
    G4cout << G4endl;

    G4std::vector<G4Material*>::const_iterator mItr = (*rItr)->GetMaterialIterator();
    size_t nMaterial = (*rItr)->GetNumberOfMaterials();
    for (size_t I=0; I<nMaterial; I++) {
      G4cout << " " << G4std::setw(19) << (*mItr)->GetName();
      for (idx=0; idx <size_display; idx++) {
        G4int coupleIndex = G4ProductionCutsTable::GetProductionCutsTable()->GetCoupleIndex((*mItr),aCut);
        if(coupleIndex<0)
        { G4cout << " ---------- "; }
        else
        { G4cout << " " << G4std::setw(11) << G4BestUnit((particle[idx]->GetEnergyCuts())[coupleIndex],"Energy"); }
      }
      G4cout << G4endl;
      mItr++;
    }
  }
 
  G4cout << "=======================================================================" << G4endl;
  G4cout << G4endl;

  G4cout.precision(prec);
  //fDisplayThreshold=0;
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
G4bool G4VUserPhysicsList::StorePhysicsTable(const G4String& directory)
{
  G4bool   ascii = fStoredInAscii;
  G4String dir   = directory;
  if (dir.isNull()) dir = directoryPhysicsTable; 
  else directoryPhysicsTable = dir; 
  
  // store CutsTable info
  if (!fCutsTable->StoreCutsTable(dir, ascii)) return false;
#ifdef G4VERBOSE  
  if (verboseLevel>2){
    G4cout << "G4VUserPhysicsList::StorePhysicsTable   ";
    G4cout << " Store material and cut values successfully" << G4endl;
  }
#endif

  G4bool success= true;  

  // loop over all particles in G4ParticleTable 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    // Store physics tables for every process for this particle type
    G4ProcessVector* pVector = (particle->GetProcessManager())->GetProcessList();
    G4int  j;
    for ( j=0; j < pVector->size(); ++j) {
      if (!(*pVector)[j]->StorePhysicsTable(particle,dir,ascii)){   
#ifdef G4VERBOSE  
	if (verboseLevel>2){
	  G4cout << "G4VUserPhysicsList::StorePhysicsTable   ";
	  G4cout << " Fail to store for ";
	  G4cout << (*pVector)[j]->GetProcessName();
	  G4cout << "(" << particle->GetParticleName() <<")" << G4endl;
	}
#endif
	success = false;
      }
    }
    // end loop over processes
  }
  // end loop over particles
  return success;
}



///////////////////////////////////////////////////////////////
 void  G4VUserPhysicsList::SetPhysicsTableRetrieved(const G4String& directory)
{
  fRetrievePhysicsTable = true;
  if(!directory.isNull()) {
    directoryPhysicsTable = directory;
  }
  fIsCheckedForRetrievePhysicsTable=false;
  fIsRestoredCutValues = false;
  ResetCuts();
}

///////////////////////////////////////////////////////////////
void G4VUserPhysicsList::RetrievePhysicsTable(G4ParticleDefinition* particle, 
					      const G4String& directory,
					      G4bool          ascii)
{
  G4int  j;
  G4bool success[100];  
  // Retrieve physics tables for every process for this particle type
  G4ProcessVector* pVector = (particle->GetProcessManager())->GetProcessList();
  for ( j=0; j < pVector->size(); ++j) {
    success[j] = 
       (*pVector)[j]->RetrievePhysicsTable(particle,directory,ascii);

    if (!success[j]) {
#ifdef G4VERBOSE  
      if (verboseLevel>2){
	G4cout << "G4VUserPhysicsList::RetrievePhysicsTable   ";
	G4cout << " Fail to retrieve for ";
        G4cout << (*pVector)[j]->GetProcessName();
        G4cout << "(" << particle->GetParticleName() <<")" << G4endl;
      }
#endif
      (*pVector)[j]->BuildPhysicsTable(*particle);
    }
  }
  for ( j=0; j < pVector->size(); ++j) {
    // temporary addition to make the integral schema
    if (!success[j]) BuildIntegralPhysicsTable((*pVector)[j], particle); 
  }
}






