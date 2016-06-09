//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
//	History
//       first version                    09 Jan 1998 by H.Kurashige
//       Added SetEnergyRange             18 Jun 1998 by H.Kurashige 
//       Change for short lived particles 27 Jun 1998 by H.Kurashige
//       G4BestUnit on output             12 nov 1998 by M.Maire
//       Added RemoveProcessManager        9 Feb 1999 by H.Kurashige
//       Fixed RemoveProcessManager       15 Apr 1999 by H.Kurashige
//       Removed ConstructAllParticles()  15 Apr 1999 by H.Kurashige
//       Modified for CUTS per REGION     10 Oct 2002 by H.Kurashige
//       Check if particle IsShortLived   18 Jun 2003 by V.Ivanchenko
//       Modify PreparePhysicsList        18 Jan 2006 by H.Kurashige
//       Added PhysicsListHelper           29 APr. 2011 H.Kurashige
//       Added default impelmentation of SetCuts 10 June 2011 H.Kurashige 
//           SetCuts is not 'pure virtual' any more 
// ------------------------------------------------------------

#include "G4VUserPhysicsList.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"
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
#include <iomanip>
#include <fstream>

////////////////////////////////////////////////////////
G4VUserPhysicsList::G4VUserPhysicsList()
  :verboseLevel(1),
   defaultCutValue(1.0 * mm),
   isSetDefaultCutValue(false),
   fRetrievePhysicsTable(false),
   fStoredInAscii(true),
   fIsCheckedForRetrievePhysicsTable(false),
   fIsRestoredCutValues(false),
   directoryPhysicsTable("."),
   fDisplayThreshold(0),
   fIsPhysicsTableBuilt(false),
   fDisableCheckParticleList(false)
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
 
  // PhysicsListHelper
  thePLHelper = G4PhysicsListHelper::GetPhysicsListHelper();
  thePLHelper->SetVerboseLevel(verboseLevel);

}

////////////////////////////////////////////////////////
G4VUserPhysicsList::~G4VUserPhysicsList()
{
  if (theMessenger != 0) {
    delete theMessenger;
    theMessenger = 0;
  }
  RemoveProcessManager();

  // invoke DeleteAllParticle
  theParticleTable->DeleteAllParticles();

}

////////////////////////////////////////////////////////
G4VUserPhysicsList::G4VUserPhysicsList(const G4VUserPhysicsList& right)
  :verboseLevel(right.verboseLevel),
   defaultCutValue(right.defaultCutValue),
   isSetDefaultCutValue(right.isSetDefaultCutValue),
   fRetrievePhysicsTable(right.fRetrievePhysicsTable),
   fStoredInAscii(right.fStoredInAscii),
   fIsCheckedForRetrievePhysicsTable(right.fIsCheckedForRetrievePhysicsTable),
   fIsRestoredCutValues(right.fIsRestoredCutValues),
   directoryPhysicsTable(right.directoryPhysicsTable),
   fDisplayThreshold(right.fDisplayThreshold),
   fIsPhysicsTableBuilt(right.fIsPhysicsTableBuilt),
   fDisableCheckParticleList(right.fDisableCheckParticleList)
{
  // pointer to the particle table
  theParticleTable = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();

  // pointer to the cuts table
  fCutsTable =  G4ProductionCutsTable::GetProductionCutsTable();

  // UI Messenger
  theMessenger = new G4UserPhysicsListMessenger(this);
 
  // PhysicsListHelper
  thePLHelper = G4PhysicsListHelper::GetPhysicsListHelper();
  thePLHelper->SetVerboseLevel(verboseLevel);

}


////////////////////////////////////////////////////////
G4VUserPhysicsList & G4VUserPhysicsList::operator=(const G4VUserPhysicsList & right)
{
  if (this != &right) {
    verboseLevel   = right.verboseLevel;
    defaultCutValue = right.defaultCutValue;
    isSetDefaultCutValue = right.isSetDefaultCutValue;
    fRetrievePhysicsTable = right.fRetrievePhysicsTable;
    fStoredInAscii = right.fStoredInAscii;
    fIsCheckedForRetrievePhysicsTable = right.fIsCheckedForRetrievePhysicsTable;
    fIsRestoredCutValues = right.fIsRestoredCutValues;
    directoryPhysicsTable = right.directoryPhysicsTable;
    fDisplayThreshold = right.fDisplayThreshold;
    fIsPhysicsTableBuilt = right.fIsPhysicsTableBuilt;
    fDisableCheckParticleList = right.fDisableCheckParticleList;
  }
  return *this;
}

////////////////////////////////////////////////////////
void G4VUserPhysicsList::AddProcessManager(G4ParticleDefinition* newParticle,
					   G4ProcessManager*     newManager)
{
  if (newParticle == 0) return;
  if (newParticle->GetProcessManager() != 0) {
#ifdef G4VERBOSE
    if (verboseLevel >1){
      G4cout << "G4VUserPhysicsList::AddProcessManager: "
	     << newParticle->GetParticleName()
	     << " already has ProcessManager " << G4endl;
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
	  G4Exception("G4VUserPhysicsList::AddProcessManager",
		      "Run0251", RunMustBeAborted,
		      "GenericIon has no ProcessMamanger"); 	
	}
      } else {
	// "GenericIon" does not exist
	newManager = new G4ProcessManager(newParticle);
	G4Exception("G4VUserPhysicsList::AddProcessManager",
		    "Run0252", RunMustBeAborted,
		    "GenericIon does not exist"); 	
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
    G4cout << "G4VUserPhysicsList::AddProcessManager: "
	   << "adds ProcessManager to "
	   << newParticle->GetParticleName() << G4endl;
    newManager->DumpInfo();
  }
#endif
  if ( fIsPhysicsTableBuilt
       && (newParticle->GetParticleType() == "nucleus")) {
    PreparePhysicsTable(newParticle);
    BuildPhysicsTable(newParticle);
  }
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

/////////////////////////////////////////////////////////
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
      G4cout << "G4VUserPhysicsList::RemoveProcessManager: "
	     << "remove ProcessManager from "
	     << particle->GetParticleName() << G4endl;
    }
#endif
  }
}


////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetCuts()
{
  if ( !isSetDefaultCutValue ){
    SetDefaultCutValue(defaultCutValue);
  }

#ifdef G4VERBOSE
  if (verboseLevel >1){
    G4cout << "G4VUserPhysicsList::SetCuts:   " << G4endl;
    G4cout << "Cut for gamma: " << GetCutValue("gamma")/mm  
	   << "[mm]" << G4endl;
    G4cout << "Cut  for e-: " << GetCutValue("e-")/mm  
	   << "[mm]" << G4endl;
    G4cout << "Cut  for e+: " << GetCutValue("e+")/mm 
	   << "[mm]" << G4endl;
    G4cout << "Cut  for proton: " << GetCutValue("proton")/mm  
	   << "[mm]" << G4endl;
  }
#endif

  // dump Cut values if verboseLevel==3
  if (verboseLevel>2) {
    DumpCutValuesTable();
  }
}


////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetDefaultCutValue(G4double value)
{
  if (value<0.0) {
#ifdef G4VERBOSE
    if (verboseLevel >0){
      G4cout << "G4VUserPhysicsList::SetDefaultCutValue: negative cut values"
	     << "  :" << value/mm << "[mm]" << G4endl;
    }
#endif
    return;
  }

  defaultCutValue = value;
  isSetDefaultCutValue = true;
  
  // set cut values for gamma at first and for e- and e+
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  SetCutValue(defaultCutValue, "proton");    

#ifdef G4VERBOSE
  if (verboseLevel >1){
    G4cout << "G4VUserPhysicsList::SetDefaultCutValue:"
	   << "default cut value is changed to   :" 
	   << defaultCutValue/mm << "[mm]" << G4endl;
  }
#endif
 }


////////////////////////////////////////////////////////
G4double G4VUserPhysicsList::GetCutValue(const G4String& name) const
{
  size_t nReg = (G4RegionStore::GetInstance())->size();
  if (nReg==0) {
#ifdef G4VERBOSE
    if (verboseLevel>0){      
      G4cout << "G4VUserPhysicsList::GetCutValue "
	     <<" : No Default Region " <<G4endl;
    }
#endif
    G4Exception("G4VUserPhysicsList::GetCutValue",
		"Run0253", FatalException,
		"No Default Region");
    return -1.*mm;
  }
  G4Region* region = (*(G4RegionStore::GetInstance()))[0];
  return region->GetProductionCuts()->GetProductionCut(name);
}

////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetCutValue(G4double aCut, const G4String& name)
{
  SetParticleCuts( aCut ,name );
}

////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetCutValue
(G4double aCut, const G4String& pname, const G4String& rname)
{
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(rname);
  if (region != 0){
    //set cut value
    SetParticleCuts( aCut ,pname, region );
  } else {
#ifdef G4VERBOSE
    if (verboseLevel>0){      
      G4cout << "G4VUserPhysicsList::SetCutValue "
	     <<" : No Region of " << rname << G4endl;
    }
#endif
  }
}


////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetCutsWithDefault()
{
  SetDefaultCutValue(defaultCutValue);
  G4VUserPhysicsList::SetCuts();
}

////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetCutsForRegion(G4double aCut, const G4String& rname)
{
  // set cut values for gamma at first and for e- and e+
  SetCutValue(aCut, "gamma", rname);
  SetCutValue(aCut, "e-", rname);
  SetCutValue(aCut, "e+", rname);
  SetCutValue(aCut, "proton", rname);
}



////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetParticleCuts( G4double cut, G4ParticleDefinition* particle, G4Region* region)
{
  SetParticleCuts(cut, particle->GetParticleName(), region);
}

////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetParticleCuts( G4double cut, const G4String& particleName, G4Region* region)
{
  if (cut<0.0) {
#ifdef G4VERBOSE
    if (verboseLevel >0){
      G4cout << "G4VUserPhysicsList::SetParticleCuts: negative cut values"
	     << "  :" << cut/mm << "[mm]" 
	     << " for "<< particleName << G4endl;
    }
#endif
    return;
  }

  if(!region){
    size_t nReg = (G4RegionStore::GetInstance())->size();
    if (nReg==0) {
#ifdef G4VERBOSE
      if (verboseLevel>0){      
	G4cout << "G4VUserPhysicsList::SetParticleCuts "
	       <<" : No Default Region " <<G4endl;
      }
#endif
      G4Exception("G4VUserPhysicsList::SetParticleCuts ",
		  "Run0254", FatalException,
		"No Default Region");
      return;
    }
    region = (*(G4RegionStore::GetInstance()))[0];
  }

  if ( !isSetDefaultCutValue ){
    SetDefaultCutValue(defaultCutValue);
  }

  G4ProductionCuts* pcuts = region->GetProductionCuts();
  pcuts->SetProductionCut(cut,particleName);
#ifdef G4VERBOSE
  if (verboseLevel>2){      
    G4cout << "G4VUserPhysicsList::SetParticleCuts: "
	   << "  :" << cut/mm << "[mm]" 
	   << " for "<< particleName << G4endl;
  }
#endif
}

///////////////////////////////////////////////////////////////
void G4VUserPhysicsList::BuildPhysicsTable()
{
  //Prepare Physics table for all particles 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    PreparePhysicsTable(particle); 
  }

  // ask processes to prepare physics table 
  if (fRetrievePhysicsTable) {
    fIsRestoredCutValues = fCutsTable->RetrieveCutsTable(directoryPhysicsTable, fStoredInAscii);
    // check if retrieve Cut Table successfully
    if (!fIsRestoredCutValues) {
#ifdef G4VERBOSE
      if (verboseLevel>0){
	G4cout << "G4VUserPhysicsList::BuildPhysicsTable"
	       << " Retrieve Cut Table failed !!" << G4endl;
      }
#endif	
      G4Exception("G4VUserPhysicsList::BuildPhysicsTable",
		  "Run0255", RunMustBeAborted,
		  "Fail to retrieve Production Cut Table");
    } else {
#ifdef G4VERBOSE
      if (verboseLevel>2){
	G4cout << "G4VUserPhysicsList::BuildPhysicsTable"
	       << "  Retrieve Cut Table successfully " << G4endl;
      }
#endif 
    }     
  } else {
#ifdef G4VERBOSE
    if (verboseLevel>2){
      G4cout << "G4VUserPhysicsList::BuildPhysicsTable"
	     << " does not retrieve Cut Table but calculate " << G4endl;
    } 
#endif	    
  }

  // Sets a value to particle
  // set cut values for gamma at first and for e- and e+
  G4String particleName;
  G4ParticleDefinition* GammaP = theParticleTable->FindParticle("gamma");
  if(GammaP) BuildPhysicsTable(GammaP);
  G4ParticleDefinition* EMinusP = theParticleTable->FindParticle("e-");
  if(EMinusP) BuildPhysicsTable(EMinusP);
  G4ParticleDefinition* EPlusP = theParticleTable->FindParticle("e+");
  if(EPlusP) BuildPhysicsTable(EPlusP);
  G4ParticleDefinition* ProtonP = theParticleTable->FindParticle("proton");
  if(ProtonP) BuildPhysicsTable(ProtonP);


  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    if( particle!=GammaP && 
	particle!=EMinusP && 
	particle!=EPlusP && 
	particle!=ProtonP   ){
      BuildPhysicsTable(particle); 
    }
  }

  // Set flag
  fIsPhysicsTableBuilt = true;

}
///////////////////////////////////////////////////////////////
void G4VUserPhysicsList::BuildPhysicsTable(G4ParticleDefinition* particle)
{
  if (fRetrievePhysicsTable) {
    if ( !fIsRestoredCutValues){
      // fail to retreive cut tables
#ifdef G4VERBOSE
      if (verboseLevel>0){
	G4cout << "G4VUserPhysicsList::BuildPhysicsTable  "
	       << "Physics table can not be retreived and will be calculated "
	       << G4endl;
      }
#endif
      fRetrievePhysicsTable = false; 

    } else {
#ifdef G4VERBOSE
      if (verboseLevel>2){
	G4cout << "G4VUserPhysicsList::BuildPhysicsTable  "
	       << " Retrieve Physics Table for "
	       << particle->GetParticleName() << G4endl;
      }
#endif
      //  Retrieve PhysicsTable from files for proccesses
      RetrievePhysicsTable(particle, directoryPhysicsTable, fStoredInAscii);
    }
  }

#ifdef G4VERBOSE
  if (verboseLevel>2){
    G4cout << "G4VUserPhysicsList::BuildPhysicsTable  "
	   << "Calculate Physics Table for " 
	   << particle->GetParticleName() << G4endl;
  }
#endif
  // Rebuild the physics tables for every process for this particle type
  // if particle is not ShortLived
  if(!particle->IsShortLived()) {
    G4ProcessManager* pManager =  particle->GetProcessManager();
    if (!pManager) {
#ifdef G4VERBOSE
      if (verboseLevel>0){
	G4cout << "G4VUserPhysicsList::BuildPhysicsTable "
	       <<" : No Process Manager for "
	       << particle->GetParticleName() << G4endl;
	G4cout << particle->GetParticleName() 
	       << " should be created in your PhysicsList" <<G4endl;
      }
#endif
      G4Exception("G4VUserPhysicsList::BuildPhysicsTable",
		  "Run0271", FatalException,  
		  "No process manager");
      return;
    }
    G4ProcessVector* pVector = pManager->GetProcessList();
    if (!pVector) {
#ifdef G4VERBOSE
      if (verboseLevel>0){      
	G4cout << "G4VUserPhysicsList::BuildPhysicsTable  "
	       <<" : No Process Vector for " 
	       << particle->GetParticleName() <<G4endl;
      }
#endif
      G4Exception("G4VUserPhysicsList::BuildPhysicsTable",
		  "Run0272", FatalException,
		  "No process Vector");
      return;
    }
    for (G4int j=0; j < pVector->size(); ++j) {
      (*pVector)[j]->BuildPhysicsTable(*particle);
    }
  }
}

///////////////////////////////////////////////////////////////
void G4VUserPhysicsList::PreparePhysicsTable(G4ParticleDefinition* particle)
{
  // Prepare the physics tables for every process for this particle type
  // if particle is not ShortLived
  if(!particle->IsShortLived()) {
    G4ProcessManager* pManager =  particle->GetProcessManager();
    if (!pManager) {
#ifdef G4VERBOSE
      if (verboseLevel>0) {
	G4cout<< "G4VUserPhysicsList::PreparePhysicsTable  "
	      << ": No Process Manager for " 
	      << particle->GetParticleName() <<G4endl;
	G4cout << particle->GetParticleName() 
	       << " should be created in your PhysicsList" <<G4endl;
      }
#endif
      G4Exception("G4VUserPhysicsList::PreparePhysicsTable",
		  "Run0273", FatalException, 
		  "No process manager");
      return;
    }
    
    G4ProcessVector* pVector = pManager->GetProcessList();
    if (!pVector) {
#ifdef G4VERBOSE
      if (verboseLevel>0) {
	G4cout << "G4VUserPhysicsList::PreparePhysicsTable  "
	       << ": No Process Vector for " 
	       << particle->GetParticleName() <<G4endl;
      }
#endif
      G4Exception("G4VUserPhysicsList::PreparePhysicsTable",
		  "Run0274", FatalException,
		  "No process Vector");
      return;
    }
    for (G4int j=0; j < pVector->size(); ++j) {
      (*pVector)[j]->PreparePhysicsTable(*particle);
    }
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
      G4cout << "G4VUserPhysicsList::BuildIntegralPhysicsTable  "
	     << " BuildPhysicsTable is invoked for "
	     << process->GetProcessName()
	     << "(" << particle->GetParticleName() << ")" << G4endl;
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
void G4VUserPhysicsList::DumpCutValuesTable(G4int flag) 
{ 
  fDisplayThreshold = flag; 
}

///////////////////////////////////////////////////////////////
void G4VUserPhysicsList::DumpCutValuesTableIfRequested()
{
  if(fDisplayThreshold==0) return;
  G4ProductionCutsTable::GetProductionCutsTable()->DumpCouples();
  fDisplayThreshold = 0;
}


///////////////////////////////////////////////////////////////
G4bool G4VUserPhysicsList::StorePhysicsTable(const G4String& directory)
{
  G4bool   ascii = fStoredInAscii;
  G4String dir   = directory;
  if (dir.isNull()) dir = directoryPhysicsTable; 
  else directoryPhysicsTable = dir; 

  // store CutsTable info
  if (!fCutsTable->StoreCutsTable(dir, ascii)) {
    G4Exception("G4VUserPhysicsList::StorePhysicsTable",
		"Run0281", JustWarning,
		"Fail to store Cut Table"); 	
    return false;
  }
#ifdef G4VERBOSE  
  if (verboseLevel>2){
    G4cout << "G4VUserPhysicsList::StorePhysicsTable   "
	   << " Store material and cut values successfully" << G4endl;
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
	G4String comment =  "Fail to store physics table for "; 
	comment += (*pVector)[j]->GetProcessName();
	comment += "(" + particle->GetParticleName()  + ")";
	G4Exception("G4VUserPhysicsList::StorePhysicsTable",
		    "Run0282", JustWarning,
		    comment); 	
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
	G4cout << "G4VUserPhysicsList::RetrievePhysicsTable   "
	       << " Fail to retrieve Physics Table for "
	       << (*pVector)[j]->GetProcessName() << G4endl;
	G4cout << "Calculate Physics Table for " 
	       << particle->GetParticleName() << G4endl;
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


///////////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetApplyCuts(G4bool value, const G4String& name)
{
#ifdef G4VERBOSE  
  if (verboseLevel>2){
    G4cout << "G4VUserPhysicsList::SetApplyCuts for " << name << G4endl;
  }
#endif
  if(name=="all") {
    theParticleTable->FindParticle("gamma")->SetApplyCutsFlag(value);
    theParticleTable->FindParticle("e-")->SetApplyCutsFlag(value);
    theParticleTable->FindParticle("e+")->SetApplyCutsFlag(value);
    theParticleTable->FindParticle("proton")->SetApplyCutsFlag(value);
  } else {
    theParticleTable->FindParticle(name)->SetApplyCutsFlag(value);
  }
}

///////////////////////////////////////////////////////////////
G4bool G4VUserPhysicsList::GetApplyCuts(const G4String& name) const
{
  return theParticleTable->FindParticle(name)->GetApplyCutsFlag();
}


////////////////////////////////////////////////////////
void G4VUserPhysicsList::CheckParticleList()
{
  if (! fDisableCheckParticleList ){
    thePLHelper->CheckParticleList();
  }
}

////////////////////////////////////////////////////////
void G4VUserPhysicsList::AddTransportation()
{   
  thePLHelper->AddTransportation();
}

////////////////////////////////////////////////////////
void G4VUserPhysicsList::UseCoupledTransportation(G4bool vl)
{ 
  thePLHelper->UseCoupledTransportation(vl);
}

////////////////////////////////////////////////////////
G4bool G4VUserPhysicsList::RegisterProcess(G4VProcess*            process,
					  G4ParticleDefinition*  particle)
{
  return thePLHelper->RegisterProcess(process, particle);
}

////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
  // set verboseLevel for G4ProductionCutsTable same as one for G4VUserPhysicsList: 
  fCutsTable->SetVerboseLevel(verboseLevel);

  thePLHelper->SetVerboseLevel(verboseLevel);

#ifdef G4VERBOSE
  if (verboseLevel >1){
    G4cout << "G4VUserPhysicsList::SetVerboseLevel  :"
	   << " Verbose level is set to " << verboseLevel << G4endl;
  }
#endif
}


///////////////////////////////////////////////////////////////
/// obsolete methods

///////////////////////////////////////////////////////////////
void G4VUserPhysicsList::ResetCuts()
{
#ifdef G4VERBOSE  
  if (verboseLevel>0){
    G4cout << "G4VUserPhysicsList::ResetCuts() is obsolete."
	   << " This method gives no effect and you can remove it. "<< G4endl;
  }
#endif
}
