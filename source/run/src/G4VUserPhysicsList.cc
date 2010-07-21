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
// $Id: G4VUserPhysicsList.cc,v 1.72 2010-07-21 14:21:19 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// ------------------------------------------------------------

#include "globals.hh"
#include "G4VUserPhysicsList.hh"
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
                 :  fDisableCheckParticleList(false),
		    verboseLevel(1),
		    fRetrievePhysicsTable(false),
		    fStoredInAscii(true),
		    fIsCheckedForRetrievePhysicsTable(false),
		    fIsRestoredCutValues(false),
                    directoryPhysicsTable("."),
                    fDisplayThreshold(0),
		    fIsPhysicsTableBuilt(false),
                    useCoupledTransportation(false)
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
  RemoveProcessManager();

  // invoke DeleteAllParticle
  theParticleTable->DeleteAllParticles();

}

////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
  // set verboseLevel for G4ProductionCutsTable same as one for G4VUserPhysicsList: 
  fCutsTable->SetVerboseLevel(value);

#ifdef G4VERBOSE
  if (verboseLevel >1){
    G4cout << "G4VUserPhysicsList::SetVerboseLevel  :";
    G4cout << " Verbose level is set to " << verboseLevel << G4endl;
  }
#endif
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
	  G4Exception("G4VUserPhysicsList::AddProcessManager","Error in GenericIon",
		RunMustBeAborted,"GenericIon has no ProcessMamanger"); 	
	}
      } else {
	// "GenericIon" does not exist
	newManager = new G4ProcessManager(newParticle);
	G4Exception("G4VUserPhysicsList::AddProcessManager","No GenericIon",
		    RunMustBeAborted,"GenericIon does not exist"); 	
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
  if ( fIsPhysicsTableBuilt
       && (newParticle->GetParticleType() == "nucleus")) {
    PreparePhysicsTable(newParticle);
    BuildPhysicsTable(newParticle);
  }
}


////////////////////////////////////////////////////////
void G4VUserPhysicsList::CheckParticleList()
{

  // skip if fDisableCheckParticleList is set  
  if (fDisableCheckParticleList) return;

  bool isElectron = false;
  bool isPositron = false;
  bool isGamma    = false;
  bool isProton   = false;
  bool isGenericIon = false;
  bool isAnyIon   = false;
  bool isAnyChargedBaryon   = false;
  bool isEmProc   = false;

  // loop over all particles in G4ParticleTable
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String name = particle->GetParticleName();
    // check if any EM process exists
    if (!isEmProc) {
      G4ProcessVector* list = particle->GetProcessManager()->GetProcessList();
      for (int idx=0; idx<list->size(); idx++){
	isEmProc = ((*list)[idx])->GetProcessType() == fElectromagnetic;
	if (isEmProc) break;
      }
    }
    
    if      ( name == "e-") isElectron = true; 
    else if ( name == "e+") isPositron = true; 
    else if ( name == "gamma") isGamma = true; 
    else if ( name == "GenericIon") isGenericIon = true; 
    else if ( name == "proton") isProton = true; 
    else if ( particle->GetParticleType() == "nucleus") isAnyIon = true;
    else if ( particle->GetParticleType() == "baryon") {
       if ( particle->GetPDGCharge() != 0.0 ) isAnyChargedBaryon = true;
    }
  }

  if (!isEmProc) return;

  // RULE 1
  //  e+, e- and gamma should exist 
  //   if one of them exist
  bool isEmBasic =  isElectron || isPositron || isGamma;
  bool isMissingEmBasic =  !isElectron || !isPositron || !isGamma;
  if (isEmBasic && isMissingEmBasic) {
    G4String missingName="";
    if (!isElectron) missingName += "e- ";
    if (!isPositron) missingName += "e+ ";
    if (!isGamma) missingName += "gamma ";

#ifdef G4VERBOSE
    if (verboseLevel >0){
      G4cout << "G4VUserPhysicsList::CheckParticleList: ";
      G4cout << missingName << " do not exist " << G4endl; 
      G4cout << " These particle are necessary for basic EM processes" << G4endl;
    }
#endif
    missingName += " should be created ";
    G4Exception("G4VUserPhysicsList::CheckParticleList","Missing EM basic particle",
		FatalException, missingName); 	
  }

  // RULE 2
  //  proton should exist 
  //   if any other charged baryon  exist
  if (!isProton && isAnyChargedBaryon) {
    G4String missingName="proton ";

#ifdef G4VERBOSE
    if (verboseLevel >0){
      G4cout << "G4VUserPhysicsList::CheckParticleList: ";
      G4cout << missingName << " does not exist "<< G4endl; 
      G4cout << " Proton is necessary for EM baryon processes" << G4endl;
    }
#endif
    missingName += " should be created ";
    G4Exception("G4VUserPhysicsList::CheckParticleList","Missing Proton",
		FatalException, missingName); 	
  }
   
  // RULE 3
  //  GenericIonn should exist 
  //   if any other ion  exist
  if (!isGenericIon && isAnyIon) {
    G4String missingName="GenericIon ";

#ifdef G4VERBOSE
    if (verboseLevel >0){
      G4cout << "G4VUserPhysicsList::CheckParticleList: ";
      G4cout << missingName << " does not exist "<< G4endl; 
      G4cout << " GenericIon should be created if any ion is necessary" << G4endl;
    }
#endif
    missingName += " should be created ";
    G4Exception("G4VUserPhysicsList::CheckParticleList","Missing GenericIon",
		FatalException, missingName); 	
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
      G4cout << "G4VUserPhysicsList::RemoveProcessManager: ";
      G4cout  << "remove ProcessManager from ";
      G4cout  << particle->GetParticleName() << G4endl;
    }
#endif
  }
}


////////////////////////////////////////////////////////
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"
#include "G4RunManagerKernel.hh"
#include "G4ScoringManager.hh"

void G4VUserPhysicsList::AddTransportation()
{
  G4int verboseLevelTransport = 0;
  G4VProcess* theTransportationProcess = 0;  // Pointer ownership handled to
                                             // G4ProcessManager !
  G4int nParaWorld = G4RunManagerKernel::GetRunManagerKernel()->GetNumberOfParallelWorld();

  if(nParaWorld || useCoupledTransportation || G4ScoringManager::GetScoringManagerIfExist())
    {
      theTransportationProcess = new G4CoupledTransportation(verboseLevelTransport);
      G4cout << "#############################################################################" << G4endl
             << " G4VUserPhysicsList::AddTransportation() --- G4CoupledTransportation is used " << G4endl
             << "#############################################################################" << G4endl;
    }
  else
    {
      theTransportationProcess = new G4Transportation(verboseLevelTransport);
    }
 
#ifdef G4VERBOSE
    if (verboseLevel >2){
      G4cout << "G4VUserPhysicsList::AddTransportation()  "<< G4endl;
    }
#endif

  // loop over all particles in G4ParticleTable
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (!particle->IsShortLived()) {
      // Add transportation process for all particles other than  "shortlived"
      if ( pmanager == 0) {
	// Error !! no process manager
	G4String particleName = particle->GetParticleName();
	G4Exception("G4VUserPhysicsList::AddTransportation","No process manager",
		    FatalException, particleName );
      } else {
	// add transportation with ordering = ( -1, "first", "first" )
	pmanager ->AddProcess(theTransportationProcess);
	pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
	pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
      }
    } else {
      // shortlived particle case
     // Add transportation process for all particles other than  "shortlived"
      if ( pmanager == 0) {
        // Error !! no process manager
        G4String particleName = particle->GetParticleName();
        G4Exception("G4VUserPhysicsList::AddTransportation","No process manager",
                    FatalException, particleName );
      } else {
        // add transportation with ordering = ( -1, "first", "first" )
        pmanager ->AddProcess(theTransportationProcess);
        pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
        pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
      }

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

   }
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
void G4VUserPhysicsList::SetCutValue
(G4double aCut, const G4String& pname, const G4String& rname)
{
  G4ParticleDefinition* particle = theParticleTable->FindParticle(pname);
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(rname);
  if (particle != 0 && region != 0){
    if (!particle->IsShortLived()) {
      //set cut value
      SetParticleCuts( aCut ,particle, region );
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
  SetCutValue(cut, "proton");

  // dump Cut values if verboseLevel==3
  if (verboseLevel>2) {
    DumpCutValuesTable();
  }
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
////////////////////////////////////////////////////////
void G4VUserPhysicsList::SetParticleCuts( G4double cut, G4ParticleDefinition* particle, G4Region* region)
{
  if(!region) region = (*(G4RegionStore::GetInstance()))[0];
  G4ProductionCuts* pcuts = region->GetProductionCuts();
  pcuts->SetProductionCut(cut,particle);
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
	G4cout << "G4VUserPhysicsList::BuildPhysicsTable";
	G4cout << " Retrieve Cut Table failed !!" << G4endl;
      }
#endif	
      G4Exception("G4VUserPhysicsList::BuildPhysicsTable","Fail to Retreive",
		  RunMustBeAborted,"Production Cut Table can not be retreived");
    } else {
#ifdef G4VERBOSE
      if (verboseLevel>2){
	G4cout << "G4VUserPhysicsList::BuildPhysicsTable";
	G4cout << "  Retrieve Cut Table successfully " << G4endl;
      }
#endif 
    }     
  } else {
#ifdef G4VERBOSE
    if (verboseLevel>2){
      G4cout << "G4VUserPhysicsList::BuildPhysicsTable";
      G4cout << " does not retrieve Cut Table but calculate " << G4endl;
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
	G4cout << "G4VUserPhysicsList::BuildPhysicsTable  ";
	G4cout << "Physics table can not be retreived and will be calculated "<< G4endl;
      }
#endif
      fRetrievePhysicsTable = false; 

    } else {
#ifdef G4VERBOSE
      if (verboseLevel>2){
	G4cout << "G4VUserPhysicsList::BuildPhysicsTable  ";
	G4cout << " Retrieve Physics Table for ";
	G4cout << particle->GetParticleName() << G4endl;
      }
#endif
      //  Retrieve PhysicsTable from files for proccesses
      RetrievePhysicsTable(particle, directoryPhysicsTable, fStoredInAscii);
    }
  }

#ifdef G4VERBOSE
  if (verboseLevel>2){
    G4cout << "G4VUserPhysicsList::BuildPhysicsTable  ";
    G4cout << "Calculate Physics Table for " << particle->GetParticleName() << G4endl;
  }
#endif
  // Rebuild the physics tables for every process for this particle type
  // if particle is not ShortLived
  if(!particle->IsShortLived()) {
    G4ProcessManager* pManager =  particle->GetProcessManager();
    if (!pManager) {
      G4cerr << "G4VUserPhysicsList::BuildPhysicsTable  : No Process Manager for " 
             << particle->GetParticleName() <<G4endl;
      G4cerr << particle->GetParticleName() << " should be created in your PhysicsList" <<G4endl;
      G4Exception("G4VUserPhysicsList::BuildPhysicsTable","No process manager",
                    FatalException,  particle->GetParticleName() );
      return;
    }
    G4ProcessVector* pVector = pManager->GetProcessList();
    if (!pVector) {
      G4cerr << "G4VUserPhysicsList::BuildPhysicsTable  : No Process Vector for " 
             << particle->GetParticleName() <<G4endl;
      G4cerr << particle->GetParticleName() << " should be created in your PhysicsList" <<G4endl;
      G4Exception("G4VUserPhysicsList::BuildPhysicsTable","No process Vector",
                    FatalException,  particle->GetParticleName() );
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
      G4cerr << "G4VUserPhysicsList::PreparePhysicsTable  : No Process Manager for " 
             << particle->GetParticleName() <<G4endl;
      G4cerr << particle->GetParticleName() << " should be created in your PhysicsList" <<G4endl;
      G4Exception("G4VUserPhysicsList::PreparePhysicsTable","No process manager",
                    FatalException,  particle->GetParticleName() );
      return;
    }

    G4ProcessVector* pVector = pManager->GetProcessList();
    if (!pVector) {
      G4cerr << "G4VUserPhysicsList::PreparePhysicsTable  : No Process Vector for " 
             << particle->GetParticleName() <<G4endl;
      G4cerr << particle->GetParticleName() << " should be created in your PhysicsList" <<G4endl;
      G4Exception("G4VUserPhysicsList::PreparePhysicsTable","No process Vector",
                    FatalException,  particle->GetParticleName() );
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
void G4VUserPhysicsList::DumpCutValuesTable(G4int nParticles) 
{ 
  fDisplayThreshold = nParticles; 
}

///////////////////////////////////////////////////////////////
void G4VUserPhysicsList::DumpCutValuesTableIfRequested()
{
  if(fDisplayThreshold==0) return;
  G4ProductionCutsTable::GetProductionCutsTable()->DumpCouples();
  fDisplayThreshold = 0;
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
  if (!fCutsTable->StoreCutsTable(dir, ascii)) {
    G4Exception("G4VUserPhysicsList::StorePhysicsTable","Faile to store ",
		JustWarning,"Cut Table can not be stored"); 	
    return false;
  }
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
	G4String comment =  "Fail to store for " + (*pVector)[j]->GetProcessName();
	comment += "(" + particle->GetParticleName()  + ")";
	G4Exception("G4VUserPhysicsList::StorePhysicsTable","Faile to store ",
		    JustWarning,comment); 	
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
	G4cout << "G4VUserPhysicsList::RetrievePhysicsTable   ";
	G4cout << " Fail to retrieve Physics Table for ";
        G4cout << (*pVector)[j]->GetProcessName() << G4endl;
	G4cout << "Calculate Physics Table for " << particle->GetParticleName() << G4endl;
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

void G4VUserPhysicsList::ResetCuts()
{
#ifdef G4VERBOSE  
  if (verboseLevel>0){
    G4cout << "G4VUserPhysicsList::ResetCuts() is obsolete.";
    G4cout << " This method gives no effect and you can remove it. "<< G4endl;
  }
#endif
}

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

G4bool G4VUserPhysicsList::GetApplyCuts(const G4String& name) const
{
  return theParticleTable->FindParticle(name)->GetApplyCutsFlag();
}




/// obsolete methods


void G4VUserPhysicsList::DumpCutValues(const G4String &particle_name)
{
  G4cerr << "WARNING !" << G4endl;
  G4cerr << " Obsolete DumpCutValues() method is invoked for " << particle_name << G4endl;
  G4cerr << " Please use DumpCutValuesTable() instead." << G4endl;
  G4cerr << " This dummy method implementation will be removed soon." << G4endl;
  DumpCutValuesTable();
}

void G4VUserPhysicsList::DumpCutValues(G4ParticleDefinition* )
{
  G4cerr << "WARNING !" << G4endl;
  G4cerr << " DumpCutValues() became obsolete." << G4endl;
  G4cerr << " Please use DumpCutValuesTable() instead." << G4endl;
  G4cerr << " This dummy method implementation will be removed soon." << G4endl;
  DumpCutValuesTable();
}
