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
// $Id: G4PhysicsListHelper.cc,v 1.72 2010-07-21 14:21:19 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
//	History
//       first version                   29 Apr 2011 by H.Kurashige
// ------------------------------------------------------------

#include "globals.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"

#include "G4ios.hh"
#include <iomanip>
#include <fstream>

////////////////////////////////////////////////////////
G4PhysicsListHelper* G4PhysicsListHelper::pPLHelper = 0;
 
////////////////////////////////////////////////////////
G4PhysicsListHelper::G4PhysicsListHelper()
  :  useCoupledTransportation(false),
     theTransportationProcess(0),
     verboseLevel(1),
     theTable(0),
     sizeOfTable(0),
     ordParamFileName("")
{
  // pointer to the particle table
  theParticleTable = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();

  ReadOrdingParameterTable();

#ifdef G4VERBOSE
  if (verboseLevel >1){
    DumpOrdingParameterTable();
  }
#endif
}

////////////////////////////////////////////////////////
G4PhysicsListHelper::~G4PhysicsListHelper()
{
  if (theTable !=0) {
    theTable->clear();
    delete theTable;
    theTable=0;
    sizeOfTable=0;
  }
  if (theTransportationProcess!=0) {
    delete theTransportationProcess;
    theTransportationProcess=0;
  }
}

////////////////////////////////////////////////////////
G4PhysicsListHelper* G4PhysicsListHelper::GetPhysicsListHelper()  
{
  static G4PhysicsListHelper thePLHelper;
  if (pPLHelper == 0) pPLHelper = &thePLHelper;
  return pPLHelper;
}

////////////////////////////////////////////////////////
void G4PhysicsListHelper::CheckParticleList() const
{
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
      G4cout << "G4PhysicsListHelper::CheckParticleList: "
	     << missingName << " do not exist " << G4endl; 
      G4cout << " These particle are necessary for basic EM processes" 
	     << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::CheckParticleList",
		"RUN003", FatalException,
		"Missing EM basic particle");
  }

  // RULE 2
  //  proton should exist 
  //   if any other charged baryon  exist
  if (!isProton && isAnyChargedBaryon) {
    G4String missingName="proton ";

#ifdef G4VERBOSE
    if (verboseLevel >0){
      G4cout << "G4PhysicsListHelper::CheckParticleList: "
	     << missingName << " does not exist "<< G4endl; 
      G4cout << " Proton is necessary for EM baryon processes" << G4endl;
    }
#endif
    missingName += " should be created ";
    G4Exception("G4PhysicsListHelper::CheckParticleList",
		"RUN003", FatalException,
		"Missing Proton");
  }
   
  // RULE 3
  //  GenericIonn should exist 
  //   if any other ion  exist
  if (!isGenericIon && isAnyIon) {
    G4String missingName="GenericIon ";

#ifdef G4VERBOSE
    if (verboseLevel >0){
      G4cout << "G4PhysicsListHelper::CheckParticleList: "
	     << missingName << " does not exist "<< G4endl; 
      G4cout << " GenericIon should be created if any ion is necessary" << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::CheckParticleList",
		"RUN003", FatalException,
		"Missing GenericIon");
  }
      
}


////////////////////////////////////////////////////////
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"
#include "G4RunManagerKernel.hh"
#include "G4ScoringManager.hh"

void G4PhysicsListHelper::AddTransportation()
{
  G4int verboseLevelTransport = 0;

#ifdef G4VERBOSE
  if (verboseLevel >2){
    G4cout << "G4PhysicsListHelper::AddTransportation()  "<< G4endl;
  }
#endif

  G4int nParaWorld = 
    G4RunManagerKernel::GetRunManagerKernel()->GetNumberOfParallelWorld();
  
  if ( nParaWorld>0 || 
       useCoupledTransportation || 
       G4ScoringManager::GetScoringManagerIfExist()) {
#ifdef G4VERBOSE
    if (verboseLevel >0) {
      G4cout << " G4PhysicsListHelper::AddTransportation()"
	     << "--- G4CoupledTransportation is used " 
	     << G4endl;
    }
#endif
    theTransportationProcess = new G4CoupledTransportation(verboseLevelTransport);    
  } else {
    theTransportationProcess = new G4Transportation(verboseLevelTransport);
  }
 
  // loop over all particles in G4ParticleTable
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    // Add transportation process for all particles 
    if ( pmanager == 0) {
      // Error !! no process manager
#ifdef G4VERBOSE
      if (verboseLevel>0){
	G4cout << "G4PhysicsListHelper::AddTransportation  "
	       <<" : No Process Manager for "
	       << particle->GetParticleName() << G4endl;
      }
#endif
      G4Exception("G4PhysicsListHelper::AddTransportation",
		  "RUN001", FatalException,
		  "No process manager");
      continue;
    } 
    // add transportation with ordering = ( -1, "first", "first" )
    pmanager ->AddProcess(theTransportationProcess);
    pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
    pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
  }
}

////////////////////////////////////////////////////////
#include "G4ProcessManager.hh"
  
void G4PhysicsListHelper::ReadOrdingParameterTable()
{
  if( getenv("G4ORDPARAMTABLE") ){
    ordParamFileName = getenv("G4ORDPARAMTABLE");
#ifdef G4VERBOSE
    if (verboseLevel >1){
      G4cout << "G4PhysicsListHelper::ReadOrdingParameterTable  :"
	     << ordParamFileName << " is assigned to Ordering Parameter Table "
	     << G4endl; 
    }
#endif
  } else {
    ordParamFileName  = getenv("G4INSTALL");
    ordParamFileName += "/source/physics_lists/builders/";
    ordParamFileName += "OrderingParameterTable";
  } 
 
  std::ifstream fIn;  
  // open input file //
  fIn.open(ordParamFileName, std::ios::in);
  // check if the file has been opened successfully 
  if (!fIn) {
#ifdef G4VERBOSE
    if (verboseLevel >0) {
      G4cout << "G4PhysicsListHelper::ReadOrdingParameterTable  "
	     << " Can not open file " << ordParamFileName << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::ReadOrdingParameterTable",
		"RUN103", JustWarning, 
		"Fail to open ordering paramter table ");
    return;
  }

  // create OrdParamTable   
  if (theTable !=0) {
    theTable->clear();
    delete theTable;
    theTable=0;
    sizeOfTable=0;
  }
  theTable = new G4OrdParamTable();
  sizeOfTable=0;
  // read in the file and fill the table 
  while(!fIn.eof()) {
    G4PhysicsListOrderingParameter tmp;
    G4int flag;
    fIn >> tmp.processTypeName >>  tmp.processType >> tmp.processSubType
	>> tmp.ordering[0] >> tmp.ordering[1] >> tmp.ordering[2] >> flag;
    tmp.isDuplicable = (flag!=0);
    theTable->push_back(tmp);
    sizeOfTable +=1;  
  }
  fIn.close();
  
  if (sizeOfTable==0){
#ifdef G4VERBOSE
    if (verboseLevel >0) {
      G4cout << "G4PhysicsListHelper::ReadOrdingParameterTable "
	     << " Empty file " << ordParamFileName << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::ReadOrdingParameterTable",
		"RUN104", JustWarning, 
		"The ordering parameter table is empty ");
    delete theTable;
    theTable=0;
    sizeOfTable=0;
  }
  return;  
}

////////////////////////////////////////////////////////
void G4PhysicsListHelper::DumpOrdingParameterTable(G4int subType) const
{
 if (theTable==0) {
#ifdef G4VERBOSE
    if (verboseLevel >0) {
      G4cout << "G4PhysicsListHelper::DumpOrdingParameterTable   "
	     << " No ordering parameter table  : " << ordParamFileName 
	     << G4endl;
    }
#endif
    return;
  }
  G4cout << "G4PhysicsListHelper::DumpOrdingParameterTable  : "
	 << ordParamFileName << G4endl; 
  G4cout << "          TypeName  " 
	 << "    ProcessType" <<  "        SubType"
	 << "         AtRest" <<  "      AlongStep" <<  "        PostStep"
	 << "     Duplicable" << G4endl;
  for (int i=0; i<sizeOfTable ; i++){
    G4PhysicsListOrderingParameter* tmp=&(theTable->at(i));
    if ((subType>=0) && (subType!=tmp->processSubType)) continue;
    G4cout << std::setw(18)     << tmp->processTypeName 
	   << std::setw(15)     << tmp->processType 
	   << std::setw(15)	<< tmp->processSubType
	   << std::setw(15)	<< tmp->ordering[0] 
	   << std::setw(15)	<< tmp->ordering[1] 
	   << std::setw(15)	<< tmp->ordering[2];
    if (tmp->isDuplicable) {
      G4cout << "  true";
    } else {
      G4cout << "  false";
    }
    G4cout <<G4endl;
  }  
}

////////////////////////////////////////////////////////
G4PhysicsListOrderingParameter G4PhysicsListHelper::GetOrdingParameter(G4int subType) const
{
  G4PhysicsListOrderingParameter value;

 if (theTable==0) {
#ifdef G4VERBOSE
    if (verboseLevel >0) {
      G4cout << "G4PhysicsListHelper::GetOrderingParameter : "
	     << " No ordering parameter table  : " << ordParamFileName 
	     << G4endl;
    }
#endif
    return value;
  }

  for (int i=0; i<sizeOfTable ; i++){
    G4PhysicsListOrderingParameter* tmp=&(theTable->at(i));
    if (subType == tmp->processSubType){
      value.processTypeName = tmp->processTypeName; 
      value.processType     = tmp->processType; 
      value.processSubType  = tmp->processSubType; 
      value.ordering[0]     = tmp->ordering[0]; 
      value.ordering[1]     = tmp->ordering[1]; 
      value.ordering[2]     = tmp->ordering[2]; 
      value.isDuplicable    = tmp->isDuplicable;
    }
  }  
  return value;
}

////////////////////////////////////////////////////////
G4bool G4PhysicsListHelper::RegisterProcess(G4VProcess*            process,
					    G4ParticleDefinition*  particle)
{
  if (theTable==0) {
#ifdef G4VERBOSE
    if (verboseLevel >0) {
      G4cout << "G4PhysicsListHelper::RegisterProcess :"
	     << " No ordering parameter table  : " << ordParamFileName 
	     << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::RegisterPorcess",
		"RUN004", FatalException, 
		"No Ordering Parameter Table");
    return false;
  }

  const G4String pName = process->GetProcessName(); 
  const G4int pType    = process->GetProcessType();
  const G4int pSubType = process->GetProcessSubType();
  
#ifdef G4VERBOSE
  if (verboseLevel >2) {
    G4cout << "G4PhysicsListHelper::RegisterProcess :"
	   << pName << " Process Type = " << pType 
	   << " SubType = "<< pSubType  
	   << " to " << particle->GetParticleName()
	   << G4endl;
  }
#endif

  // Check Process Type/SubType
  if ((pType <1)||(pSubType<1)) {
#ifdef G4VERBOSE
    if (verboseLevel >0) {
      G4cout << "G4PhysicsListHelper::RegisterProcess :"
	     << pName << " for " << particle->GetParticleName()
	     << " has illegal Process Type = " << pType 
	     << " SubType = "<< pSubType << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::RegisterPorcess",
		"RUN003", FatalException, 
		"No Matching process Type/SubType");
    return false;
  }
  
  G4bool isFound = false;
  G4int  ord[3];
  G4bool duplicable = false;
  for (int i=0; i<sizeOfTable ; i++){
    G4PhysicsListOrderingParameter* tmp=&(theTable->at(i));
    if ((tmp->processType==pType)&&(tmp->processSubType==pSubType)){
      ord[0] = tmp->ordering[0]; 
      ord[1] = tmp->ordering[1]; 
      ord[2] = tmp->ordering[2];
      duplicable = tmp->isDuplicable;
      isFound = true;
      break;
    }
  } 
  if (!isFound) {
#ifdef G4VERBOSE
    if (verboseLevel >0) {
      G4cout << "G4PhysicsListHelper::RegisterProcess :"
	     << pName << " for " << particle->GetParticleName()
	     << " with  type/subtype =" 
	     << pType << "/" << pSubType 
	     << "  is not reigstered in OrdingParameterTable  "
	     << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::RegisterPorcess",
		"RUN003", FatalException, 
		"No Matching process Type/SubType");
    return false;
  }

  // Check Process Manager
  G4ProcessManager* pManager = particle->GetProcessManager();
  if ( pManager == 0) {
      // Error !! no process manager
#ifdef G4VERBOSE
    if (verboseLevel>0){
      G4cout << "G4PhysicsListHelper::RegisterProcess "
	     <<" : No Process Manager for "
	     << particle->GetParticleName() << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::RegisterProcess   ",
		"RUN001", FatalException,
		"No process manager");
    return false;
  }

  // Check Duplication
  if (!duplicable){
    G4bool duplicated = false;
    G4ProcessVector* pList = pManager->GetProcessList();
    for (G4int idx=0; idx<pList->size(); idx++) {
      const G4VProcess* p = (*pList)[idx];
      if ((p->GetProcessType()== pType)  && 
	  (p->GetProcessSubType()== pSubType)){
	duplicated = true;
#ifdef G4VERBOSE
	if (verboseLevel >0) {
	  G4cout << "G4PhysicsListHelper::RegisterProcess :"
		 << pName << " for " << particle->GetParticleName()
		 << " with  type/subtype =" 
		 << pType << "/" << pSubType 
		 << "  is has same subType as "
		 << p->GetProcessName()
		 << " for " << particle->GetParticleName()
		 << G4endl;
	  G4cout << "It will not be added !!" << G4endl;
	}
#endif
	G4Exception("G4PhysicsListHelper::RegisterPorcess",
		    "RUN105", JustWarning, 
		    "Duplication of processes");
      }
    }
    if (duplicated) return false;
  }

  // Add Process
  G4int code = pManager ->AddProcess(process);
  if (code <0) return false;

  // Set Ordering Parameter
  for(G4int idx=0; idx<3; idx++){
    G4ProcessVectorDoItIndex idxOrd = static_cast<G4ProcessVectorDoItIndex>(idx);
    if (ord[idx]<0) {
      // Do Nothing because NO DOIT
    } else if (ord[idx]==0) {
      pManager->SetProcessOrderingToFirst( process, idxOrd );
    } else if (ord[idx]<9999) {
      pManager->SetProcessOrdering( process, idxOrd , ord[idx]);
    } else {
      pManager->SetProcessOrderingToLast( process, idxOrd );
    } 
  } 
#ifdef G4VERBOSE
  if (verboseLevel >1) {
    G4cout << "G4PhysicsListHelper::RegisterProcess :"
	   << pName << " for " << particle->GetParticleName()
	   << " with  type/subtype =" 
	   << pType << "/" << pSubType 
	   << " is sucessfully registered with ordering parameters "
	   << ord[0] << ":" << ord[1] << ":" << ord[2]
	   << G4endl;
	}
#endif
  return true;
}
