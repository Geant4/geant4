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
// $Id: G4PhysicsTableHelper.cc 70369 2013-05-29 14:59:24Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// Class Description
//  G4PhysicsTableHelper is a static utility class 
//  for helping proceeses to build their physics table
//
// ------------------------------------------------------------
//   First Implementation          20 Aug. 2004   H.Kurashige
//
// ------------------------------------------------------------

#include "G4PhysicsTableHelper.hh" 
#include  "G4ProductionCutsTable.hh"

G4ThreadLocal G4int G4PhysicsTableHelper::verboseLevel = 1; 

G4PhysicsTableHelper::G4PhysicsTableHelper()
{
}

G4PhysicsTableHelper::~G4PhysicsTableHelper()
{
}

G4PhysicsTableHelper::G4PhysicsTableHelper(const G4PhysicsTableHelper&)
{
}

G4PhysicsTableHelper& G4PhysicsTableHelper::operator=(const G4PhysicsTableHelper&)
{
  return *this;
}


G4PhysicsTable* G4PhysicsTableHelper::PreparePhysicsTable(G4PhysicsTable* physTable)
{
  G4ProductionCutsTable* cutTable = G4ProductionCutsTable::GetProductionCutsTable();  
  size_t numberOfMCC = cutTable->GetTableSize(); 

  if ( physTable !=0) {
    // compare size of physics table and number of material-cuts-couple
    if ( physTable->size() < numberOfMCC) {
      // enlarge physcis table
      physTable->resize(numberOfMCC, (G4PhysicsVector*)(0));
#ifdef G4VERBOSE  
      if (verboseLevel>2) {
	G4cerr << "G4PhysicsTableHelper::PreparePhysicsTable  ";
	G4cerr << "Physics Table "<< physTable ;
	G4cerr << " is resized to " << numberOfMCC << G4endl;
      }
#endif 
    } else if ( physTable->size() > numberOfMCC){
      // ERROR: this situation should not occur  
      //  size of physics table is shorter than  number of material-cuts-couple
      physTable->resize(numberOfMCC);
#ifdef G4VERBOSE  
      if (verboseLevel>0) {
	G4cerr << "G4PhysicsTableHelper::PreparePhysicsTable  ";
	G4cerr << "Physics Table "<< physTable ;
	G4cerr << " is longer than number of material-cuts-couple " << G4endl;
      }
#endif 
      G4Exception( "G4PhysicsTableHelper::PreparePhysicsTable()",
		   "ProcCuts001", FatalException, 
		   "Physics Table is inconsistent with  material-cuts-couple");
    } 
  } else {
    // create PhysicsTable is given poitner is null
    physTable = new G4PhysicsTable(numberOfMCC);
    if (physTable!=0) {
      physTable->resize(numberOfMCC, (G4PhysicsVector*)(0));
    } else {
      G4Exception( "G4PhysicsTableHelper::PreparePhysicsTable()",
		   "ProcCuts002", FatalException, 
		   "Can't create Physics Table");
    }
  }

#ifdef G4VERBOSE  
  if (verboseLevel>2) {
    if ( physTable !=0) { 
      G4cerr << "Physics Table size "<< physTable->size();
    } else {
      G4cerr << "Physics Table does not exist   ";
    }
    G4cerr << ": number of material-cuts-couple " << numberOfMCC << G4endl;
  }
#endif 

  // Reset recal-needed flag for all physics vectors
  physTable->ResetFlagArray();

  for (size_t idx = 0; idx <numberOfMCC; idx +=1){
    const G4MaterialCutsCouple* mcc = cutTable->GetMaterialCutsCouple(idx);
    //check if re-calculation of the physics vector is needed 
    // MCC is not used
    if ( !mcc->IsUsed() ) physTable->ClearFlag(idx);

    // RecalcNeeded flag of MCC is not asserted 
    if ( !mcc->IsRecalcNeeded() ) physTable->ClearFlag(idx);
  }
  
  return physTable;
}



G4bool G4PhysicsTableHelper::RetrievePhysicsTable(G4PhysicsTable* physTable,
						  const G4String& fileName,
						  G4bool ascii              )
{
  if (physTable == 0) return false;
  
  // retrieve physics table from the given file
  G4PhysicsTable* tempTable = new G4PhysicsTable();
  if (! tempTable->RetrievePhysicsTable(fileName,ascii) ){
#ifdef G4VERBOSE  
    if (verboseLevel>1) {
      G4cerr << "G4PhysicsTableHelper::RetrievePhysicsTable  ";
      G4cerr << "Fail to retreive from "<< fileName << G4endl;
    }
#endif 
    G4Exception( "G4ProductionCutsTable::RetrievePhysicsTable()",
		 "ProcCuts105",
		 JustWarning, "Can not retrieve physics tables from file");
    delete tempTable;
    return false;
  } 

  G4ProductionCutsTable* cutTable = G4ProductionCutsTable::GetProductionCutsTable();  
  const G4MCCIndexConversionTable* converter = cutTable->GetMCCIndexConversionTable();

  // check physics table size
  if ( tempTable->size() != converter->size()){
#ifdef G4VERBOSE  
    if (verboseLevel>0) {
      G4cerr << "G4PhysicsTableHelper::RetrievePhysicsTable  ";
      G4cerr << "Size of the physics table in "<< fileName;
      G4cerr << "( size =" << tempTable->size() << ")";
      G4cerr << " is inconsistent with material-cut info";
      G4cerr << "( size =" << converter->size() << ")";
      G4cerr << G4endl;
    }
#endif
    G4Exception( "G4ProductionCutsTable::RetrievePhysicsTable()",
		 "ProcCuts106",
		 JustWarning, "Retrived file is inconsistent with current physics tables ");
    delete tempTable;
    return false;
  }
  
  // fill the given physics table with retrived physics vectors 
  for (size_t idx=0; idx<converter->size(); idx++){
    if (converter->IsUsed(idx)){
      if (converter->GetIndex(idx)<0) continue;
      size_t i = converter->GetIndex(idx);
      G4PhysicsVector* vec = (*physTable)[i];
       if (vec !=0 ) delete vec;
      (*physTable)[i] =  (*tempTable)[idx];
      physTable->ClearFlag(i);
    }
  }
  tempTable->clear();
  delete tempTable;

  return true;
}


void G4PhysicsTableHelper::SetPhysicsVector(G4PhysicsTable* physTable,
					    size_t idx,
					    G4PhysicsVector* vec)
{
  if ( physTable ==0) {  return;  }

  if ( physTable->size() <= idx) {
#ifdef G4VERBOSE  
    if (verboseLevel>0) {
      G4cerr << "G4PhysicsTableHelper::SetPhysicsVector   ";
      G4cerr << "Given index (" << idx << ")  exceeds ";
      G4cerr << "size of the physics table ";
      G4cerr << "( size =" << physTable->size()<< ")";
      G4cerr << G4endl;
    }
#endif
    G4Exception( "G4ProductionCutsTable::SetPhysicsVector()",
		 "ProcCuts107",
		 JustWarning, "Illegal index ");
    return;
  } 

  // set physics vector 
  (*physTable)[idx] = vec;
  // clear flag
  physTable->ClearFlag(idx);
 

}


void  G4PhysicsTableHelper::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

G4int G4PhysicsTableHelper::GetVerboseLevel()
{
  return verboseLevel;
}



