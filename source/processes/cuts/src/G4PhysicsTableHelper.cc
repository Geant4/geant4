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
// $Id: G4PhysicsTableHelper.cc,v 1.3 2004/12/15 15:42:37 gunter Exp $
// GEANT4 tag $Name: geant4-07-01 $
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

G4int G4PhysicsTableHelper::verboseLevel = 1; 

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
    } 
  } else {
    // create PhysicsTable is given poitner is null
    physTable = new G4PhysicsTable(numberOfMCC);
    physTable->resize(numberOfMCC, (G4PhysicsVector*)(0));

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
    return false;
  }
  
  // fill the given physics table with retrived physics vectors 
  for (size_t idx=0; idx<converter->size(); idx++){
    if (converter->IsUsed(idx)){
      size_t i = converter->GetIndex(idx);
      G4PhysicsVector* vec = (*physTable)[i];
      if (vec !=0 ) delete vec;
      (*physTable)[i] =  (*tempTable)[idx];
      physTable->ClearFlag(i);
    }
  }

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
    return;
  } 

  // set physics vector 
  (*physTable)[idx] = vec;
  // clear flag
  physTable->ClearFlag(idx);
 

}




