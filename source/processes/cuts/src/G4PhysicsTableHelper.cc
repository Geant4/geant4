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
// G4PhysicsTableHelper class implementation 
//
// Author: H.Kurashige, 20 August 2004 - First implementation
// --------------------------------------------------------------------

#include "G4PhysicsTableHelper.hh" 
#include "G4ProductionCutsTable.hh"
#include "G4MCCIndexConversionTable.hh"
#include "G4Threading.hh"
#include "G4ios.hh"

G4int G4PhysicsTableHelper::verboseLevel = 1; 

// --------------------------------------------------------------------
G4PhysicsTableHelper::G4PhysicsTableHelper()
{
}

// --------------------------------------------------------------------
G4PhysicsTableHelper::~G4PhysicsTableHelper()
{
}

// --------------------------------------------------------------------
G4PhysicsTable*
G4PhysicsTableHelper::PreparePhysicsTable(G4PhysicsTable* physTable)
{
  G4ProductionCutsTable* cutTable
    = G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numberOfMCC = cutTable->GetTableSize(); 

  if ( physTable != nullptr )
  {
    // compare size of physics table and number of material-cuts-couple
    if ( physTable->size() < numberOfMCC )
    {
#ifdef G4VERBOSE  
      if (verboseLevel>2)
      {
        G4cout << "G4PhysicsTableHelper::PreparePhysicsTable: "
               << " the table " << physTable << " size="
	       << physTable->size()
               << " will be is resized to " << numberOfMCC << G4endl;
      }
#endif 
      // enlarge physics table
      physTable->resize(numberOfMCC, nullptr);
    }
    else if ( physTable->size() > numberOfMCC )
    {
      // ERROR: this situation should not occur  
      // size of physics table is larger than number of material-cuts-couple
      G4ExceptionDescription ed;
      ed << "table " << physTable << " size=" << physTable->size()
	 << " is longer than number of material-cuts-couple " << numberOfMCC; 
      G4Exception( "G4PhysicsTableHelper::PreparePhysicsTable()",
                   "ProcCuts001", FatalException, ed);
    } 
  }
  else
  {
    // create PhysicsTable is given poitner is null
    physTable = new G4PhysicsTable();
    physTable->resize(numberOfMCC, nullptr);
  }

#ifdef G4VERBOSE  
  if (verboseLevel>2)
  {
    G4cout << "G4PhysicsTableHelper::PreparePhysicsTable: "
	   << " the table "<< physTable
	   << " size=" << numberOfMCC << G4endl;
  }
#endif 

  // Reset recal-needed flag for all physics vectors
  physTable->ResetFlagArray();

  for (std::size_t idx = 0; idx <numberOfMCC; ++idx)
  {
    const G4MaterialCutsCouple* mcc = cutTable->GetMaterialCutsCouple((G4int)idx);

    // check if re-calculation of the physics vector is needed 
    // MCC is not used
    if ( !mcc->IsUsed() ) physTable->ClearFlag(idx);

    // RecalcNeeded flag of MCC is not asserted 
    if ( !mcc->IsRecalcNeeded() ) physTable->ClearFlag(idx);
  }
  
  return physTable;
}

// --------------------------------------------------------------------
G4bool G4PhysicsTableHelper::RetrievePhysicsTable(G4PhysicsTable* physTable,
                                                  const G4String& fileName,
                                                  G4bool ascii, G4bool spline)
{
  if (physTable == nullptr ) return false;
  
  // retrieve physics table from the given file
  G4PhysicsTable* tempTable = new G4PhysicsTable();
  if (! tempTable->RetrievePhysicsTable(fileName,ascii,spline) )
  {
    G4ExceptionDescription ed;
    ed << "Cannot retrieve physics table from the file <" << fileName << ">";
    G4Exception( "G4ProductionCutsTable::RetrievePhysicsTable()",
                 "ProcCuts105", JustWarning, ed);
    delete tempTable;
    return false;
  } 

  G4ProductionCutsTable* cutTable
    = G4ProductionCutsTable::GetProductionCutsTable();  
  const G4MCCIndexConversionTable* converter
    = cutTable->GetMCCIndexConversionTable();

  // check physics table size
  if ( tempTable->size() != converter->size())
  {
    G4ExceptionDescription ed;
    ed << "Physics table in " << fileName
       << "\n   size=" << tempTable->size() << " "
       << " is inconsistent with material-cut-couple "
       << "size=" << converter->size() << " the table is not retrieved!";
    G4Exception("G4ProductionCutsTable::RetrievePhysicsTable()",
                "ProcCuts106", JustWarning, ed);
    delete tempTable;
    return false;
  }
  
  // fill the given physics table with retrieved physics vectors 
  for (std::size_t idx=0; idx<converter->size(); ++idx)
  {
    if (converter->IsUsed(idx))
    {
      G4int i = converter->GetIndex(idx);
      if(i < 0) 
      {
        tempTable->clearAndDestroy();
	delete tempTable;
	return false;
      }	
      G4PhysicsVector* vec = (*physTable)[i];
      if (vec != nullptr ) delete vec;
      (*physTable)[i] = (*tempTable)[idx];
      physTable->ClearFlag(i);
    }
  }
  tempTable->clear();
  delete tempTable;

  return true;
}

// --------------------------------------------------------------------
void G4PhysicsTableHelper::SetPhysicsVector(G4PhysicsTable* physTable,
                                            std::size_t idx,
                                            G4PhysicsVector* vec)
{
  if ( physTable == nullptr) {  return;  }

  if ( physTable->size() <= idx)
  {
    G4ExceptionDescription ed;
    ed << "Given index (" << idx << ")  exceeds "
       << "the size of the physics table "
       << "( size =" << physTable->size() << ") the vector is not added!";
    G4Exception("G4ProductionCutsTable::SetPhysicsVector()",
                "ProcCuts107",
                JustWarning, ed);
    return;
  } 

  // set physics vector 
  (*physTable)[idx] = vec;
  // clear flag
  physTable->ClearFlag(idx);
}

// --------------------------------------------------------------------
void G4PhysicsTableHelper::SetVerboseLevel(G4int value)
{
  if( !G4Threading::IsWorkerThread() ) verboseLevel = value;
}

// --------------------------------------------------------------------
G4int G4PhysicsTableHelper::GetVerboseLevel()
{
  return verboseLevel;
}
