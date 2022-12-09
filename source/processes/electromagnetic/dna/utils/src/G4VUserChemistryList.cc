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
/*
 * G4VUserChemistryList.cc
 *
 *  Created on: 23 oct. 2013
 *      Author: kara
 */

#include "G4VUserChemistryList.hh"

#include <G4VScheduler.hh>
#include "G4MoleculeTable.hh"
#include "G4MoleculeDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4DNAChemistryManager.hh"

G4VUserChemistryList::G4VUserChemistryList(G4bool flag) :
fIsPhysicsConstructor(flag)
{
  verboseLevel = 1;
}

G4VUserChemistryList::~G4VUserChemistryList()
{
  G4DNAChemistryManager* chemMan = G4DNAChemistryManager::GetInstanceIfExists();
  if (chemMan != nullptr)
  {
    chemMan->Deregister(*this);
  }
}

void G4VUserChemistryList::RegisterTimeStepModel(G4VITStepModel* timeStepModel,
                                                 G4double startingTime)
{
  G4VScheduler::Instance()->RegisterModel(timeStepModel, startingTime);
}

void G4VUserChemistryList::BuildPhysicsTable()
{
  G4MoleculeTable* theMoleculeTable = G4MoleculeTable::Instance();

  G4MoleculeDefinitionIterator iterator =
      theMoleculeTable->GetDefintionIterator();

  iterator.reset();
  while (iterator())
  {
    G4MoleculeDefinition* moleculeDef = iterator.value();
    BuildPhysicsTable(moleculeDef);
  }
}

void G4VUserChemistryList::BuildPhysicsTable(G4MoleculeDefinition* moleculeDef)
{
  //Get processes from master thread;
  G4ProcessManager* pManager = moleculeDef->GetProcessManager();

  if (pManager == nullptr)
  {
#ifdef G4VERBOSE
    if (verboseLevel > 0)
    {
      G4cout << "G4VUserPhysicsList::BuildPhysicsTable "
             << " : No Process Manager for " << moleculeDef->GetParticleName()
             << G4endl;
      G4cout << moleculeDef->GetParticleName()
      << " should be created in your PhysicsList" <<G4endl;
    }
#endif
    G4Exception("G4VUserChemistryList::BuildPhysicsTable",
        "Run0271", FatalException,
        "No process manager");
    return;
  }

  G4ProcessManager* pManagerShadow = moleculeDef->GetMasterProcessManager();
  G4ProcessVector* pVector = pManager->GetProcessList();
  if (pVector == nullptr)
  {
#ifdef G4VERBOSE
    if (verboseLevel > 0)
    {
      G4cout << "G4VUserChemistryList::BuildPhysicsTable  "
             << " : No Process Vector for " << moleculeDef->GetParticleName()
             << G4endl;
    }
#endif
    G4Exception("G4VUserChemistryList::BuildPhysicsTable",
        "Run0272", FatalException,
        "No process Vector");
    return;
  }
#ifdef G4VERBOSE
  if (verboseLevel > 2)
  {
    G4cout << "G4VUserChemistryList::BuildPhysicsTable %%%%%% "
           << moleculeDef->GetParticleName() << G4endl;
    G4cout << " ProcessManager : " << pManager
           << " ProcessManagerShadow : " << pManagerShadow << G4endl;
    for(G4int iv1=0;iv1<(G4int)pVector->size();++iv1)
    {
      G4cout << "  " << iv1 << " - " << (*pVector)[iv1]->GetProcessName()
          << G4endl;
    }
    G4cout << "--------------------------------------------------------------"
        << G4endl;
    G4ProcessVector* pVectorShadow = pManagerShadow->GetProcessList();

    for(G4int iv2=0;iv2<(G4int)pVectorShadow->size();++iv2)
    {
      G4cout << "  " << iv2 << " - " << (*pVectorShadow)[iv2]->GetProcessName()
      << G4endl;
    }
  }
#endif
  for (G4int j = 0; j < (G4int)pVector->size(); ++j)
  {
    //Andrea July 16th 2013 : migration to new interface...
    //Infer if we are in a worker thread or master thread
    //Master thread is the one in which the process manager
    // and process manager shadow pointers are the same
    if (pManagerShadow == pManager)
    {
      (*pVector)[j]->BuildPhysicsTable(*moleculeDef);
    }
    else
    {
      (*pVector)[j]->BuildWorkerPhysicsTable(*moleculeDef);
    }

  }
}
